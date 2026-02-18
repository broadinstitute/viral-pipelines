version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_kallisto_multi {
    meta {
        description: "Runs multiple FASTQ/BAM files through kb classify process"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        Array[File]+    reads_bams

        Int             kmer_size=31
        Int             threshold=1

        String          technology
        String          parity

        Boolean         h5ad=false
        Boolean         loom=false
        Boolean         protein=false

        File            kallisto_index
        File            t2g    
    }

    parameter_meta {
        reads_bams: {
          description: "Set of files containing reads to process. Single-end or paired-end.",
          patterns: ["*.bam", "*.fastq", "*.fastq.gz"]
        }
        kmer_size: {
          description: "K-mer size to use for classification. Default is 31."
        }
        threshold: {
          description: "When extracting hit/gene ID's from an h5ad before extract, minimum read threshold to filter on. Default is 1"
        }
        technology: {
          description: "Technology used for sequencing (e.g., '10xv2', '10xv3')."
        }
        parity: {
          description: "Parity of the reads ('single' or 'paired')."
        }
        h5ad: {
          description: "Output matrix in HDF5 format. Default is false."
        }
        loom: {
          description: "Output a loom file. Default is false"
        }
        protein: {
          description: "Whether the extraction process is dealing with amino-acid sequences.",
          options: ["true", "false"]
        }
        kallisto_index: {
          description: "Kallisto index file for the reference transcriptome.",
          patterns: ["*.idx", "*.idx.gz", "*.idx.zst"]
        }
        t2g: {
          description: "Transcript-to-gene mapping file.",
          patterns: ["*.tsv", "*.tsv.gz", "*.tsv.zst"]
        }
    }

    # Process each BAM file and lookup corresponding h5ad
    scatter(reads_bam in reads_bams) {
        # Strip common read file extensions to get sample basename
        String bam_basename = basename(basename(basename(reads_bam, ".bam"), ".fastq.gz"), ".fastq")
        
        call metagenomics.kallisto as classify_kallisto_single {
            input: 
                reads_bam = reads_bam,
                kmer_size = kmer_size,
                technology = technology,
                parity = parity,
                h5ad = h5ad,
                loom = loom,
                protein = protein,
                kb_index = kallisto_index,
                t2g = t2g,
        }
    }

    ## Now let's merge back our h5ad files
    call metagenomics.kallisto_merge_h5ads as merge_h5ads {
        input:
            in_count_tars = classify_kallisto_single.kb_count_tar,
            out_basename = "merged_kallisto_classify",
    }

    # Now we need to go ahead and extract the reads that were classified by kallisto
    # We'll need the individual h5ad files for this.
    # Dual scatter over reads_bams and count_tars in parallel
    scatter(pair in zip(reads_bams, classify_kallisto_single.kb_count_tar)) {
        File in_bam = pair.left
        File count_tar = pair.right
        
        call metagenomics.kallisto_extract as kallisto_extract_single {
            input:
                reads_bam = in_bam,
                h5ad_file = count_tar,
                protein = protein,
                kb_index = kallisto_index,
                t2g = t2g,
                threshold = threshold
        }
    }

    ## TODO: Re-tag all of our sequence IDs from extract reads with our sample names + kallisto DB IDs
    
    output {
        Array[File] kb_classify_reads   = classify_kallisto_single.kb_count_tar
        Array[File] kb_extracted_reads  = kallisto_extract_single.kb_extract_tar
        File        merged_h5ad         = merge_h5ads.kb_merged_h5ad
    }
}
