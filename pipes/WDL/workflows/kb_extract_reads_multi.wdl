version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow kb_extract_reads {
    meta {
        description: "Runs multiple FASTQ/BAM files through kb extract process"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        Array[File]+    reads_bams
        Array[File]+    h5ad_files

        Boolean         protein=false
        Int             threshold=1

        File            kallisto_index
        File            t2g    
    }

    parameter_meta {
        reads_bams: {
          description: "Reads to process. Unmapped, single-end or interleaved. Corresponding h5ad files are matched by basename (e.g., sample1.bam matches sample1.h5ad, sample1.fastq.gz matches sample1.h5ad).",
          patterns: ["*.bam", "*.fastq", "*.fastq.gz"]
        }
        h5ad_files: {
          description: "h5ad files containing hits for each read. Files are matched to BAMs by basename. Order does not matter.",
          patterns: ["*.h5ad"]
        }
        protein: {
          description: "Whether the extraction process is dealing with amino-acid sequences.",
          options: ["true", "false"]
        }
        threshold: {
          description: "Minimum number of reads that must map to a given transcript for it to be included in the output. (default = 1)"
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

    # Process each BAM file with its matched h5ad
    scatter(reads_bam in reads_bams) {
        # Extract basename from reads file
        String reads_basename = sub(sub(sub(sub(basename(reads_bam), "\\.gz$", ""), "\\.fastq$", ""), "\\.R[12]$", ""), "\\.[12]$", "")
        
        # Find matching h5ad file by comparing basenames
        scatter(h5ad in h5ad_files) {
            String h5ad_basename = sub(basename(h5ad), "\\.h5ad$", "")
            Boolean is_match = (reads_basename == h5ad_basename)
            File? matched_file = if is_match then h5ad else None
        }
        
        # Select the matching h5ad (will have exactly one non-None value)
        Array[File] matched = select_all(matched_file)
        File matched_h5ad = matched[0]
        
        call metagenomics.kb_extract as kb_extract_single {
            input: 
                reads_bam = reads_bam,
                h5ad_file = matched_h5ad,
                protein = protein,
                kb_index = kallisto_index,
                t2g = t2g,
                threshold = threshold
        }
    }

    ## TODO: Re-tag all of our sequence IDs from extract reads with our sample names + kallisto DB IDs
    ## TODO: Parse kraken2 read summaries to create a new file containing taxonomic ID's resolve + read IDs
    ## TODO: Run kraken2 - kallisto read classifier script to create a classification table
    

    output {
        Array[File] kb_extracted_reads = kb_extract_single.kb_extract_tar
    }
}
