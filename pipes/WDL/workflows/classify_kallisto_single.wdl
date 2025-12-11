version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_kallisto_single {
    meta {
        description: "Run a single FASTQ/BAM file through the kb classify process"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File            reads_bam

        Int             kmer_size=31
        Int             threshold=1

        String          technology
        String          parity

        Boolean         h5ad=false
        Boolean         loom=false
        Boolean         protein=false

        File            kallisto_index
        File            id_to_taxa_map  
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
          description: "When extracting hit/gene ID's from an a5ad before extract, minimum read threshold to filter on. Default is 1"
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
        id_to_taxa_map: {
          description: "Mapping file from transcript IDs to taxonomic IDs.",
          patterns: ["*.csv", "*.csv.gz", "*.csv.zst"]
        }
        t2g: {
          description: "Transcript-to-gene mapping file.",
          patterns: ["*.tsv", "*.tsv.gz", "*.tsv.zst"]
        }
    }


    call metagenomics.kb as classify_kb_single {
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

    call metagenomics.kb_extract as kb_extract_single {
        input: 
            reads_bam = reads_bam,
            h5ad_file = classify_kb_single.kb_count_tar,
            protein = protein,
            kb_index = kallisto_index,
            t2g = t2g,
            threshold = threshold
    }

    call metagenomics.report_primary_kb_taxa {
        input:
            kb_count_tar = classify_kb_single.kb_count_tar,
            id_to_taxon_map = id_to_taxa_map
    }

    ## TODO: Re-tag all of our sequence IDs from extract reads with our sample names + kallisto DB IDs
    ## TODO: Parse kallisto read summaries to create a new file containing taxonomic ID's resolve + read IDs
    ## TODO: Run kallisto - kallisto read classifier script to create a classification table
    
    output {
        File  kb_classify_reads                 = classify_kb_single.kb_count_tar
        File  kb_extracted_reads                = kb_extract_single.kb_extract_tar
        File   kallisto_top_taxa_report         = report_primary_kb_taxa.ranked_focal_report
        String kallisto_focal_taxon_name        = report_primary_kb_taxa.focal_tax_name
        Int    kallisto_focal_total_reads       = report_primary_kb_taxa.total_focal_reads
        String kallisto_top_taxon_id            = report_primary_kb_taxa.tax_id
        String kallisto_top_taxon_name          = report_primary_kb_taxa.tax_name
        String kallisto_top_taxon_rank          = report_primary_kb_taxa.tax_rank
        Int    kallisto_top_taxon_num_reads     = report_primary_kb_taxa.num_reads
        Float  kallisto_top_taxon_pct_of_focal  = report_primary_kb_taxa.percent_of_focal
    }
}
