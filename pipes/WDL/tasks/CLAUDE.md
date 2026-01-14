# CLAUDE.md - WDL Tasks

Navigation index for WDL task modules.

| File | WHAT | WHEN |
|------|------|------|
| tasks_16S_amplicon.wdl | 16S rRNA amplicon analysis tasks | Working with bacterial 16S sequencing data |
| tasks_assembly.wdl | De novo and reference-based assembly tasks | Assembling viral genomes from reads |
| tasks_demux.wdl | Sample demultiplexing utilities | Processing multiplexed sequencing runs |
| tasks_interhost.wdl | Phylogenetics and multi-sample analysis (includes BEAST GPU task) | Building phylogenetic trees, population genetics |
| tasks_intrahost.wdl | Variant calling and iSNV analysis | Detecting within-host viral variants |
| tasks_megablast.wdl | BLAST-based sequence search | Identifying sequences via homology |
| tasks_metagenomics.wdl | Taxonomic classification and filtering (Kraken2, Kaiju, kb, BLASTx, VirNucPro). Now includes classify_virnucpro task: GPU-accelerated viral classifier using DNABERT_S + ESM2-3B models | Classifying reads, metagenomic analysis, GPU-based viral sequence classification |
| tasks_ncbi.wdl | GenBank/SRA submission workflows | Submitting data to NCBI repositories |
| tasks_ncbi_tools.wdl | NCBI data retrieval utilities | Downloading sequences from GenBank/SRA |
| tasks_nextstrain.wdl | Nextstrain/Augur integration | Creating Nextstrain visualizations |
| tasks_read_utils.wdl | BAM/FASTQ processing utilities | Manipulating sequencing data files |
| tasks_reports.wdl | QC and analysis reports | Generating MultiQC and custom reports |
| tasks_sarscov2.wdl | SARS-CoV-2 specific tools | SARS-CoV-2 lineage calling, variant analysis |
| tasks_taxon_filter.wdl | Host/contaminant filtering | Removing non-viral reads |
| tasks_terra.wdl | Terra platform utilities | Terra-specific data management |
| tasks_utils.wdl | General utility tasks | Common operations across workflows |
| README.md | VirNucPro GPU task architecture and invariants | Understanding GPU task design decisions |
