# CLAUDE.md - WDL Workflows

Navigation index for WDL workflow definitions.

| File | WHAT | WHEN |
|------|------|------|
| align_and_count.wdl | Alignment with read counting | Quantifying reads mapping to reference |
| align_and_count_multiple_report.wdl | Multi-sample alignment reporting | Batch alignment statistics |
| align_and_plot.wdl | Alignment with coverage plots | Visualizing genome coverage |
| assemble_denovo.wdl | De novo assembly with SPAdes | Assembling without reference genome |
| assemble_denovo_metagenomic.wdl | Metagenomic de novo assembly | Assembling mixed viral samples |
| assemble_refbased.wdl | Reference-based consensus calling | Primary workflow for viral genome assembly |
| augur_export_only.wdl | Nextstrain export step only | Converting data for Nextstrain visualization |
| augur_from_assemblies.wdl | Nextstrain phylogenetic analysis from assemblies | Building Nextstrain trees from consensus sequences |
| augur_from_beast_mcc.wdl | Nextstrain from BEAST MCC tree | Converting BEAST output to Nextstrain |
| augur_from_mltree.wdl | Nextstrain from ML tree | Building Nextstrain from existing phylogeny |
| augur_from_msa.wdl | Nextstrain from multiple sequence alignment | Complete Nextstrain workflow from MSA |
| augur_from_msa_with_subsampler.wdl | Nextstrain with subsampling | Large dataset Nextstrain with downsampling |
| bam_to_qiime.wdl | BAM to QIIME format conversion | Preparing data for QIIME2 analysis |
| bams_multiqc.wdl | MultiQC report for BAM files | Aggregated QC metrics across samples |
| beast_gpu.wdl | BEAST phylogenetics with GPU | Bayesian phylogenetic inference on GPU |
| blastoff.wdl | BLAST-based classification | Homology-based sequence identification |
| calc_bam_read_depths.wdl | Read depth calculation | Computing coverage statistics |
| classify_kaiju.wdl | Kaiju taxonomic classification | Protein-based read classification |
| classify_kallisto_single.wdl | Kallisto single-sample quantification | RNA-seq abundance estimation |
| classify_kb.wdl | kb-python taxonomic classification | Fast RNA-seq pseudoalignment classification |
| classify_virnucpro_single.wdl | Single-sample VirNucPro workflow | Classifying one BAM file, testing VirNucPro integration |
| classify_virnucpro_multi.wdl | Multi-sample VirNucPro workflow with scatter (1 GPU instance per BAM) | Batch classifying multiple BAM files in parallel, production VirNucPro runs |
