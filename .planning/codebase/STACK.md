# Technology Stack

**Analysis Date:** 2026-03-31

## Languages

**Primary:**
- WDL (Workflow Description Language) - All workflow and task definitions in `pipes/WDL/`
- Python 3.10+ - Used in Python scripts within WDL command blocks and documentation generation
- Shell (Bash) - Command execution in WDL tasks and CI/CD scripts

**Secondary:**
- YAML - Configuration files (`.dockstore.yml`, `.readthedocs.yml`, `.github/workflows/build.yml`)
- JSON - Test input/output definitions and DNAnexus configuration

## Runtime

**Environment:**
- Python 3.14.2+ (as minimum, <3.15) - Defined in `pixi.toml`
- Docker - All tasks execute within containerized environments
- miniWDL - Local workflow execution engine
- Cromwell - Multi-platform workflow execution engine (local, HPC, cloud backends)

**Package Manager:**
- Pixi - For Python environment management (`pixi.toml`)
- pip - For Python dependencies in documentation build (`docs/requirements.txt`)
- Docker image registries: GHCR (GitHub Container Registry), Quay.io, Biocontainers

## Frameworks

**Core Workflow:**
- WDL 1.0 - Workflow Description Language for task and workflow definitions
- Cromwell - Workflow execution and orchestration engine
- miniWDL - Lightweight local WDL executor

**Testing:**
- GitHub Actions - CI/CD pipeline (`.github/workflows/build.yml`)
- miniWDL test runner - For local integration testing
- Cromwell test runner - For cloud platform testing

**Build/Dev:**
- Pixi - Environment management and dependency resolution
- ShellCheck - Shell script linting (installed in CI)
- miniwdl check - WDL syntax and semantic validation

**Documentation:**
- Sphinx - Documentation generation
- ReadTheDocs - Documentation hosting and automated builds
- wdl-aid - WDL documentation aid library

## Key Dependencies

**Critical Docker Images (from `requirements-modules.txt`):**
- `broadinstitute/viral-ngs:3.0.6` (multiple flavors) - Core bioinformatics toolkit
  - `3.0.6-core` - Base tools (samtools, bwa, minimap2, picard, etc.)
  - `3.0.6-assemble` - Assembly tools (SPAdes, reference-based assembly)
  - `3.0.6-classify` - Classification tools (Kraken2, Kaiju, Centrifuge)
  - `3.0.6-phylo` - Phylogenetics tools (Augur, MAFFT, RAxML)
  - `3.0.6-baseimage` - Base image for utility tasks

- `broadinstitute/read-qc-tools:1.0.1` - Read quality control utilities
- `broadinstitute/py3-bio:0.1.3` - Python 3 bioinformatics libraries (Biopython, pysam)
- `broadinstitute/beast-beagle-cuda:1.10.5pre` - Bayesian phylogenetics with GPU support
- `broadinstitute/ncbi-tools:2.11.1` - NCBI submission tools
- `broadinstitute/viral-classify:2.1.33.0` - Classification pipelines
- `broadinstitute/virnucpro-cuda:1.0.9` - VirNucPro viral nucleotide profiling (GPU)

**External Bioinformatics Tools (via Docker):**
- `nextstrain/nextclade:3.18.1` - Nextclade viral clade assignment
- `nextstrain/base:build-20240318T173028Z` - Nextstrain augur phylogenetics
- `quay.io/staphb/pangolin:4.3.3-pdata-1.36` - SARS-CoV-2 lineage assignment (Pangolin)
- `andersenlabapps/ivar:1.3.1` - iVar intrahost variant calling
- `quay.io/biocontainers/trimal:1.4.1` - Alignment trimming
- `quay.io/biocontainers/bcftools:1.10.2` - VCF tools
- `biocontainers/krona:v2.7.1_cv1` - Interactive metagenomic visualization
- `quay.io/biocontainers/multiqc:1.32` - Quality control aggregation
- `mirror.gcr.io/staphb/vadr:1.6.4` - Viral annotation and detection
- `quay.io/broadinstitute/gisaid-cli:3.0` - GISAID sequence submission
- `quay.io/broadinstitute/sc2-rmd:0.1.25` - SARS-CoV-2 report generation
- `quay.io/broadinstitute/qiime2:latest` - 16S amplicon analysis
- `quay.io/broadinstitute/polyphonia:latest` - Intrahost variant calling
- `quay.io/broadinstitute/reconstructr:main` - Phylogenetic reconstruction
- `quay.io/broadinstitute/subsampler` - Sequence subsampling
- `python:slim` - Lightweight Python image for utility tasks
- `ubuntu` - Base Ubuntu for generic tasks

**Documentation Dependencies (from `docs/requirements.txt`):**
- Sphinx==7.4.7
- sphinx-rtd-theme>=2.0.0
- sphinx-argparse==0.5.2
- jinja2==3.1.4
- PyYAML==6.0.1
- recommonmark
- wdl-aid==1.0.0

## Configuration

**Environment:**
- Pixi (`pixi.toml`) - Declares Python version and channels (conda-forge)
- Dockstore (`.dockstore.yml`) - Registers 100+ workflows for platform access
- ReadTheDocs (`.readthedocs.yml`) - Sphinx documentation build configuration

**Build:**
- GitHub Actions (`build.yml`) - CI jobs for validation, testing, documentation build, and deployment
- Cromwell config files (`pipes/cromwell/`) - Backend configurations for local and GCP execution
- DNAnexus defaults (`pipes/dnax/dx-defaults-*.json`) - Platform-specific default values
- Docker version pinning (`requirements-modules.txt`) - Exact Docker image versions for reproducibility

**Version Management:**
- `requirements-modules.txt` - Single source of truth for all Docker image versions
- Scripts to validate (`github_actions_ci/check-wdl-runtimes.sh`) and update (`github_actions_ci/version-wdl-runtimes.sh`) WDL runtime declarations against this file

## Platform Requirements

**Development:**
- Python 3.10+ (for documentation)
- Docker (for task execution)
- miniWDL or Cromwell (for workflow execution)
- Git (for version control)
- Shell environment (bash/zsh)
- Optional: Java 11+ (for Cromwell/womtool, not required for miniWDL)

**Production/Execution:**
- Docker runtime - All compute tasks containerized
- Can run on:
  - Local machines with Docker and miniWDL/Cromwell
  - GCP (Google Cloud Platform) - Via Cromwell backends
  - Terra (academic cloud platform) - Workflows deployable via Dockstore
  - DNAnexus (commercial platform) - Via dx-launcher applet
  - AWS - Via Cromwell CromwellOnAzure
  - Azure - Via Microsoft's CromwellOnAzure
  - HPC systems - Via Cromwell job scheduler backends

---

*Stack analysis: 2026-03-31*
