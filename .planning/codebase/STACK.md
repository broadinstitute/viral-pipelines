# Technology Stack

**Analysis Date:** 2026-02-11

## Languages

**Primary:**
- WDL (Workflow Description Language) 1.0 - Bioinformatics workflow definitions in `pipes/WDL/workflows/` and `pipes/WDL/tasks/`
- Python 3.10+ - Documentation generation, CI/CD scripts, and utilities

**Secondary:**
- Bash/Shell - CI/CD automation scripts in `github_actions_ci/`
- Java 17 - Required for workflow execution engines (Cromwell, womtool)

## Runtime

**Environment:**
- Python 3.10 (specified in `.readthedocs.yml` and GitHub Actions workflows)
- Java 17 (Temurin distribution, specified in `.github/workflows/build.yml`)
- Docker - All computational tasks execute within Docker containers

**Package Manager:**
- Pixi - Modern conda package manager (configured in `pixi.toml`)
- Lockfile: `pixi.lock` present
- pip - Python package management for documentation dependencies

## Frameworks

**Core:**
- **WDL Execution Engines:**
  - miniWDL >=1.13.1,<2 - Local workflow execution and validation (specified in `pixi.toml`)
  - Cromwell 88 - Cloud and local workflow execution (specified in `pixi.toml` and installed via `github_actions_ci/install-wdl.sh`)
  - womtool 88 - WDL validation tool (installed via `github_actions_ci/install-wdl.sh`)

- **Cloud/Platform Integration:**
  - dxWDL v1.50 - DNAnexus platform integration (installed via `github_actions_ci/install-wdl.sh`)
  - dxCompiler 2.11.6 - DNAnexus compiler (installed via `github_actions_ci/install-wdl.sh`)
  - dx-toolkit v0.311.0 - DNAnexus CLI tools (installed via `github_actions_ci/install-wdl.sh`)

**Testing:**
- miniWDL - Workflow testing framework (tests in `test/input/WDL/miniwdl-local/`)
- Cromwell - Integration testing (tests executed via `github_actions_ci/tests-cromwell.sh`)

**Build/Dev:**
- Sphinx 7.4.7 - Documentation generation (specified in `docs/requirements.txt`)
- sphinx-rtd-theme >=2.0.0 - ReadTheDocs theme
- wdl-aid 1.0.0 - WDL documentation helper

## Key Dependencies

**Critical:**
- **viral-ngs Docker images** - Core bioinformatics toolkit
  - ghcr.io/broadinstitute/viral-ngs:3.0.4-* (multiple variants: assemble, classify, core, phylo, baseimage)
  - Specified in `requirements-modules.txt` and referenced throughout WDL tasks
  - Versions managed centrally and validated via `github_actions_ci/check-wdl-runtimes.sh`

**Infrastructure:**
- PyYAML 6.0.1 - Configuration parsing
- Jinja2 3.1.4 - Template rendering for documentation
- sphinx-argparse 0.5.2 - CLI documentation
- recommonmark - Markdown support in Sphinx

**Bioinformatics Tools (via Docker):**
- broadinstitute/read-qc-tools:1.0.1
- broadinstitute/py3-bio:0.1.3
- broadinstitute/beast-beagle-cuda:1.10.5pre
- broadinstitute/ncbi-tools:2.11.1
- nextstrain/base:build-20240318T173028Z
- nextstrain/nextclade:3.18.1
- andersenlabapps/ivar:1.3.1
- quay.io/staphb/pangolin:4.3.3-pdata-1.36
- quay.io/broadinstitute/viral-classify:2.1.33.0

All bioinformatics tool versions are pinned in `requirements-modules.txt`.

## Configuration

**Environment:**
- No `.env` files detected - configuration via WDL input JSON files
- Workflow inputs specified per-workflow in `test/input/WDL/` directory
- Platform-specific configurations in `pipes/cromwell/` and `pipes/dnax/`

**Build:**
- `pixi.toml` - Pixi environment configuration
- `.readthedocs.yml` - ReadTheDocs build configuration
- `.dockstore.yml` - Dockstore registry configuration (107 workflows registered)
- `docs/conf.py` - Sphinx documentation configuration
- `pipes/cromwell/*.conf` - Cromwell execution engine configurations
- `pipes/dnax/*.json` - DNAnexus platform configurations

## Platform Requirements

**Development:**
- Linux-64 platform (specified in `pixi.toml`)
- Docker runtime for workflow execution
- 8-32GB RAM (workflow-dependent, specified in WDL runtime blocks)
- Local SSD storage (sizes calculated dynamically in WDL tasks)

**Production:**
- Multi-platform deployment target:
  - **Terra/GCP** - Google Cloud Platform with PAPI v2 or GCP Batch backend (configuration in `pipes/cromwell/cromwell.gcid-viral-seq.conf`)
  - **DNAnexus** - Cloud platform deployment (configuration in `pipes/dnax/`)
  - **AWS/Azure** - Via Cromwell backends (supported via WDL runtime parameters)
  - **Local HPC** - Via miniWDL or Cromwell
- Docker registry access: ghcr.io, quay.io, docker.io
- Google Cloud Storage (GCS) for Terra workflows
- Terra workspace integration for data table operations

---

*Stack analysis: 2026-02-11*
