# Codebase Structure

**Analysis Date:** 2026-02-11

## Directory Layout

```
viral-pipelines/
├── pipes/                      # Core WDL pipeline definitions
│   ├── WDL/
│   │   ├── workflows/          # End-to-end analysis workflows (90+ files)
│   │   └── tasks/              # Reusable task modules (16 files)
│   ├── cromwell/               # Cromwell-specific configurations
│   └── dnax/                   # DNAnexus deployment utilities
├── test/                       # Test data and configuration
│   └── input/                  # Test inputs organized by workflow
│       ├── WDL/                # JSON test input files
│       │   └── miniwdl-local/  # miniWDL-specific test configs (symlinks)
│       ├── genbank/            # Sample GenBank files
│       ├── TestSplitcodeDemuxFastqs/  # Demux test data
│       └── zika-tutorial/      # Tutorial data
├── github_actions_ci/          # CI/CD validation scripts
├── docs/                       # Sphinx documentation source
├── .github/                    # GitHub Actions workflows
├── .planning/                  # GSD planning documents
│   └── codebase/               # Auto-generated codebase analysis
├── requirements-modules.txt    # Docker image version pins
├── .dockstore.yml              # Dockstore registry configuration
├── AGENTS.md                   # AI agent guidance and conventions
├── CLAUDE.md                   # Entry point for Claude Code
└── README.md                   # Project overview
```

## Directory Purposes

**pipes/WDL/workflows/:**
- Purpose: End-to-end analysis pipelines for viral genomics
- Contains: 90+ WDL workflow files
- Key files:
  - `assemble_refbased.wdl` (250 lines) - Primary reference-based assembly
  - `classify_single.wdl` (206 lines) - Taxonomic classification and depletion
  - `sarscov2_illumina_full.wdl` (525 lines) - Complete SARS-CoV-2 analysis
  - `augur_from_assemblies.wdl` (203 lines) - Phylogenetic analysis with Nextstrain
  - `demux_deplete.wdl` (317 lines) - Demultiplexing and host depletion

**pipes/WDL/tasks/:**
- Purpose: Reusable task modules organized by functional domain
- Contains: 16 task module WDL files
- Key files:
  - `tasks_assembly.wdl` (1361 lines) - Assembly, alignment, scaffolding tasks
  - `tasks_nextstrain.wdl` (1898 lines) - Nextstrain/Augur phylogenetics
  - `tasks_ncbi.wdl` (1836 lines) - GenBank/SRA submission tasks
  - `tasks_read_utils.wdl` (578 lines) - BAM/FASTQ manipulation
  - `tasks_taxon_filter.wdl` (282 lines) - Kraken2 classification and filtering
  - `tasks_intrahost.wdl` (401 lines) - Variant calling (lofreq, iSNV)
  - `tasks_interhost.wdl` (603 lines) - Multi-sample phylogenetics
  - `tasks_sarscov2.wdl` (597 lines) - SARS-CoV-2 specific tools (Pangolin, Nextclade)
  - `tasks_reports.wdl` (951 lines) - QC metrics, plots, MultiQC
  - `tasks_utils.wdl` (1299 lines) - General utilities (file ops, data tables)
  - `tasks_demux.wdl` (1165 lines) - Illumina demultiplexing
  - `tasks_metagenomics.wdl` (954 lines) - Metagenomic classification tools
  - `tasks_terra.wdl` (784 lines) - Terra platform integrations
  - `tasks_ncbi_tools.wdl` (575 lines) - Additional NCBI utilities
  - `tasks_megablast.wdl` (410 lines) - BLAST operations
  - `tasks_16S_amplicon.wdl` (425 lines) - 16S amplicon analysis

**test/input/WDL/:**
- Purpose: Test input configurations for workflow validation
- Contains: JSON files with test parameters
- Key files:
  - `test_inputs-<workflow_name>-local.json` - Cromwell test configs
  - `miniwdl-local/test_inputs-<workflow_name>-local.json` - miniWDL test configs (symlinks)
  - `test_outputs-<workflow_name>-local.json` - Expected output validation files
  - `*-dnanexus.dx.json` - DNAnexus platform test configs

**github_actions_ci/:**
- Purpose: CI/CD scripts for validation, testing, and deployment
- Contains: Bash scripts invoked by GitHub Actions
- Key files:
  - `validate-wdl-womtool.sh` - WDL validation with Cromwell's womtool
  - `tests-miniwdl.sh` - Run workflow tests with miniWDL
  - `tests-cromwell.sh` - Run workflow tests with Cromwell
  - `build-dx.sh` - Deploy workflows to DNAnexus
  - `check-wdl-runtimes.sh` - Verify Docker versions match requirements
  - `version-wdl-runtimes.sh` - Update Docker versions in WDL files
  - `install-wdl.sh` - Install Cromwell and womtool for CI
  - `build-docs.sh` - Build Sphinx documentation

**docs/:**
- Purpose: ReadTheDocs documentation source
- Contains: Sphinx documentation files
- Key files: `conf.py` (Sphinx config)

**.github/workflows/:**
- Purpose: GitHub Actions CI/CD pipeline definitions
- Contains: `build.yml` - Main CI workflow

**pipes/cromwell/:**
- Purpose: Cromwell execution engine configurations
- Contains: Backend configurations for HPC/cloud execution

**pipes/dnax/:**
- Purpose: DNAnexus platform deployment utilities
- Contains: dxWDL launcher scripts

## Key File Locations

**Entry Points:**
- `pipes/WDL/workflows/*.wdl`: All workflow entry points (invoked via `miniwdl run` or Cromwell submit)

**Configuration:**
- `requirements-modules.txt`: Docker image version specifications
- `.dockstore.yml`: Workflow registry publication config
- `test/input/WDL/test_inputs-*.json`: Test parameter files
- `.github/workflows/build.yml`: CI/CD pipeline definition

**Core Logic:**
- `pipes/WDL/tasks/tasks_assembly.wdl`: Assembly and alignment operations
- `pipes/WDL/tasks/tasks_read_utils.wdl`: BAM/FASTQ manipulation
- `pipes/WDL/tasks/tasks_reports.wdl`: Metrics and QC reporting

**Testing:**
- `test/input/WDL/miniwdl-local/*.json`: miniWDL test configurations
- `github_actions_ci/tests-miniwdl.sh`: Test execution script
- `github_actions_ci/tests-cromwell.sh`: Cromwell test execution

**Documentation:**
- `AGENTS.md`: Comprehensive guidance for AI coding agents
- `README.md`: User-facing project overview
- `docs/`: Sphinx documentation source

## Naming Conventions

**Files:**
- Workflows: `<use_case>.wdl` (e.g., `assemble_refbased.wdl`, `classify_single.wdl`)
- Tasks: `tasks_<domain>.wdl` (e.g., `tasks_assembly.wdl`, `tasks_nextstrain.wdl`)
- Test inputs: `test_inputs-<workflow_name>-<platform>.json`
- Test outputs: `test_outputs-<workflow_name>-<platform>.json`
- CI scripts: `<action>-<tool>.sh` (e.g., `validate-wdl-womtool.sh`, `tests-miniwdl.sh`)

**Directories:**
- lowercase with underscores: `read_utils`, `WDL` (exception - uppercase for acronym)
- Platform-specific: `miniwdl-local`, `dnanexus.dx`

**WDL Entities:**
- Workflows: `snake_case` (e.g., `workflow assemble_refbased`)
- Tasks: `snake_case` (e.g., `task align_reads`)
- Variables: `snake_case` (e.g., `reads_unmapped_bam`)
- Namespaces: `snake_case` (e.g., `import "../tasks/tasks_assembly.wdl" as assembly`)

**Docker Images:**
- Format: `registry/org/name:version-flavor`
- Example: `ghcr.io/broadinstitute/viral-ngs:3.0.4-assemble`
- Pinned in: `requirements-modules.txt`

## Where to Add New Code

**New Workflow:**
- Primary code: `pipes/WDL/workflows/<use_case>.wdl`
- Test config: `test/input/WDL/test_inputs-<use_case>-local.json`
- Dockstore registration: Add entry to `.dockstore.yml`
- Import pattern: `import "../tasks/tasks_<domain>.wdl" as <namespace>`

**New Task:**
- Add to existing module: `pipes/WDL/tasks/tasks_<domain>.wdl`
- Or create new module: `pipes/WDL/tasks/tasks_<new_domain>.wdl`
- Specify Docker image from `requirements-modules.txt`
- Follow pattern: `task <name> { input {...} command <<< ... >>> output {...} runtime {...} }`

**New Test:**
- Input JSON: `test/input/WDL/test_inputs-<workflow>-local.json`
- Expected outputs: `test/input/WDL/test_outputs-<workflow>-local.json`
- miniWDL symlink: `test/input/WDL/miniwdl-local/test_inputs-<workflow>-local.json`

**New CI Script:**
- Location: `github_actions_ci/<action>-<tool>.sh`
- Make executable: `chmod +x github_actions_ci/<script>.sh`
- Invoke from: `.github/workflows/build.yml`

**New Documentation:**
- Sphinx docs: `docs/<topic>.rst`
- Agent guidance: Update `AGENTS.md`
- User README: Update `README.md`

## Special Directories

**.pixi/**
- Purpose: Pixi package manager environment
- Generated: Yes (by pixi install)
- Committed: No (in .gitignore)

**.planning/codebase/**
- Purpose: Auto-generated codebase analysis for GSD AI agents
- Generated: Yes (by `/gsd:map-codebase`)
- Committed: Yes (ephemeral planning state)

**test/input/WDL/miniwdl-local/**
- Purpose: miniWDL-specific test configurations
- Generated: Partially (mostly symlinks to parent directory)
- Committed: Yes (symlinks committed)

**pipes/dnax/dx-launcher/**
- Purpose: DNAnexus workflow launcher applet
- Generated: No
- Committed: Yes

## Import Patterns

**Workflow to Task Import:**
```wdl
import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_read_utils.wdl" as read_utils
```
- Pattern: Relative path from `workflows/` to `tasks/`
- All imports use `..` to go up one level

**Workflow to Workflow Import:**
```wdl
import "demux_deplete.wdl"
import "assemble_refbased.wdl"
```
- Pattern: Same-directory relative import (no path prefix)
- Used for composing sub-workflows

**Task Module Organization:**
- Each module contains multiple related tasks
- Tasks within a module can reference each other
- External references require import with namespace

## File Type Distribution

**Total WDL Files:** 109 (90 workflows + 16 task modules + 3 other)

**Total Test Configs:** 40+ JSON files across platforms

**Docker Images:** 10 distinct images with pinned versions

**CI Scripts:** 14 shell scripts

## Platform-Specific Considerations

**miniWDL:**
- Test inputs: `test/input/WDL/miniwdl-local/*.json`
- Execution: Local Docker, outputs to `_LAST/` symlink
- Validation: `miniwdl check`

**Cromwell:**
- Test inputs: `test/input/WDL/test_inputs-*-local.json`
- Execution: Configurable backends (local, GCP, AWS, HPC)
- Validation: `womtool validate`

**Terra:**
- Import: Via Dockstore TRS API from `.dockstore.yml`
- Execution: Cromwell on GCP
- Data tables: Integration via `tasks_terra.wdl`

**DNAnexus:**
- Deployment: Via `github_actions_ci/build-dx.sh` and dxWDL compiler
- Test inputs: `test/input/WDL/*-dnanexus.dx.json`
- Runtime: `dx_instance_type` in task runtime blocks

---

*Structure analysis: 2026-02-11*
