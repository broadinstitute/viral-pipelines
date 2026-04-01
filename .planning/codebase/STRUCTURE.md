# Codebase Structure

## Directory Layout

```
viral-pipelines/
├── pipes/                        # Core pipeline definitions
│   ├── WDL/
│   │   ├── workflows/            # 101 top-level workflow WDL files
│   │   └── tasks/                # 16 task module WDL files
│   ├── cromwell/                 # Cromwell backend configs
│   │   ├── cromwell.local-github_actions.conf
│   │   ├── cromwell.local-travis.conf
│   │   └── cromwell.gcid-viral-seq.conf
│   └── dnax/                     # DNAnexus-specific configs
├── test/
│   └── input/
│       ├── WDL/
│       │   ├── miniwdl-local/    # 19 miniWDL test input/output JSON files
│       │   └── cromwell-local/   # Cromwell test input JSON files
│       ├── genbank/              # GenBank submission test data
│       ├── TestSplitcodeDemuxFastqs/
│       └── zika-tutorial/
│       [+ various reference FASTA, BAM, FASTQ, JSON files]
├── github_actions_ci/            # CI/CD shell scripts
│   ├── tests-miniwdl.sh          # miniWDL integration test runner
│   ├── tests-cromwell.sh         # Cromwell integration test runner
│   ├── tests-dx.sh               # DNAnexus test runner
│   ├── install-wdl.sh            # WDL toolchain installer
│   ├── install-pip-docs.sh       # Docs dependency installer
│   ├── validate-wdl-womtool.sh   # WDL validation via Womtool
│   ├── check-wdl-runtimes.sh     # Docker version checker
│   ├── build-dx.sh               # DNAnexus build script
│   ├── build-docs.sh             # Documentation builder
│   ├── dockstoreyml.sh           # Dockstore YAML generator
│   ├── list-docker-tags.sh       # Docker tag lister
│   └── relative-wdl-paths.sh    # Path validator for WDL imports
├── docs/                         # ReadTheDocs documentation source
├── .github/
│   └── workflows/
│       └── build.yml             # GitHub Actions CI definition (16.6K)
├── pixi.toml                     # Pixi environment definition
├── pixi.lock                     # Locked pixi dependencies
├── requirements-modules.txt      # Docker module versions (10 containers)
├── .dockstore.yml                # Dockstore workflow registry config (12.5K)
├── .coveragerc                   # Python coverage config
├── AGENTS.md                     # Codebase conventions for AI agents (13KB)
├── CLAUDE.md                     # Claude Code instructions
└── README.md
```

## Key Locations

### Workflow WDL Files
`pipes/WDL/workflows/` — 101 top-level workflow files.

**Assembly:**
- `pipes/WDL/workflows/assemble_refbased.wdl` — Reference-based genome assembly
- `pipes/WDL/workflows/assemble_denovo.wdl` — De novo assembly
- `pipes/WDL/workflows/assemble_denovo_metagenomic.wdl` — Metagenomic de novo

**Classification:**
- `pipes/WDL/workflows/classify_kraken2.wdl`
- `pipes/WDL/workflows/classify_krakenuniq.wdl`
- `pipes/WDL/workflows/classify_kaiju.wdl`
- `pipes/WDL/workflows/classify_kallisto_single.wdl` / `classify_kallisto_multi.wdl`
- `pipes/WDL/workflows/classify_virnucpro_single.wdl` / `classify_virnucpro_multi.wdl`
- `pipes/WDL/workflows/classify_multi.wdl` — Aggregated classifier

**Demultiplexing / Sample Prep:**
- `pipes/WDL/workflows/demux_plus.wdl`
- `pipes/WDL/workflows/load_illumina_fastqs.wdl`
- `pipes/WDL/workflows/load_illumina_fastqs_deplete.wdl`

**SARS-CoV-2 Specific:**
- `pipes/WDL/workflows/sarscov2_illumina_full.wdl` — End-to-end SC2 pipeline
- `pipes/WDL/workflows/sarscov2_lineages.wdl`
- `pipes/WDL/workflows/sarscov2_nextstrain.wdl`
- `pipes/WDL/workflows/sarscov2_genbank.wdl`
- `pipes/WDL/workflows/sarscov2_data_release.wdl`
- `pipes/WDL/workflows/sarscov2_nextclade_multi.wdl`

**Phylogenetics / Nextstrain:**
- `pipes/WDL/workflows/augur_from_assemblies.wdl`
- `pipes/WDL/workflows/augur_from_msa.wdl`
- `pipes/WDL/workflows/augur_from_beast_mcc.wdl`

**Alignment:**
- `pipes/WDL/workflows/align_and_count.wdl`
- `pipes/WDL/workflows/align_and_plot.wdl`
- `pipes/WDL/workflows/align_and_generate_PAF.wdl`

**NCBI Submissions:**
- `pipes/WDL/workflows/submit_biosample.wdl`
- `pipes/WDL/workflows/submit_genbank.wdl`
- `pipes/WDL/workflows/submit_sra.wdl`

### Task Module WDL Files
`pipes/WDL/tasks/` — 16 task modules (~10-76KB each):

| File | Domain |
|------|--------|
| `tasks_assembly.wdl` (55.4K) | Assembly, consensus, coverage |
| `tasks_nextstrain.wdl` (76.7K) | Augur, auspice, phylogenetics |
| `tasks_ncbi.wdl` (72.7K) | GenBank, SRA, BioSample submissions |
| `tasks_metagenomics.wdl` (46.7K) | Kraken2, KrakenUniq, Kaiju, Virnucpro |
| `tasks_demux.wdl` (44.1K) | Illumina demultiplexing, FASTQ processing |
| `tasks_utils.wdl` (43.4K) | General utilities |
| `tasks_terra.wdl` (32.7K) | Terra platform integration |
| `tasks_sarscov2.wdl` (28.6K) | SC2-specific tasks |
| `tasks_reports.wdl` (38.3K) | MultiQC, reporting |
| `tasks_ncbi_tools.wdl` (22.2K) | NCBI tools |
| `tasks_interhost.wdl` (19.5K) | Multi-sample analysis |
| `tasks_read_utils.wdl` (21.3K) | BAM/FASTQ utilities |
| `tasks_megablast.wdl` (22.0K) | BLAST-based classification |
| `tasks_intrahost.wdl` (14.6K) | Within-host variant analysis |
| `tasks_taxon_filter.wdl` (10.8K) | Depletion and filtering |
| `tasks_16S_amplicon.wdl` (14.3K) | 16S amplicon processing |

### CI/CD Entry Points
`github_actions_ci/tests-miniwdl.sh` — Primary local test runner
`github_actions_ci/tests-cromwell.sh` — Cromwell-based test runner
`.github/workflows/build.yml` — GitHub Actions CI configuration

### Test Data
`test/input/WDL/miniwdl-local/` — 19 JSON files (12 input, 5 output)
`test/input/WDL/cromwell-local/` — Cromwell-specific input JSONs

## Naming Conventions

### WDL Files
- **Workflows:** snake_case verbs — `assemble_refbased.wdl`, `classify_kraken2.wdl`
- **Task modules:** `tasks_<domain>.wdl` prefix — `tasks_assembly.wdl`
- **Test inputs:** `test_inputs-<workflow_name>-<platform>.json`
- **Test outputs:** `test_outputs-<workflow_name>-<platform>.json`

### WDL Identifiers
- **Workflow names:** PascalCase matching filename — `AssembleRefbased`, `ClassifyKraken2`
- **Task names:** PascalCase — `AssembleRefbasedWithViral`, `RunPangolin`
- **Input/output parameters:** snake_case — `reads_bam`, `assembly_fasta`
- **Runtime variables:** snake_case — `docker`, `memory`, `cpu`, `disk_size`

### Docker Image References
Defined in `requirements-modules.txt` and referenced in task runtime blocks:
```
broadinstitute/viral-ngs=3.0.6
broadinstitute/read-qc-tools=1.0.1
nextstrain/nextclade=3.18.1
```

## Import Pattern

WDL workflows import task modules using relative paths:
```wdl
import "../tasks/tasks_assembly.wdl" as tasks_assembly
import "../tasks/tasks_reports.wdl" as tasks_reports
```

The `github_actions_ci/relative-wdl-paths.sh` script validates that all imports use relative paths (required for multi-platform compatibility with Terra, DNAnexus, Cromwell).
