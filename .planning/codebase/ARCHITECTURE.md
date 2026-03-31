# Architecture

**Analysis Date:** 2026-03-31

## Pattern Overview

**Overall:** Modular task-based composition using WDL (Workflow Description Language)

**Key Characteristics:**
- Workflows decompose complex bioinformatics pipelines into reusable task modules
- Tasks are grouped by functional domain (assembly, classification, reports, etc.)
- Workflows import task modules as namespaced imports and orchestrate execution
- Parallel execution through WDL scatter/gather patterns for multi-sample and multi-input processing
- Execution platform agnostic: runs on miniWDL, Cromwell, Terra, DNAnexus, HPC systems
- Docker containerization for all compute tasks with pinned image versions

## Layers

**Workflow Layer:**
- Purpose: Orchestrate multi-step bioinformatics analyses by composing tasks and managing data flow
- Location: `pipes/WDL/workflows/`
- Contains: ~100 WDL workflow files (e.g., `assemble_refbased.wdl`, `classify_kraken2.wdl`, `demux_deplete.wdl`)
- Depends on: Task modules, imported workflows, WDL runtime
- Used by: End users via Dockstore, Terra, DNAnexus, or local execution with miniWDL/Cromwell

**Task Module Layer:**
- Purpose: Define atomic computation units (tools, scripts, commands) that operate on inputs and produce outputs
- Location: `pipes/WDL/tasks/`
- Contains: 16 task definition files organized by function:
  - `tasks_assembly.wdl` - Genome assembly and scaffolding (SPAdes, minimap2, bwa, novoalign, reference-based consensus)
  - `tasks_read_utils.wdl` - BAM/FASTQ utilities (merging, reheadering, duplicate removal, grouping)
  - `tasks_taxon_filter.wdl` - Taxonomic filtering and depletion (bmtagger, BLAST, bwa depletion)
  - `tasks_intrahost.wdl` - Intra-host variant detection (LoFreq iSNV analysis)
  - `tasks_interhost.wdl` - Multi-sample phylogenetic analysis (MAFFT, SNP calling, tree building)
  - `tasks_metagenomics.wdl` - Taxonomic classification (Kraken2, Kaiju, QIIME)
  - `tasks_ncbi.wdl` - GenBank submission preparation and metadata handling
  - `tasks_ncbi_tools.wdl` - SRA/BioSample submission tools
  - `tasks_nextstrain.wdl` - Nextstrain/Augur integration for phylogenetic visualization
  - `tasks_sarscov2.wdl` - SARS-CoV-2 specific tools (lineage assignment, VADR QC)
  - `tasks_reports.wdl` - QC metrics and visualization (coverage plots, alignment metrics)
  - `tasks_demux.wdl` - Illumina demultiplexing (BCL to FASTQ)
  - `tasks_terra.wdl` - Terra workspace integration utilities
  - `tasks_utils.wdl` - General utilities (TSV operations, file operations, archive handling)
  - `tasks_16S_amplicon.wdl` - 16S rRNA amplicon analysis
  - `tasks_megablast.wdl` - BLAST-based sequence search and classification
- Depends on: Docker images (viral-ngs flavors, bioinformatics tools)
- Used by: Workflows that import task modules via relative paths

**Integration Layer:**
- Purpose: Orchestrate complex multi-workflow pipelines for common analysis scenarios
- Location: `pipes/WDL/workflows/` (distinguished by composition patterns)
- Examples: `demux_deplete.wdl`, `assemble_denovo.wdl`, `sarscov2_illumina_full.wdl`
- Pattern: Import other workflows as sub-workflows (e.g., `assemble_refbased.wdl` imported into `assemble_denovo.wdl`)

**Testing/Configuration Layer:**
- Purpose: Define test inputs, expected outputs, and CI/CD validation
- Location: `test/input/WDL/`, `github_actions_ci/`
- Contains:
  - Test input JSON files (format: `test_inputs-{workflow_name}-local.json`)
  - Optional expected output JSON files
  - CI/CD scripts for validation and testing
  - Docker version pins in `requirements-modules.txt`

## Data Flow

**Reference-Based Assembly (assemble_refbased.wdl) Flow:**

1. **Alignment Phase** (scatter across input BAMs):
   - Call `assembly.align_reads` → align reads to reference genome
   - Call `assembly.ivar_trim` → trim PCR primers if bed file provided
   - Merge aligned BAMs across inputs if multiple BAMs

2. **Variant Detection Phase**:
   - Call `intrahost.lofreq` → detect intra-host SNVs against reference
   - Call `assembly.run_discordance` → compute concordance metrics
   - Call `reports.alignment_metrics` → calculate coverage and alignment stats

3. **Consensus Calling Phase**:
   - Call `assembly.refine_assembly_with_aligned_reads` → call consensus sequence
   - Filters by minimum coverage and allele frequency thresholds

4. **Self-Alignment Phase**:
   - Call `assembly.align_reads` → re-align reads to newly called assembly
   - Call `intrahost.lofreq` → detect variants against new assembly
   - Merge aligned BAMs if multiple input BAMs

5. **Output Generation**:
   - Emits: assembly_fasta, variant VCF, coverage metrics, alignment plots

**De Novo Assembly (assemble_denovo.wdl) Flow:**

1. **Read Processing** (scatter across input BAMs):
   - Optional: rename sample identifiers via `read_utils.merge_and_reheader_bams`
   - Optional: deplete host/contaminant reads via `taxon_filter.deplete_taxa`
   - Optional: filter to taxon-of-interest via `taxon_filter.filter_to_taxon`
   - Deduplicate reads via `read_utils.rmdup_ubam`

2. **Assembly Phase**:
   - Merge deduplicated reads across inputs
   - Call `assembly.assemble` → SPAdes de novo assembly
   - Call `assembly.scaffold` → scaffold contigs against reference genome(s)

3. **Refinement Phase**:
   - Compose sub-workflow `assemble_refbased` → polish assembly with aligned reads

4. **Header Annotation**:
   - Optional: rename fasta headers via `ncbi.rename_fasta_header`

**Demux + Depletion (demux_deplete.wdl) Flow:**

1. **Demultiplexing** (scatter across flowcell lanes):
   - Extract BCL files from tarball, demultiplex to FASTQ per sample
   - Collapse barcode duplicates, generate demultiplexing metrics
   - Convert FASTQ to unaligned BAM format (with read group headers)

2. **QC & Depletion**:
   - Generate FastQC per sample
   - Deplete host/contaminants via spike-in detection and depletion strategies
   - Merge MultiQC reports across all samples

3. **Metadata Management**:
   - Map samples through rename tables, merge biosample attributes
   - Generate SRA submission metadata TSV

**SARS-CoV-2 Full Analysis (sarscov2_illumina_full.wdl) Flow:**

1. Demultiplex and deplete reads via `demux_deplete` sub-workflow
2. Group BAMs by sample via `read_utils.group_bams_by_sample`
3. Parallel per-sample analysis (scatter):
   - Assemble via `assemble_refbased` (with amplicon trimming if applicable)
   - Fetch biosample metadata
   - Call variants, assess genome quality
   - Perform lineage assignment via `sarscov2_lineages` sub-workflow
4. Generate data release package with all samples

**State Management:**
- No persistent state within workflows; all state flows through file I/O
- BAM files carry sample metadata in read group headers (SM, LB, PL tags)
- Tab-separated values (TSV) used for metadata tables and mappings
- JSON used for structured metadata (demux sample maps, Terra workspace integration)

## Key Abstractions

**Task Abstraction:**
- Purpose: Encapsulate a bioinformatics tool or computational step
- Examples: `pipes/WDL/tasks/tasks_assembly.wdl` contains `assemble`, `scaffold`, `align_reads`, `refine_assembly_with_aligned_reads` tasks
- Pattern: Each task declares inputs (with type and metadata), command block (bash/Python), and outputs (with type and name)
- Docker runtime specified per task; resource requirements (memory, CPU, disk) parameterized

**Workflow Abstraction:**
- Purpose: Compose tasks into executable pipeline; manage parallelization and conditional execution
- Examples: `pipes/WDL/workflows/assemble_refbased.wdl`, `pipes/WDL/workflows/classify_kraken2.wdl`
- Pattern: Import task modules, scatter over inputs, call tasks, produce workflow outputs
- Metadata includes description, author, email; parameter_meta provides documentation

**Module Import Pattern:**
- Purpose: Namespace task definitions and enable code reuse
- Usage: `import "../tasks/tasks_assembly.wdl" as assembly`
- Pattern: Enables qualified task calls like `assembly.align_reads` and `assembly.assemble`
- Relative imports: Workflows reference tasks via relative paths (e.g., `../tasks/`)

**Scatter/Gather Pattern:**
- Purpose: Parallelize computation across multiple inputs (samples, lanes, reference genomes)
- Example in `assemble_refbased.wdl`:
  ```wdl
  scatter(reads_unmapped_bam in reads_unmapped_bams) {
      call assembly.align_reads as align_to_ref { input: reads_unmapped_bam = reads_unmapped_bam }
  }
  File aligned_bam = select_first([merge_align_to_ref.out_bam, ivar_trim.aligned_trimmed_bam[0]])
  ```
- Enables processing multiple samples in a single workflow submission

**Conditional Execution Pattern:**
- Purpose: Skip steps based on user input or data presence
- Example in `assemble_refbased.wdl`:
  ```wdl
  if(length(reads_unmapped_bams)>1) {
      call read_utils.merge_and_reheader_bams as merge_align_to_ref { ... }
  }
  ```
- Allows optional trimming, depletion, metadata transformation

**Workflow Composition Pattern:**
- Purpose: Nest workflows as sub-workflows for complex multi-stage pipelines
- Example in `assemble_denovo.wdl`: `call assemble_refbased.assemble_refbased as refine { ... }`
- Enables progressive refinement and layered analysis

## Entry Points

**Workflow Entry Points:**
- Location: All `.wdl` files in `pipes/WDL/workflows/`
- Triggers: Submitted via Dockstore, Terra, DNAnexus, or `miniwdl run` command
- Responsibilities: Declare inputs, orchestrate task calls, produce outputs

**Top-Level Workflow Examples:**
- `pipes/WDL/workflows/assemble_refbased.wdl` - Consensus assembly from aligned reads
  - Primary inputs: reads_unmapped_bams (Array[File]), reference_fasta (File)

- `pipes/WDL/workflows/classify_kraken2.wdl` - Taxonomic classification
  - Primary input: kraken2 database and reads

- `pipes/WDL/workflows/demux_deplete.wdl` - Demultiplexing and depletion
  - Primary input: flowcell_tgz (BCL tarball), samplesheets (Array[File])

- `pipes/WDL/workflows/sarscov2_illumina_full.wdl` - Full SARS-CoV-2 end-to-end pipeline
  - Primary inputs: flowcell_tgz, reference_fasta, amplicon_bed_prefix

**CI/CD Entry Points:**
- Location: `github_actions_ci/`
- Scripts:
  - `validate-wdl-womtool.sh` - Validate syntax with Cromwell's womtool
  - `tests-miniwdl.sh` - Run workflow tests with miniWDL
  - `tests-cromwell.sh` - Run workflow tests with Cromwell
  - `check-wdl-runtimes.sh` - Verify Docker image versions match requirements

## Error Handling

**Strategy:** Graceful degradation with optional success parameters and task retry logic

**Patterns:**
- **Task-level retries:** `maxRetries: 2` in task runtimes (allows transient GCP preemption recovery)
- **Always-succeed tasks:** `always_succeed=true` parameter in assembly tasks (outputs empty file if assembly fails)
- **Conditional task execution:** Skip non-critical steps if inputs missing (e.g., optional primer trimming)
- **Optional outputs:** Use `select_first()` to pick from conditional task outputs with fallback
- **Validation in commands:** Shell scripts use `set -e -o pipefail` to fail on command errors
- **Metadata validation:** Samplesheet parsing validates required columns; illumina_demux validates barcode structure

**Error Detection:**
- No centralized error handling framework; errors propagate through task exit codes
- stderr/stdout logs available for post-mortem analysis
- Terra MCP provides `get_batch_job_status` for infrastructure-level failure diagnostics (preemption, image pull failures, resource exhaustion)

## Cross-Cutting Concerns

**Logging:** Uses shell `set -ex` for command tracing; Python logging in embedded Python scripts; viral-ngs tools emit DEBUG level logs; output captured in task stderr

**Validation:**
- **Samplesheet validation:** Tasks validate TSV format and required columns
- **BAM validation:** `samtools view` used to verify BAM structure and read group headers
- **FASTA validation:** File extensions and format verification before processing
- **Genome quality:** Minimum coverage and unambiguous base thresholds control assembly quality gates

**Authentication:**
- No authentication within WDL; reliant on runtime platform:
  - **Terra:** Uses workspace service accounts and workspace-scoped GCS bucket access
  - **GCP (Cromwell):** Uses service account credentials for GCS access
  - **DNAnexus:** Uses platform-native auth for file staging
  - **Local (miniWDL):** Direct filesystem access

**Resource Management:**
- **Disk:** Parameterized per task; calculated as multiple of input file sizes + overhead
- **Memory:** Parameterized (defaults conservative, e.g., 32GB for assembly); tasks read cgroup limits to detect available memory at runtime
- **CPU:** Parameterized per task (defaults 1-8 depending on tool parallelism)
- **Runtime platform profiles:** dx_instance_type for DNAnexus, standard Cromwell resource declarations for others

**File I/O & Staging:**
- Implicit file staging/localization handled by WDL runtime (miniWDL, Cromwell, Terra, DNAnexus)
- `localization_optional: true` used for large BAM files to avoid unnecessary copies (stream mode)
- GCS bucket operations use `gcloud storage` CLI for efficient batch queries

---

*Architecture analysis: 2026-03-31*
