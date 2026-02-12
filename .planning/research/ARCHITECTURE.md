# Architecture Research

**Domain:** WDL bioinformatics pipeline for genomad integration
**Researched:** 2026-02-12
**Confidence:** HIGH

## Standard Architecture

### System Overview

```
┌─────────────────────────────────────────────────────────────┐
│                     Workflow Layer                           │
│  ┌──────────────────┐         ┌──────────────────┐          │
│  │  genomad_single  │         │  genomad_multi   │          │
│  │   (per-sample)   │         │  (batch scatter) │          │
│  └────────┬─────────┘         └────────┬─────────┘          │
│           │                            │                     │
│           └────────────┬───────────────┘                     │
├────────────────────────┴─────────────────────────────────────┤
│                      Task Layer                              │
│  ┌─────────────────────────────────────────────────────┐    │
│  │  tasks_metagenomics.wdl (new genomad tasks)         │    │
│  │  - genomad_end_to_end                               │    │
│  │  - summarize_genomad_results (optional)             │    │
│  └─────────────────────────────────────────────────────┘    │
├──────────────────────────────────────────────────────────────┤
│                   Container Layer                            │
│  ┌──────────────────┐  ┌──────────────────┐                 │
│  │ genomad Docker   │  │  viral-ngs       │                 │
│  │ (primary tool)   │  │  (BAM utilities) │                 │
│  └──────────────────┘  └──────────────────┘                 │
├──────────────────────────────────────────────────────────────┤
│                    Data Layer                                │
│  ┌─────────────┐  ┌─────────────┐  ┌──────────────┐         │
│  │  Input      │  │  Database   │  │   Output     │         │
│  │  FASTA/BAM  │  │  (~5GB)     │  │  TSV/FASTA   │         │
│  └─────────────┘  └─────────────┘  └──────────────┘         │
└──────────────────────────────────────────────────────────────┘
```

### Component Responsibilities

| Component | Responsibility | Typical Implementation |
|-----------|----------------|------------------------|
| tasks_metagenomics.wdl | WDL task definitions for genomad execution | Docker-based tasks with resource specs |
| genomad_single.wdl | Single-sample workflow orchestration | Imports tasks, calls genomad on one sample |
| genomad_multi.wdl | Multi-sample batch workflow | Scatter-gather pattern over sample array |
| Test fixtures | Validation and CI/CD integration | JSON input files in test/input/WDL/ |
| Docker containers | Isolated execution environments | genomad official image + viral-ngs for utilities |

## Recommended Project Structure

```
pipes/WDL/
├── tasks/
│   └── tasks_metagenomics.wdl       # Add genomad tasks here (existing file)
│       ├── task genomad_end_to_end  # Primary analysis task
│       └── task summarize_genomad   # Optional: aggregate results across samples
├── workflows/
│   ├── genomad_single.wdl           # Single-sample workflow
│   └── genomad_multi.wdl            # Multi-sample workflow (optional Phase 2)
test/input/WDL/
├── test_inputs-genomad_single-local.json      # miniwdl test inputs
├── test_outputs-genomad_single-local.json     # Expected outputs (optional)
├── genomad_db-minitest.tar.zst                # Small test database (~100MB)
└── test.contigs.fasta                         # Small test input file
```

### Structure Rationale

- **tasks_metagenomics.wdl:** Genomad fits naturally alongside kraken2, kaiju, and other classification tools. All are metagenomics/viral classification utilities operating on similar inputs (FASTA/BAM).
- **workflows/genomad_single.wdl:** Follows established pattern (classify_single.wdl, nextclade_single.wdl). Single-sample workflows are simpler, easier to test, and sufficient for initial adoption.
- **workflows/genomad_multi.wdl:** Deferred to Phase 2. Multi-sample workflows add scatter-gather complexity and are built after single-sample validation.
- **test/input/WDL/:** Tests co-locate with existing workflow tests. Test database must be tiny for CI performance (1000x smaller than production 5GB database).

## Architectural Patterns

### Pattern 1: Task Module Organization

**What:** Group related bioinformatics tools into task modules by functional category (assembly, classification, read processing, etc.).

**When to use:** When adding a new tool that shares a functional domain with existing tasks.

**Trade-offs:**
- PRO: Easy to find related tasks, clear imports in workflows
- PRO: Consistent with existing codebase (15 task modules organized this way)
- CON: Large task files can become unwieldy (tasks_metagenomics.wdl already has 10+ tasks)

**Example:**
```wdl
# pipes/WDL/tasks/tasks_metagenomics.wdl
version 1.0

task kraken2 { ... }
task krakenuniq { ... }
task kaiju { ... }
task genomad_end_to_end { ... }  # Add genomad here
```

### Pattern 2: Single/Multi Workflow Split

**What:** Separate workflows for single-sample (simple) and multi-sample (scatter-gather) execution.

**When to use:** When a tool needs to run on both individual samples and batches.

**Trade-offs:**
- PRO: Single-sample workflows are simpler, easier to debug, better for Terra data tables
- PRO: Multi-sample workflows batch database loading (critical for genomad's 5GB DB)
- CON: Code duplication between single and multi variants
- CON: Two workflows to maintain and test

**Example:**
```wdl
# genomad_single.wdl: processes one sample
workflow genomad_single {
  input {
    File   contigs_fasta
    File   genomad_db_tgz
  }
  call metagenomics.genomad_end_to_end {
    input: contigs = contigs_fasta, database = genomad_db_tgz
  }
  output {
    File virus_summary = genomad_end_to_end.virus_summary_tsv
  }
}

# genomad_multi.wdl: scatter over sample array
workflow genomad_multi {
  input {
    Array[File]+ contigs_fastas
    File         genomad_db_tgz
  }
  scatter(contigs in contigs_fastas) {
    call metagenomics.genomad_end_to_end {
      input: contigs = contigs, database = genomad_db_tgz
    }
  }
  output {
    Array[File] virus_summaries = genomad_end_to_end.virus_summary_tsv
  }
}
```

### Pattern 3: Database Decompression in Task Command

**What:** Decompress tarball databases inside task command blocks, not as separate WDL tasks.

**When to use:** Always, for databases that are tarballs (kraken2, genomad, etc.).

**Trade-offs:**
- PRO: Matches existing pattern used by kraken2, krakenuniq (consistency)
- PRO: Avoids call-caching issues and extra task overhead
- PRO: Database decompression time included in task timing metrics
- CON: Slightly more complex command blocks

**Example:**
```wdl
task genomad_end_to_end {
  command <<<
    DB_DIR=$(mktemp -d --suffix _db)

    # Decompress database tarball
    read_utils extract_tarball \
      "~{genomad_db_tgz}" $DB_DIR \
      --loglevel=DEBUG

    # Run genomad with decompressed database
    genomad end-to-end \
      "~{contigs_fasta}" \
      output_dir \
      $DB_DIR \
      --threads ~{cpu}
  >>>
}
```

### Pattern 4: Resource Autoscaling

**What:** Calculate CPU, memory, and disk requirements based on input file sizes and database sizes.

**When to use:** For tasks with highly variable resource needs depending on input size.

**Trade-offs:**
- PRO: Avoids over-provisioning for small inputs (saves cost)
- PRO: Avoids under-provisioning for large inputs (prevents failures)
- PRO: Matches patterns in deplete_taxa, kraken2 tasks (consistency)
- CON: More complex resource calculation logic in WDL

**Example:**
```wdl
task genomad_end_to_end {
  input {
    File   contigs_fasta
    File   genomad_db_tgz
    Int?   cpu
    Int?   machine_mem_gb
  }

  # Autoscale CPU based on input size (2 CPUs per GB, min 4, max 32)
  Float input_size_gb = size(contigs_fasta, "GB")
  Int   cpu_actual = select_first([cpu,
    if input_size_gb < 2.0 then 4
    else if input_size_gb > 16.0 then 32
    else ceil(input_size_gb * 2)])

  # Memory scales with CPU at 2x ratio
  Int machine_mem_gb_actual = select_first([machine_mem_gb, cpu_actual * 2])

  # Disk: 10x input + 3x database + 50GB overhead, rounded to 750GB (SSD pairs)
  Int disk_size = ceil((10 * input_size_gb + 3 * size(genomad_db_tgz, "GB") + 50) / 750.0) * 750

  runtime {
    cpu: cpu_actual
    memory: "~{machine_mem_gb_actual} GB"
    disks: "local-disk ~{disk_size} LOCAL"
  }
}
```

## Data Flow

### Single-Sample Flow

```
Input FASTA/BAM
    ↓
[Convert BAM to FASTA if needed] (optional preprocessing)
    ↓
[genomad_end_to_end task]
    ├─ Decompress database tarball → tmpdir
    ├─ Execute: genomad end-to-end INPUT OUTPUT_DIR DATABASE
    └─ Extract summary files from output_dir/PREFIX_summary/
    ↓
Output Files:
    ├─ virus_summary.tsv     (primary classification results)
    ├─ plasmid_summary.tsv   (plasmid results)
    ├─ virus.fna             (identified virus sequences)
    ├─ plasmid.fna           (identified plasmid sequences)
    └─ taxonomy.tsv          (taxonomic assignments in PREFIX_annotate/)
```

### Multi-Sample Flow (Phase 2)

```
Array[File] contigs_fastas
    ↓
scatter(contigs in contigs_fastas) {
    ├─ genomad_end_to_end(contigs[0])
    ├─ genomad_end_to_end(contigs[1])
    └─ genomad_end_to_end(contigs[N])
} → Array[outputs]
    ↓
[Optional: aggregate summary reports]
    ↓
Output Arrays:
    ├─ Array[File] virus_summaries
    ├─ Array[File] plasmid_summaries
    └─ File        merged_summary (if aggregated)
```

### Key Data Flows

1. **FASTA input:** Genomad requires FASTA format. If workflow receives BAM, must convert using `samtools fasta` or `read_utils.py bam_to_fasta`.

2. **Database handling:** 5GB database tarball is localized once per task, decompressed to tmpdir, used for execution, then discarded (not part of output).

3. **Output directory structure:** Genomad creates multi-level directory structure (PREFIX_annotate, PREFIX_summary, etc.). Task must extract specific files from subdirectories for WDL outputs.

4. **File globbing:** Use WDL `glob()` to extract multiple related files from output directory, or explicit file paths for key outputs.

## Component Boundaries

### Task Module: tasks_metagenomics.wdl

**Boundary:** Low-level WDL task definitions that wrap individual tool executions.

**Responsibilities:**
- Define task inputs/outputs with parameter metadata
- Specify Docker container and resource requirements
- Write command block to execute genomad
- Extract output files from genomad directory structure
- Calculate resource autoscaling logic

**Does NOT:**
- Orchestrate multi-step workflows
- Handle conditional logic based on results
- Aggregate results across samples

**Interface:**
```wdl
task genomad_end_to_end {
  input {
    File   contigs_fasta           # FASTA format sequences
    File   genomad_db_tgz          # Database tarball (~5GB)
    Int?   min_score               # Optional: score threshold
    Int?   threads                 # Optional: override CPU count
    String docker                  # Container image
  }
  output {
    File   virus_summary_tsv       # Main results: virus hits
    File   plasmid_summary_tsv     # Plasmid hits
    File   virus_sequences_fna     # Extracted virus sequences
    File   plasmid_sequences_fna   # Extracted plasmid sequences
    File   taxonomy_tsv            # Taxonomic assignments
    String genomad_version         # Tool version string
    Int    max_ram_gb              # Peak memory usage
    Int    runtime_sec             # Execution time
  }
}
```

### Workflow: genomad_single.wdl

**Boundary:** Single-sample orchestration and optional pre/post-processing.

**Responsibilities:**
- Accept high-level inputs (may be FASTA or BAM)
- Convert BAM to FASTA if needed (conditional block)
- Call genomad_end_to_end task with appropriate inputs
- Expose user-friendly outputs with clear naming
- Provide workflow metadata (description, author)

**Does NOT:**
- Process multiple samples (that's genomad_multi.wdl)
- Aggregate or compare results across runs

**Interface:**
```wdl
workflow genomad_single {
  input {
    File   sequences_bam_or_fasta  # Accept either format
    File   genomad_db_tgz
    Int?   min_virus_score
  }
  output {
    File   virus_summary
    File   plasmid_summary
    Int    num_viruses_found       # Parsed from summary
    Int    num_plasmids_found      # Parsed from summary
    String genomad_version
  }
}
```

### Workflow: genomad_multi.wdl (Phase 2)

**Boundary:** Multi-sample batch processing with scatter-gather pattern.

**Responsibilities:**
- Accept array of input files
- Scatter genomad_end_to_end task over all samples in parallel
- Gather array outputs for downstream use
- Optionally aggregate summary statistics across samples

**Does NOT:**
- Process single samples (use genomad_single.wdl for that)
- Perform cross-sample analysis or comparison

**Interface:**
```wdl
workflow genomad_multi {
  input {
    Array[File]+ sequences_fastas  # Array of FASTA files
    File         genomad_db_tgz
  }
  output {
    Array[File] virus_summaries    # One per input sample
    Array[File] plasmid_summaries
    Array[Int]  num_viruses_found
    File?       merged_summary     # Optional aggregation
  }
}
```

### Test Configuration Files

**Boundary:** CI/CD integration and workflow validation.

**Responsibilities:**
- Provide minimal test inputs that execute quickly (<2 minutes)
- Use tiny test database (~100MB, not production 5GB)
- Verify workflow executes successfully (exit code 0)
- Optionally validate output structure (test_outputs-*.json)

**Does NOT:**
- Validate biological correctness of results
- Test all genomad features exhaustively

**Files:**
```
test/input/WDL/test_inputs-genomad_single-local.json
{
  "genomad_single.sequences_bam_or_fasta": "test/input/test.contigs.fasta",
  "genomad_single.genomad_db_tgz": "test/input/genomad_db-minitest.tar.zst"
}
```

## Build Order

### Phase 1: Core Task + Single-Sample Workflow

**Order:** Task → Workflow → Tests → Validation

1. **Task implementation** (`tasks_metagenomics.wdl`)
   - Add `task genomad_end_to_end` to existing file
   - Implement command block with database decompression
   - Define inputs/outputs with parameter_meta
   - Test task in isolation: `miniwdl run --task genomad_end_to_end ...`

2. **Single-sample workflow** (`genomad_single.wdl`)
   - Create new workflow file
   - Import tasks_metagenomics.wdl
   - Implement BAM-to-FASTA conversion (if needed)
   - Call genomad_end_to_end task
   - Wire up workflow inputs/outputs

3. **Test fixtures**
   - Create test input JSON
   - Acquire or create miniature test database (<100MB)
   - Create tiny test FASTA file (<1MB, 10-100 sequences)
   - Validate: `miniwdl run genomad_single.wdl -i test_inputs-genomad_single-local.json`

4. **WDL validation**
   - Run `miniwdl check` on both task and workflow files
   - Fix any syntax or type errors
   - Verify imports resolve correctly

5. **Integration testing**
   - Add test to CI/CD scripts (github_actions_ci/tests-miniwdl.sh)
   - Run full CI suite locally
   - Verify test completes in <2 minutes

### Phase 2: Multi-Sample Workflow (Optional)

**Order:** Multi workflow → Multi tests → Aggregation task (optional)

6. **Multi-sample workflow** (`genomad_multi.wdl`)
   - Copy genomad_single.wdl as template
   - Change input to `Array[File]+ sequences_fastas`
   - Add scatter block around genomad_end_to_end call
   - Update outputs to be arrays

7. **Multi-sample tests**
   - Create test_inputs-genomad_multi-local.json
   - Use array of 2-3 small test FASTA files
   - Validate scatter-gather execution

8. **Optional: Summary aggregation task**
   - Add `task summarize_genomad` to tasks_metagenomics.wdl
   - Aggregate stats across Array[File] summaries
   - Produce single merged report

### Dependency Rationale

- **Task before workflow:** Workflow imports and calls task, so task must exist first.
- **Single before multi:** Multi-sample workflow reuses same task, but adds complexity. Validate simple case first.
- **Tests after implementation:** Can't test what doesn't exist, but tests should be added before declaring a component "done".
- **Aggregation task last:** Optional nice-to-have. Core functionality works without it.

### Parallel Work Opportunities

These can be developed in parallel once task is complete:
- Single-sample workflow + its tests
- Docker container updates (if custom image needed)
- Documentation updates

These must be sequential:
- Task → Workflow (workflow depends on task)
- Implementation → Tests (tests depend on implementation)
- Single workflow → Multi workflow (multi builds on single pattern)

## Scaling Considerations

| Scale | Architecture Adjustments |
|-------|--------------------------|
| 1-10 samples | Use genomad_single.wdl, run sequentially or as separate Terra submissions. No scatter needed. |
| 10-100 samples | Use genomad_multi.wdl. Scatter is efficient. Database localized once per task (inefficient, but acceptable). |
| 100-1000 samples | Optimize database handling: consider pre-decompressing database to a shared bucket location, mount as read-only volume instead of re-decompressing per task. |
| 1000+ samples | Consider database sharding strategies or alternative architectures (e.g., persistent Cromwell server with shared disk). Database I/O becomes bottleneck. |

### Scaling Priorities

1. **First bottleneck: Database re-decompression overhead**
   - **What breaks:** For 100+ samples in scatter, each task re-downloads and decompresses 5GB database independently.
   - **How to fix:** Pre-decompress database to GCS bucket, pass decompressed directory path instead of tarball. Requires Docker image with gcsfuse or similar.
   - **When to optimize:** After validating basic functionality, if running >50 samples regularly.

2. **Second bottleneck: Disk I/O for large input assemblies**
   - **What breaks:** Genomad generates large intermediate files (gene predictions, marker matches). Large assemblies (>10MB) may exceed autoscaled disk.
   - **How to fix:** Tune disk autoscaling formula based on real usage patterns. Add `--cleanup` flags to genomad to remove intermediates.
   - **When to optimize:** If seeing disk quota failures in production.

3. **Third bottleneck: Memory for large databases**
   - **What breaks:** Genomad loads database into memory. Custom larger databases (>10GB) may OOM.
   - **How to fix:** Increase `machine_mem_gb` parameter, or use database subsetting strategies.
   - **When to optimize:** Only if using custom databases larger than default.

## Anti-Patterns

### Anti-Pattern 1: Creating Separate Task Module for Single Tool

**What people do:** Create `tasks_genomad.wdl` as a new standalone file for just genomad tasks.

**Why it's wrong:**
- Violates existing organizational pattern (tools grouped by function, not one file per tool)
- Creates proliferation of task files (would eventually have 50+ task files)
- Genomad is functionally a metagenomics/classification tool, belongs with kraken2/kaiju

**Do this instead:** Add genomad tasks to existing `tasks_metagenomics.wdl` file, following precedent of kraken2, krakenuniq, kaiju.

### Anti-Pattern 2: Passing Decompressed Database as Workflow Input

**What people do:** Decompress genomad database outside WDL, pass decompressed directory path as input.

**Why it's wrong:**
- WDL cannot localize directories, only individual files
- Would require passing 10,000+ individual database files as Array[File] (unwieldy)
- Breaks portability (assumes specific filesystem structure)

**Do this instead:** Pass database as tarball, decompress inside task command block (matches kraken2 pattern).

### Anti-Pattern 3: Hardcoding Docker Image Versions in Task Definitions

**What people do:** Put specific version in task: `String docker = "genomad/genomad:1.8.0"`

**Why it's wrong:**
- Makes version updates require editing WDL files
- Inconsistent with viral-pipelines pattern (Docker param has default but can be overridden)
- Version management should be in requirements-modules.txt + version-wdl-runtimes.sh

**Do this instead:**
```wdl
task genomad_end_to_end {
  input {
    String docker = "genomad/genomad:1.8.0"  # Default, but can override
  }
  runtime {
    docker: docker  # Use input param, not hardcoded
  }
}
```

### Anti-Pattern 4: Glob All Output Files

**What people do:**
```wdl
output {
  Array[File] all_outputs = glob("output_dir/**/*")
}
```

**Why it's wrong:**
- Captures log files, temp files, and other non-essential outputs
- Makes call-caching less effective (cache invalidated by irrelevant file changes)
- Downstream workflows don't know which files are which

**Do this instead:** Explicitly specify key output files:
```wdl
output {
  File virus_summary_tsv    = "output_dir/PREFIX_summary/PREFIX_virus_summary.tsv"
  File plasmid_summary_tsv  = "output_dir/PREFIX_summary/PREFIX_plasmid_summary.tsv"
  File virus_sequences_fna  = "output_dir/PREFIX_summary/PREFIX_virus.fna"
  # etc. - only files that downstream workflows actually need
}
```

## Integration Points

### External Services

| Service | Integration Pattern | Notes |
|---------|---------------------|-------|
| Terra/Cromwell | WDL workflow execution | Primary deployment target. Workflows must validate with `miniwdl check` and `womtool validate`. |
| DNAnexus | WDL workflow execution via dxWDL | Workflows published via .dockstore.yml. Use `dx_instance_type` in runtime block. |
| Docker Hub / GHCR | Container image registry | Genomad official images on Docker Hub. May need to mirror to GHCR for consistency. |
| GCS / S3 | Input/output file storage | Cromwell localizes files automatically. Use GCS URIs in Terra, local paths in miniwdl. |
| GitHub Actions | CI/CD testing | tests-miniwdl.sh runs on every PR. Must complete in <10 minutes (GHA timeout). |

### Internal Boundaries

| Boundary | Communication | Notes |
|----------|---------------|-------|
| tasks ↔ workflows | WDL import + call | Workflows import tasks via relative path: `import "../tasks/tasks_metagenomics.wdl"` |
| workflow ↔ test config | JSON input file | Test JSON provides input values. Convention: `test_inputs-{workflow_name}-local.json` |
| task ↔ Docker container | Docker runtime + command block | Task command block runs inside container. Container must have genomad installed. |
| WDL ↔ requirements-modules.txt | Version pinning | Docker image versions specified in requirements-modules.txt, synced to WDL via version-wdl-runtimes.sh |
| genomad ↔ viral-ngs tools | Sequential execution in command block | If BAM input, use viral-ngs `read_utils.py bam_to_fasta` before calling genomad |

## Sources

- [geNomad official documentation](https://portal.nersc.gov/genomad/)
- [geNomad GitHub repository](https://github.com/apcamargo/genomad)
- [geNomad pipeline documentation](https://portal.nersc.gov/genomad/pipeline.html)
- [geNomad Nature Biotechnology paper](https://www.nature.com/articles/s41587-023-01953-y)
- [nf-core genomad_endtoend module](https://nf-co.re/modules/genomad_endtoend)
- viral-pipelines codebase (tasks_metagenomics.wdl, classify_single.wdl, classify_multi.wdl)

---
*Architecture research for: genomad WDL pipeline integration*
*Researched: 2026-02-12*
