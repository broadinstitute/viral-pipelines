# Architecture

**Analysis Date:** 2026-02-11

## Pattern Overview

**Overall:** WDL (Workflow Description Language) pipeline composition pattern

**Key Characteristics:**
- Declarative workflow definition using WDL 1.0 specification
- Task-based modular architecture with workflows composing reusable tasks
- Container-based execution with all computation in Docker images
- Platform-agnostic design supporting local (miniWDL), cloud (Cromwell/Terra), and HPC execution
- Data-flow oriented with explicit file inputs/outputs and type safety

## Layers

**Workflows Layer:**
- Purpose: Define end-to-end analysis pipelines for viral genomics use cases
- Location: `pipes/WDL/workflows/`
- Contains: 90+ workflow definitions orchestrating multiple tasks
- Depends on: Task modules (via relative imports)
- Used by: End users via miniWDL, Cromwell, Terra, DNAnexus, etc.

**Tasks Layer:**
- Purpose: Encapsulate individual bioinformatics operations as reusable units
- Location: `pipes/WDL/tasks/`
- Contains: 16 task module files organized by domain (assembly, taxonomy, variants, etc.)
- Depends on: Docker runtime images (viral-ngs, ncbi-tools, etc.)
- Used by: Workflow layer via namespaced imports

**Runtime Layer:**
- Purpose: Provide bioinformatics software environment and execution isolation
- Location: Docker images specified in `requirements-modules.txt`
- Contains: Pre-built Docker images from ghcr.io/broadinstitute and quay.io registries
- Depends on: External container registries
- Used by: Task `runtime` blocks via `docker` parameter

**Configuration Layer:**
- Purpose: Parameterize workflows for specific analyses
- Location: `test/input/WDL/` (test configs), user-provided JSON (production)
- Contains: JSON input files specifying workflow parameters
- Depends on: Workflow input schema
- Used by: Workflow execution engines

## Data Flow

**Reference-Based Assembly Flow (assemble_refbased.wdl):**

1. **Input Stage**: Receive unmapped BAM files and reference genome
2. **Parallel Alignment**: Scatter across input BAMs, align to reference (novoalign/bwa/minimap2)
3. **Primer Trimming**: Apply ivar_trim if BED coordinates provided
4. **Merging**: Merge aligned BAMs if multiple inputs
5. **Variant Calling**: Call variants against reference with lofreq
6. **Consensus Calling**: Generate consensus sequence from aligned reads
7. **Self-Alignment**: Align original reads to new consensus
8. **QC/Metrics**: Generate coverage plots, alignment metrics, variant calls
9. **Output Stage**: Emit consensus FASTA, VCF, BAM, metrics, plots

**State Management:**
- Stateless task execution - all state in files
- Workflow engines (Cromwell/miniWDL) manage task scheduling and intermediate files
- No shared mutable state between tasks
- File localization and delocalization handled by execution backend

## Key Abstractions

**Workflow:**
- Purpose: Compose tasks into analysis pipeline with control flow
- Examples: `pipes/WDL/workflows/assemble_refbased.wdl`, `pipes/WDL/workflows/classify_single.wdl`
- Pattern: Import task modules, define inputs, use scatter/if for parallelism/conditionals, call tasks, define outputs

**Task:**
- Purpose: Encapsulate single bioinformatics operation with inputs, command, outputs, resource requirements
- Examples: `tasks_assembly.wdl::assemble`, `tasks_assembly.wdl::align_reads`, `tasks_intrahost.wdl::lofreq`
- Pattern: Define typed inputs/outputs, specify Docker image, write bash command block, declare runtime resources

**Import with Namespace:**
- Purpose: Organize and scope task references
- Examples: `import "../tasks/tasks_assembly.wdl" as assembly`
- Pattern: Relative path imports from workflows to tasks, namespace prevents name collisions

**Scatter-Gather:**
- Purpose: Parallelize operations across collections
- Examples: Aligning multiple input BAMs in parallel in `assemble_refbased.wdl` (line 87)
- Pattern: `scatter(item in array) { call task { input: x = item } }`

**Conditional Execution:**
- Purpose: Execute tasks only when conditions met
- Examples: Merge BAMs only if multiple inputs (line 110 in assemble_refbased.wdl)
- Pattern: `if(condition) { call task }` with `select_first` to handle optional outputs

## Entry Points

**miniWDL Execution:**
- Location: Command line via `miniwdl run <workflow.wdl>`
- Triggers: User invocation on local machine or shared server
- Responsibilities: Validates WDL, localizes inputs, executes tasks in Docker, manages outputs

**Cromwell Execution:**
- Location: Cromwell server (local, HPC, or cloud)
- Triggers: Workflow submission via API or command line
- Responsibilities: Schedules tasks on backend (GCP, AWS, Azure, HPC), handles retries, manages storage

**Terra Platform:**
- Location: https://app.terra.bio
- Triggers: User launches workflow in workspace via Dockstore TRS import
- Responsibilities: Provides UI, data tables, workspace storage, executes via Cromwell on GCP

**DNAnexus Platform:**
- Location: DNAnexus project (continuous deployment to CI project)
- Triggers: Workflow compiled from WDL to DNAnexus applets via dxWDL
- Responsibilities: Translates WDL to platform-native format, executes on DNAnexus infrastructure

**Dockstore Registry:**
- Location: https://dockstore.org/organizations/BroadInstitute/collections/pgs
- Triggers: Git tags/releases trigger automatic sync via `.dockstore.yml`
- Responsibilities: TRS API endpoint for workflow discovery and import to platforms

## Error Handling

**Strategy:** Fail-fast with configurable retries at task level

**Patterns:**
- Tasks exit with non-zero code on failure (bash `set -ex`)
- Workflow engines detect failure via exit code
- Retries configured in task `runtime` block (`maxRetries: 2`)
- Optional graceful degradation via `always_succeed` input parameters (e.g., assembly task)
- Conditional outputs with `select_first` handle optional task execution

## Cross-Cutting Concerns

**Logging:**
- Each task writes stdout/stderr to execution directory
- Python-based tasks use viral-ngs logging framework (timestamps, levels)
- Workflow engines capture and store logs (GCS for Terra, local dirs for miniWDL)

**Validation:**
- WDL syntax validation via `miniwdl check` and `womtool validate`
- Type checking enforced by WDL spec
- Integration tests verify end-to-end execution with test inputs
- CI/CD runs validation and tests on all PRs

**Resource Management:**
- Each task declares memory, CPU, disk requirements in `runtime` block
- Platform-specific fields support multiple backends (disks for GCP, disk for TES, dx_instance_type for DNAnexus)
- Tasks measure actual resource usage (MEM_BYTES, runtime_sec, cpu_load_15min)

**Version Control:**
- Docker image versions pinned in `requirements-modules.txt`
- Workflows import tasks via relative paths (executed versions travel together)
- Task outputs include version strings (viralngs_version, ivar_version, etc.)
- Git tags correlate with Dockstore releases

**Parallelism:**
- Scatter-gather pattern for data parallelism (multiple samples, multiple input files)
- WDL execution engines determine task concurrency based on dependencies
- Independent tasks can execute in parallel automatically

**Data Types:**
- WDL primitive types: String, Int, Float, Boolean, File
- Composite types: Array[T], Map[K,V], Pair[L,R], Struct
- Optional types: T? for nullable values
- Type coercion and validation enforced at workflow parse time

---

*Architecture analysis: 2026-02-11*
