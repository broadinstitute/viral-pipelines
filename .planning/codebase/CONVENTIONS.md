# Coding Conventions

**Analysis Date:** 2026-03-31

## Language & Format

**Primary Language:** WDL (Workflow Description Language)

**Version:** WDL 1.0 (specified at top of all files: `version 1.0`)

**File Organization:** All files use consistent structure with version declaration, imports, then definitions.

## Naming Patterns

**Files:**
- Task files: `tasks_<category>.wdl` (e.g., `tasks_assembly.wdl`, `tasks_read_utils.wdl`, `tasks_intrahost.wdl`)
  - Located in `pipes/WDL/tasks/`
- Workflow files: `<workflow_name>.wdl` (e.g., `assemble_refbased.wdl`, `classify_kraken2.wdl`)
  - Located in `pipes/WDL/workflows/`
- Test input files: `test_inputs-<workflow_name>-local.json`
- Test output files: `test_outputs-<workflow_name>-local.json`

**Tasks:**
- Use snake_case for task names: `task align_reads { ... }`, `task merge_and_reheader_bams { ... }`
- Descriptive names indicating the action: `refine_assembly_with_aligned_reads`, `get_bam_samplename`

**Workflows:**
- Use snake_case for workflow names: `workflow assemble_refbased { ... }`, `workflow classify_kraken2 { ... }`

**Input Variables:**
- Use snake_case: `reads_unmapped_bams`, `reference_fasta`, `trim_clip_db`, `skip_mark_dupes`
- Boolean inputs use meaningful names: `always_succeed`, `skip_mark_dupes`, `compare_direct_neighbors`
- Optional parameters use `?`: `File?`, `String?`, `Array[File]?`

**Output Variables:**
- Use snake_case: `contigs_fasta`, `subsampBam`, `subsample_read_count`
- Some legacy inconsistency exists (e.g., `subsampBam` mixing camelCase with underscore prefix)
- Descriptive names indicating content: `aligned_only_reads_bam`, `primer_trimmed_read_count`

**Runtime Variables:**
- Disk sizes: `Int disk_size = <value>` (always specified before command block)
- Memory/CPU defaults: `Int machine_mem_gb`, `Int cpu` (optional inputs with defaults in runtime section)
- Docker image: `String docker = "image:version"` (specified as input with default value)

**Internal Variables:**
- Temporary files in commands use UPPERCASE: `SAMPLE_NAME`, `VERSION`, `UPTIME_SEC`, `MEM_BYTES`, `LOAD_15M`
- This convention separates temporary working values from persistent outputs

## Command Blocks & Scripts

**Bash Wrapper:**
- All command blocks begin with error handling: `set -e`, `set -ex`, or `set -ex -o pipefail`
  - `set -e`: Exit on first error
  - `set -x`: Print commands as executed
  - `set -o pipefail`: Fail if any command in pipe fails
- Example from `tasks_assembly.wdl`:
```wdl
command <<<
    set -ex -o pipefail

    # find 90% memory
    mem_in_mb=$(/opt/viral-ngs/scripts/calc_mem.py mb 90)
    mem_in_gb=$(/opt/viral-ngs/scripts/calc_mem.py gb 90)

    assembly --version | tee VERSION

    # command invocation
    assembly assemble_spades \
      ~{reads_unmapped_bam} \
      ...
>>>
```

**Inline Scripts:**
- Python scripts embedded with heredoc: `python3 << CODE` ... `CODE`
- Perl one-liners: `perl -lane 'if (/pattern/) { print ... }'`
- Shell parameter expansion via WDL: `~{variable_name}` and `~{sep="*" array_variable}`

**Logging & Version Tracking:**
- All tasks capture tool versions: `tool --version | tee VERSION`
- Output version as string: `String viralngs_version = read_string("VERSION")`
- Enable debug logging with `--loglevel DEBUG` flag to viral-ngs tools
- Resource monitoring: capture `UPTIME_SEC`, `LOAD_15M`, `MEM_BYTES` for performance analysis

## Input Organization

**Input Block Structure:**
1. Main input files/arrays first (e.g., `reads_unmapped_bams`, `reference_fasta`)
2. Optional files next (e.g., `File? trim_clip_db`)
3. Configuration parameters (integers, floats, strings for tool options)
4. Derived/computed defaults: `String sample_name = basename(...)`
5. Runtime parameters last (docker image, memory, cpu)

**Parameter Metadata:**
- Comprehensive `parameter_meta` block describing all inputs
- Each parameter includes:
  - `description`: What this parameter is
  - `patterns`: File patterns (if applicable) - e.g., `["*.bam"]`, `["*.fasta", "*.fasta.gz"]`
  - `category`: Classification like "required", "common", "other"
- Example from `tasks_assembly.wdl`:
```wdl
parameter_meta{
  reads_unmapped_bam: {
    description: "Unaligned reads in BAM format.",
    patterns: ["*.bam"],
    category: "required"
  }
  spades_options: {
    description: "Display additional options to pass the SPAdes assembler.",
    category: "other"
  }
}
```

**Metadata Block:**
- Tasks can include `meta` block with `description`:
```wdl
meta {
  description: "Merge and/or reheader bam files using a mapping table..."
}
```
- Workflows include `meta` block with author, email, and description:
```wdl
meta {
    description: "Reference-based microbial consensus calling..."
    author: "Broad Viral Genomics"
    email:  "viral-ngs@broadinstitute.org"
    allowNestedInputs: true
}
```

## Output Specification

**Output Block Organization:**
- Primary outputs declared first (assembly files, BAM files)
- Derived outputs next (counts, metrics)
- Version strings last
- All outputs explicitly typed (File, Array[File], String, Int, Map, etc.)

**Reading Output Values:**
- From files on disk: `File my_output = "~{basename}.fasta"`
- From stdout: `String my_output = read_string(stdout())`
- From text files: `Int count = read_int("filename")`
- From line-separated files: `Array[String] samples = read_lines("filename")`
- From TSV: `Array[Array[String]] table = read_tsv("filename.tsv")`
- From JSON: `Map[String,String] metadata = read_json("metadata.json")`
- Glob patterns: `Array[File] figures = glob("figs/*")`

## Runtime Block

**Standard Structure:**
```wdl
runtime {
    docker: docker
    memory: "~{select_first([machine_mem_gb, 32])} GB"
    cpu:    select_first([cpu, 8])
    disks: "local-disk ~{disk_size} SSD"
    disk: "~{disk_size} GB" # TES (Task Execution Service)
    dx_instance_type: "mem1_ssd1_v2_x8" # DNAnexus-specific
    maxRetries: 2
}
```

**Docker Images:**
- Pinned to specific versions: `ghcr.io/broadinstitute/viral-ngs:3.0.6-core`
- Flavor-based images for specialized tools: `-assemble`, `-phylo`, `-core`, `-classify`
- Version consistency enforced via `requirements-modules.txt` and CI validation
- Alternative registries: `quay.io/broadinstitute/`, `python:slim`, `nextstrain/nextclade:3.18.1`

**Memory & CPU:**
- Use `select_first([variable, default])` for optional resource overrides
- Common defaults: 32 GB memory, 8 CPU for compute tasks
- Lighter tasks: 1 GB memory, 1 CPU
- Always allow optional overrides via `Int? machine_mem_gb`, `Int? cpu`

**Disk:**
- Both `disks` (WDL standard) and `disk` (TES) specified
- TES format for compatibility: `disk: "~{disk_size} GB"`
- Disk types: SSD (faster) or HDD (cheaper)
- Disk size calculated based on input size or fixed estimate

**Platform-Specific Configuration:**
- DNAnexus `dx_instance_type`: maps to machine types (e.g., `mem1_ssd1_v2_x8`)
- Always include for DNAnexus deployment compatibility

## Call Invocation Pattern

**Standard Pattern in Workflows:**
```wdl
call task_name as alias_name {
    input:
        param1 = value1,
        param2 = value2
}
```

**Key Conventions:**
- Use `as alias_name` when calling same task multiple times: `call assembly.align_reads as align_to_ref`
- Alias describes the specific usage context
- Multi-line call inputs with one parameter per line
- No commas after last parameter

**Conditional Execution:**
```wdl
if(defined(filter_to_taxon_db)) {
    call taxon_filter.filter_to_taxon {
        input:
            reads_unmapped_bam = reads_depleted_bams,
            lastal_db_fasta    = select_first([filter_to_taxon_db])
    }
}
```

**Scatter Pattern:**
```wdl
scatter(reads_unmapped_bam in reads_unmapped_bams) {
    call assembly.align_reads {
        input:
            reads_unmapped_bam = reads_unmapped_bam
    }
}
```

## Error Handling

**Command-Level Error Handling:**
- Use `set -e` (fail on error) as default
- Use `set -o pipefail` when piping commands
- Pattern: `set -ex -o pipefail` for maximum safety

**Graceful Failures:**
- Some tasks support graceful degradation via boolean flags
- Example from `tasks_assembly.wdl`: `always_succeed = false` allows SPAdes failure to output empty FASTA instead of failing the task
- Use `if [ $? -eq 0 ]` for explicit exit code checking when needed

**Conditional Output Selection:**
- Use `select_first([option1, option2])` to pick first available result
- Example: `File reads_taxfilt_bams = select_first([filter_to_taxon.taxfilt_bam, reads_depleted_bams])`
- This allows optional processing steps with fallback to unprocessed input

## Imports & Composition

**Import Syntax:**
```wdl
import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_read_utils.wdl" as read_utils
```

**Conventions:**
- Use relative paths with `../` to navigate from workflows to tasks
- Namespace with short alias: `as assembly`, `as read_utils`, `as taxon_filter`
- Call tasks via alias: `call assembly.align_reads`
- Imports appear at top of file before workflow/task definition

**Workflow Composition:**
- Workflows compose tasks into pipelines
- No direct code reuse between workflows (only task composition)
- Workflows can import and call other workflows (example: `assemble_denovo.wdl` imports `assemble_refbased.wdl`)

## Comments

**In-Command Documentation:**
- Bash comments with `#` for explanation of operations
- Example from `tasks_read_utils.wdl`:
```bash
# filename must be <samplename>.l<xxx>.bam
assert bam.endswith('.bam'), "filename does not end in .bam: {}".format(bam)
```

**Metadata Over Comments:**
- Prefer `parameter_meta` and `meta` blocks over inline comments
- Parameter descriptions should be thorough and explain both what and why

**Complex Logic:**
- Inline Python/Perl scripts include explanatory comments
- Example pattern: `# WDL arrays to python arrays` before conversion logic

## Special Conventions

**Basename Derivation:**
- Use `basename()` to extract sample names from BAM files
- Pattern: `String sample_name = basename(basename(reads_unmapped_bam, ".bam"), ".taxfilt")`
- Double `basename()` to strip multiple extensions

**String Interpolation:**
- WDL variable expansion: `"~{variable}"`
- Array joining: `~{sep=' ' array_variable}` (space-separated)
- Alternative separators: `~{sep='*' list}` for parsing in Python

**Resource Calculation:**
- Memory calculated as percentage of available: `/opt/viral-ngs/scripts/calc_mem.py gb 90`
- Enables tasks to use 90% of provisioned memory safely

**Skip-Global-Version-Pin Comment:**
- Comments `#skip-global-version-pin` on version pinning lines allow exceptions
- Used when specific version overrides are intentional

---

*Conventions analysis: 2026-03-31*
