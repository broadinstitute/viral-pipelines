# Coding Conventions

**Analysis Date:** 2026-02-11

## Naming Patterns

**WDL Files:**
- Workflows: `{purpose}.wdl` (e.g., `assemble_refbased.wdl`, `classify_kraken2.wdl`)
- Tasks: `tasks_{category}.wdl` (e.g., `tasks_assembly.wdl`, `tasks_read_utils.wdl`)
- Naming uses snake_case throughout

**WDL Tasks:**
- Task names: snake_case (e.g., `align_reads`, `merge_and_reheader_bams`, `ivar_trim`)
- Input parameters: snake_case (e.g., `reads_unmapped_bam`, `reference_fasta`, `trim_coords_bed`)
- Output variables: snake_case (e.g., `aligned_only_reads_bam`, `contigs_fasta`, `subsample_read_count`)

**WDL Variables:**
- Sample names derived from file basenames: `basename(basename(file, ".bam"), ".taxfilt")`
- Default values use parameter references: `sample_name = basename(reads_unmapped_bams[0], '.bam')`

**Bash Scripts:**
- Script names: kebab-case (e.g., `check-wdl-runtimes.sh`, `tests-miniwdl.sh`, `validate-wdl-womtool.sh`)
- Shell variables: UPPER_SNAKE_CASE for environment variables, lower_snake_case for local variables

**Python (inline in WDL):**
- Variables: snake_case (e.g., `bam_uris`, `sample_to_bams`, `meta_cols`)
- Module imports: Standard library style (e.g., `import os.path`, `import csv`, `import json`)

## Code Style

**WDL Version:**
- All WDL files use `version 1.0` as the first line

**WDL Formatting:**
- Indentation: 2 or 4 spaces (inconsistent across files)
- Command blocks use heredoc syntax: `command <<<` ... `>>>`
- Runtime blocks use consistent structure with standard keys
- Parameter blocks aligned with consistent spacing

**WDL Docker Images:**
- Always specified as input parameter: `String docker = "ghcr.io/broadinstitute/viral-ngs:3.0.4-{flavor}"`
- Versions must match `requirements-modules.txt`
- Flavored images for different purposes: `-core`, `-assemble`, `-classify`, `-phylo`
- Optional version pin skip via comment: `#skip-global-version-pin`

**Bash in WDL:**
- Always start with strict error handling: `set -ex -o pipefail` (or `set -e` for simpler scripts)
- Exception: `set -ex` without pipefail when using grep (which exits 1 on no match)
- Disable pipefail explicitly when needed: `set +o pipefail`

**Resource Metrics:**
- Tasks capture runtime metrics: `UPTIME_SEC`, `LOAD_15M`, `MEM_BYTES`
- Memory calculation: `/opt/viral-ngs/scripts/calc_mem.py mb 90` (find 90% of available memory)
- Conditional memory detection for cgroups v1 vs v2:
  ```bash
  { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak;
    elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak;
    elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes;
    else echo "0"; fi; } > MEM_BYTES
  ```

## Import Organization

**WDL Imports:**
- Tasks imported with relative paths: `import "../tasks/tasks_assembly.wdl" as assembly`
- Workflows can import other workflows: `import "assemble_refbased.wdl" as assemble_refbased`
- Aliases use category names: `as assembly`, `as reports`, `as read_utils`
- Standard import order: tasks first, then workflows

**Python Imports (in WDL command blocks):**
- Standard library first
- Third-party packages after
- Local modules last (e.g., `from viral_ngs.core import file as util_file`)

## Error Handling

**WDL Task-Level:**
- Optional `always_succeed` parameter in some tasks to emit empty output instead of failing
- Validation via assertions in Python blocks: `assert bam.endswith('.bam'), "filename does not end in .bam: {}".format(bam)`
- `maxRetries: 2` in runtime blocks for transient failures

**Bash Scripts:**
- Exit immediately on error: `set -e`
- Exit on undefined variables (when using `-u`)
- Pipefail to catch errors in pipes: `set -o pipefail`
- Trap cleanup on exit: `trap cleanup EXIT SIGINT SIGQUIT SIGTERM`

**Error Output:**
- Tools run with `--loglevel=DEBUG` or `--verbose` flags
- Version logging: `assembly --version | tee VERSION`
- Stderr redirection for counts: `samtools view -c file.bam | tee subsample_read_count >&2`

## Logging

**Framework:** Inline in bash commands, no structured logging framework

**Patterns:**
- Capture tool versions: `tool --version | tee VERSION`
- Echo progress: `echo "Executing $workflow_name using miniWDL on local instance"`
- Use `tee` to capture output to files and display to console
- Date stamps in test scripts: `date` before and after operations
- Redirect counts to stderr: `| tee count >&2`

## Comments

**WDL Metadata:**
- `meta` blocks describe task/workflow purpose
- `parameter_meta` blocks document each input/output:
  ```wdl
  parameter_meta {
    input_name: {
      description: "Human-readable description"
      patterns: ["*.bam"]
      category: "required" | "common" | "other"
    }
  }
  ```

**Inline Comments:**
- Explain non-obvious logic: `# do this in two steps in case the input doesn't actually have "taxfilt" in the name`
- Document WDL features: `# https://wdl-aid.readthedocs.io/en/latest/usage.html`
- Explain bash patterns: `# find 90% memory`
- Skip rules: `#skip-global-version-pin`

**Comment Style:**
- Bash: `# comment`
- WDL: `# comment` (not `//`)
- Python: `# comment`

## Function Design

**WDL Tasks:**
- Single responsibility per task
- Clear input/output contracts via parameter_meta
- Resource requirements specified: `machine_mem_gb`, `cpu`, `disk_size`
- Platform-specific hints: `dx_instance_type` for DNAnexus

**WDL Workflows:**
- Compose tasks via scatter/gather patterns: `scatter(reads_unmapped_bam in reads_unmapped_bams) { ... }`
- Use Maps for configuration: `Map[String,String] align_to_ref_options`
- Default values in input block, not parameter declarations

**Bash Functions:**
- Use functions for cleanup: `function cleanup() { ... }`
- Absolute path resolution: `function absolute_path() { ... }`

## Module Design

**WDL Task Modules:**
- Organized by functional area:
  - `tasks_assembly.wdl` - Assembly operations
  - `tasks_read_utils.wdl` - BAM/FASTQ utilities
  - `tasks_taxon_filter.wdl` - Classification
  - `tasks_intrahost.wdl` - Variant calling
  - `tasks_interhost.wdl` - Phylogenetics
  - `tasks_reports.wdl` - Reporting and QC
  - `tasks_utils.wdl` - General utilities
  - `tasks_ncbi.wdl`, `tasks_ncbi_tools.wdl` - NCBI submission
  - `tasks_sarscov2.wdl` - SARS-CoV-2 specific tools
  - `tasks_nextstrain.wdl` - Nextstrain integration

**Exports:**
- Tasks are imported and called with module namespace: `assembly.align_reads`, `reports.MultiQC`

**Workflow Composition:**
- Workflows can import other workflows as reusable components
- Nested workflows maintain their own namespace

## Version Management

**Docker Version Pinning:**
- All versions defined in `requirements-modules.txt` using format: `owner/image=version`
- Script validates versions: `github_actions_ci/check-wdl-runtimes.sh`
- Script updates versions: `github_actions_ci/version-wdl-runtimes.sh`
- Flavored variants allowed: `module:version-flavor` where flavor is `-core`, `-assemble`, etc.

**WDL Version:**
- All files use WDL version 1.0: `version 1.0`

## Platform Compatibility

**Runtime Specifications:**
- Docker/Cromwell: `docker`, `memory`, `cpu`, `disks`
- TES: `disk` (in addition to `disks`)
- DNAnexus: `dx_instance_type`
- Example runtime block:
  ```wdl
  runtime {
    docker: docker
    memory: "~{select_first([machine_mem_gb, 32])} GB"
    cpu: select_first([cpu, 8])
    disks: "local-disk ~{disk_size} SSD"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x8"
    maxRetries: 2
  }
  ```

**Validation:**
- Validate with both miniwdl and womtool
- Maintain compatibility across Cromwell, miniWDL, Terra, DNAnexus, AWS

---

*Convention analysis: 2026-02-11*
