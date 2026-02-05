# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

viral-pipelines is a collection of WDL (Workflow Description Language) workflows and tasks for analyzing viral NGS data. The repository contains workflows and tasks that can be executed on various platforms including local machines (miniWDL/Cromwell), cloud platforms (Terra, DNAnexus, AWS, GCP, Azure), and HPC systems.

## Repository Structure

- `pipes/WDL/workflows/` - Main WDL workflow definitions
- `pipes/WDL/tasks/` - Reusable WDL task modules organized by function:
  - `tasks_assembly.wdl` - De novo and reference-based assembly tasks
  - `tasks_read_utils.wdl` - BAM/FASTQ processing utilities
  - `tasks_taxon_filter.wdl` - Taxonomic classification and filtering
  - `tasks_intrahost.wdl` - Variant calling and iSNV analysis
  - `tasks_interhost.wdl` - Phylogenetics and multi-sample analysis
  - `tasks_ncbi.wdl`, `tasks_ncbi_tools.wdl` - GenBank/SRA submission
  - `tasks_sarscov2.wdl` - SARS-CoV-2 specific tools
  - `tasks_nextstrain.wdl` - Nextstrain/Augur integration
  - `tasks_reports.wdl`, `tasks_utils.wdl`, `tasks_terra.wdl`, etc.
- `test/input/` - Test data and input JSON files
- `github_actions_ci/` - CI/CD scripts for validation and testing
- `docs/` - Sphinx documentation (published to ReadTheDocs)
- `requirements-modules.txt` - Docker image versions for all dependencies
- `.dockstore.yml` - Dockstore registry configuration

## Development Commands

### WDL Validation

```bash
# Validate all workflows with miniwdl
miniwdl check pipes/WDL/workflows/*.wdl

# Validate with womtool (requires Java 11 and Cromwell installed)
github_actions_ci/install-wdl.sh  # First-time setup
github_actions_ci/validate-wdl-womtool.sh
```

### Testing Workflows

```bash
# Test workflows locally with miniwdl
# Test inputs are in test/input/WDL/miniwdl-local/
miniwdl run pipes/WDL/workflows/assemble_refbased.wdl \
  -i test/input/WDL/miniwdl-local/test_inputs-assemble_refbased-local.json

# Run all miniwdl tests (as done in CI)
github_actions_ci/tests-miniwdl.sh

# Run all Cromwell tests (as done in CI)
github_actions_ci/install-wdl.sh  # Install Cromwell first
github_actions_ci/tests-cromwell.sh
```

### Testing Individual Tasks

```bash
# Test a single task in isolation (useful for debugging)
miniwdl run --task task_name pipes/WDL/tasks/tasks_assembly.wdl \
  input1=value1 \
  input2=value2
```

### List Workflow Inputs

```bash
# List required and optional inputs for any workflow
miniwdl run pipes/WDL/workflows/assemble_refbased.wdl
```

### Documentation

```bash
# Build documentation locally
github_actions_ci/install-pip-docs.sh  # First-time setup
github_actions_ci/build-docs.sh
```

### Docker Image Management

The `requirements-modules.txt` file specifies exact Docker image versions for all dependencies. When updating task files:

```bash
# Check if WDL runtime Docker versions match requirements-modules.txt
github_actions_ci/check-wdl-runtimes.sh

# Update versions in WDL files to match requirements-modules.txt
github_actions_ci/version-wdl-runtimes.sh
```

## Development Guidelines for Claude

### Syntax Validation
- Use `miniwdl check` for WDL syntax validation
- Run validation after any WDL edits before committing

### Testing Requirements
- Use `miniwdl run` for integration testing workflows
- Use `miniwdl run --task` to test individual tasks in isolation when debugging
- After editing any WDL task, verify all existing integration tests for workflows that depend on that task
- Any newly created tasks or workflows must include short/basic integration test cases in the existing testing framework:
  - Add test input JSON to `test/input/WDL/miniwdl-local/`
  - Follow the naming convention: `test_inputs-{workflow_name}-local.json`
  - Optionally add expected outputs: `test_outputs-{workflow_name}-local.json`

### Workflow Testing After Task Changes
When editing a task file (e.g., `tasks_assembly.wdl`):
1. Identify workflows that import the edited task
2. Run integration tests for those workflows with `miniwdl run`
3. Verify outputs match expected results

## WDL Architecture

### Workflow Organization

Workflows import tasks from the `tasks/` modules and compose them into pipelines. Common workflow patterns:

1. **Single-sample assembly**: `assemble_refbased.wdl`, `assemble_denovo.wdl`
2. **Classification**: `classify_kraken2.wdl`, `classify_kaiju.wdl`
3. **Phylogenetics**: `augur_from_assemblies.wdl`, `mafft_and_snp.wdl`
4. **SARS-CoV-2 pipelines**: `sarscov2_illumina_full.wdl`, `sarscov2_nextstrain.wdl`
5. **Data submission**: `submit_genbank.wdl`, `submit_sra.wdl`

### Task Structure

Tasks in `tasks/*.wdl` files define:
- Docker images (specified in `docker` input parameter)
- Resource requirements (`machine_mem_gb`, `cpu`, `disk`)
- Command execution (in `command` blocks)
- Input/output parameters with metadata

### Import Paths

WDL imports use relative paths from workflows to tasks:
```wdl
import "../tasks/tasks_assembly.wdl" as assembly
```

## Testing

Test configurations are in `test/input/WDL/`:
- `test_inputs-{workflow_name}-local.json` - Cromwell test inputs
- `miniwdl-local/test_inputs-{workflow_name}-local.json` - miniWDL test inputs
- `test_outputs-{workflow_name}-local.json` - Expected outputs for validation

Not all workflows have test configurations; only those with corresponding JSON files in the test directory are tested in CI.

## CI/CD

GitHub Actions (`.github/workflows/build.yml`) runs on all PRs and pushes:
1. `validate_wdl_miniwdl` - Validates WDL syntax with miniwdl
2. `validate_wdl_womtool` - Validates WDL with Cromwell's womtool
3. `test_docs` - Builds Sphinx documentation
4. `test_miniwdl` - Runs workflows with test inputs using miniwdl
5. `test_cromwell` - Runs workflows with test inputs using Cromwell
6. `deploy_dnanexus` - Deploys to DNAnexus (on master branch or releases)

## Key Workflows

- **assemble_refbased.wdl**: Reference-based consensus calling from BAM files
  - Aligns reads to reference, trims primers (optional), calls consensus
  - Supports novoalign, bwa, or minimap2 aligners
  - Primary workflow for viral genome assembly

- **assemble_denovo_metagenomic.wdl**: De novo metagenomic assembly with SPAdes

- **classify_single.wdl**: Taxonomic classification and depletion pipeline

- **nextclade_single.wdl**: Nextclade analysis for single samples

- **genbank_single.wdl**: GenBank submission preparation for single samples

- **augur_from_assemblies.wdl**: Nextstrain phylogenetic analysis from assemblies

## Docker Images

All compute tasks run in Docker containers. The primary runtime images are multiflavored builds of `viral-ngs` accessed from `ghcr.io/broadinstitute/viral-ngs:{version}-{flavor}` (also mirrored to `quay.io`). Dockerfiles and build definitions live in `github.com/broadinstitute/viral-ngs`. Other images (ncbi-tools, py3-bio, etc.) are accessed from `quay.io/broadinstitute/{image-name}`.

Flavored images include:
- `viral-ngs:{version}-core` - Core bioinformatics tools
- `viral-ngs:{version}-assemble` - Assembly tools (SPAdes, etc.)
- `viral-ngs:{version}-classify` - Classification tools (Kraken2, etc.)
- `viral-ngs:{version}-phylo` - Phylogenetics tools (Augur, etc.)
- `ncbi-tools` - NCBI submission tools (separate image on quay.io)

Image versions are pinned in `requirements-modules.txt` and must be kept in sync with WDL files.

## Dockstore Integration

Workflows are registered on Dockstore for easy import to Terra, DNAnexus, and other platforms. The `.dockstore.yml` file defines all published workflows and their test parameter files.

## Terra Performance Analysis

When analyzing workflow performance from Terra submissions, use the Terra MCP tools for structure/status queries and direct GCS access for log analysis.

### Timing Methodology for WDL Tasks

**Preferred method - use `get_batch_job_status`:**

The Terra MCP's `get_batch_job_status` tool returns timing data directly from the Google Batch API:

```
get_batch_job_status(
  workspace_namespace="<namespace>",
  workspace_name="<workspace>",
  submission_id="<submission-uuid>",
  workflow_id="<workflow-uuid>",
  task_name="<task_name>",
  shard_index=<optional>,
  attempt=<optional>
)
```

Returns timing in the `batch_job.timing` field:
- **run_duration**: Actual task execution time (what you usually want for performance analysis)
- **pre_run_duration**: Queue and setup time (VM provisioning, Docker pull, etc.)

This is more accurate than log-based methods because it captures the complete execution including post-script I/O operations.

**Alternative method - log-based timing (for detailed analysis):**

When you need finer-grained timing within a task (e.g., timing individual steps):

1. **Start time**: Use first Python log timestamp in stderr
   - Pattern: `^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}),\d+`

2. **End time**: Use GCS file modification timestamp of stderr
   - Get via: `gcloud storage ls -l <path>/stderr`
   - This captures ALL execution including post-Python BAM I/O

3. **Why not use Python log end time?**
   - Many tasks run external tools (Java, pysam) after Python logging ends
   - Python logs don't capture full execution time

### Efficient GCS Queries with Wildcards

**Always use `gcloud storage` instead of `gsutil`** - it's faster, more reliable, and the preferred CLI for GCS operations.

Use wildcards to batch GCS queries instead of iterating:
```bash
# Get all stderr files from a submission with timestamps in one query
gcloud storage ls -l "gs://bucket/submissions/<sub_id>/classify_single/*/call-deplete/stderr"
gcloud storage ls -l "gs://bucket/submissions/<sub_id>/classify_single/*/call-deplete/attempt-*/stderr"
```

### Handling Preemption Retries

When a task is preempted, Cromwell creates `attempt-*` directories:
```
call-deplete/
  stderr           # First attempt (may be incomplete)
  attempt-2/       # Second attempt
    stderr         # Final successful run
```

**Always use the final (highest-numbered) attempt** for performance analysis - preemption time shouldn't count against code performance.

### Sample Identification

To identify which workflow corresponds to which sample:
1. Read first few KB of stderr from each workflow
2. Look for sample name in BAM file paths (e.g., `/S20.l1.xxxx.bam`)
3. Cache the sample-to-workflow mapping for reuse

### Debugging Infrastructure-Level Failures

Some workflow failures have errors that aren't visible in standard stderr logs. These include:
- Docker pull failures (rate limits, image not found, auth errors)
- VM provisioning failures
- Preemption before task execution started
- Network connectivity issues during container setup

**Signs you need Batch logs instead of stderr:**
- Batch reports exit code 0 (success) but task is marked as failed ("GCP Batch task exited with Success(0)")
- Error message says "The job was stopped before the command finished"
- stderr is empty or very short
- Error message says "Executor error" without details
- Task failed instantly (0 seconds runtime)
- `get_job_metadata` summary shows failure but no useful error message

**Use `get_batch_job_status` to diagnose infrastructure failures:**

The Terra MCP provides `get_batch_job_status` which queries the Google Batch API directly:

```
get_batch_job_status(
  workspace_namespace="<namespace>",
  workspace_name="<workspace>",
  submission_id="<submission-uuid>",
  workflow_id="<workflow-uuid>",
  task_name="<task_name>",
  shard_index=<optional>,  # For scattered tasks
  attempt=<optional>       # For retried tasks
)
```

The tool returns:
- **Batch job status**: QUEUED, SCHEDULED, RUNNING, SUCCEEDED, or FAILED
- **Timing**: run_duration and pre_run_duration (queue/setup time)
- **Resources**: machine_type, CPU, memory, disk sizes
- **Status events**: State transitions with timestamps
- **Detected issues**: Auto-detected problems with severity and suggestions
- **Cloud Logging query**: Ready-to-use gcloud command for deeper debugging

**Recommended debugging workflow:**
1. `get_submission_status` → identify failed workflows
2. `get_job_metadata` (summary mode) → identify failed tasks and error messages
3. `get_workflow_logs` → check stderr for application errors
4. `get_batch_job_status` → check infrastructure issues if logs don't explain failure

**Common failure patterns detected:**
- `"Failed to pull image"` - Check image name, tag, and registry auth
- `"429 Too Many Requests"` - Registry rate limit, retry later
- `"manifest unknown"` - Image tag doesn't exist
- `"unauthorized"` - Service account lacks permission to pull from registry
- `"PREEMPTED"` - VM was preempted, usually retried automatically
- `"exit code 137"` - OOM killed (out of memory)
- `"exit code 1"` - Application error in the task script
