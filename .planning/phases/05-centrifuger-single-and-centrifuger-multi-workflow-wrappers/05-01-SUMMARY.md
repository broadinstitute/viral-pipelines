---
phase: 05-centrifuger-single-and-centrifuger-multi-workflow-wrappers
plan: 01
subsystem: wdl
tags: [wdl, centrifuger, metagenomics, classification, bash-loop, array-inputs]

# Dependency graph
requires:
  - phase: 04-core-centrifuger-task
    provides: "centrifuger WDL task (scalar BAM input, single-sample)"
provides:
  - "Modified centrifuger task with Array[File] reads_bams and krakenuniq-style bash loop"
  - "centrifuger_single.wdl workflow (scalar BAM wrapped into 1-element array)"
  - "centrifuger_multi.wdl workflow (array of BAMs, no scatter, DB loaded once)"
affects: [06-centrifuger-test-inputs-and-dockstore]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Array[File] + bash for-loop inside task command block (krakenuniq pattern) for amortizing large DB load"
    - "1-element array literal [reads_bam] to wrap scalar into array-accepting task"
    - "Array[File]+ non-empty array type in workflow inputs for batch workflows"
    - "glob() outputs (classification_tsvs, kreports, centrifuger_logs) for array task results"

key-files:
  created:
    - pipes/WDL/workflows/centrifuger_single.wdl
    - pipes/WDL/workflows/centrifuger_multi.wdl
  modified:
    - pipes/WDL/tasks/tasks_metagenomics.wdl

key-decisions:
  - "ONE task (centrifuger), TWO workflow wrappers (single/multi) — same krakenuniq approach"
  - "centrifuger_multi.wdl has NO scatter block — loop is inside task command block"
  - "centrifuger_multi defaults: machine_mem_gb=480, cpu=16 for NT-scale batch processing"
  - "Sample names derived from basename(bam, '.bam') inside bash loop — no sample_name input"
  - "FASTQ inputs dropped — BAM-only interface for array mode; picard SamToFastq converts per-sample"

patterns-established:
  - "Workflow call alias pattern: call metagenomics.centrifuger as run_centrifuger (avoids name collision)"
  - "Both workflows: meta { allowNestedInputs: true } for Terra-compatible JSON inputs"
  - "Bash loop FASTQ cleanup: rm -f per-sample FQ files after each iteration to free disk"

requirements-completed: [CFGR-03, CFGR-04]

# Metrics
duration: 3min
completed: 2026-04-02
---

# Phase 05 Plan 01: centrifuger_single and centrifuger_multi Workflow Wrappers Summary

**centrifuger task refactored to Array[File] reads_bams with krakenuniq bash loop; two thin workflow wrappers (centrifuger_single.wdl, centrifuger_multi.wdl) wrapping the same task, both passing miniwdl check**

## Performance

- **Duration:** ~3 min
- **Started:** 2026-04-02T14:48:10Z
- **Completed:** 2026-04-02T14:50:48Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- Modified `centrifuger` task in `tasks_metagenomics.wdl`: scalar `File? reads_bam` + FASTQ inputs replaced with `Array[File] reads_bams`; krakenuniq-style bash loop processes each BAM with picard SamToFastq + centrifuger + centrifuger-kreport; outputs changed to Array[File] via glob
- Created `centrifuger_single.wdl`: thin wrapper that accepts scalar `File reads_bam` and passes `[reads_bam]` (1-element array literal) to the task — Terra scatter-friendly
- Created `centrifuger_multi.wdl`: batch wrapper with `Array[File]+ reads_bams`, no scatter block, machine_mem_gb=480/cpu=16 for amortizing 200+ GB database load across all samples on one node

## Task Commits

Each task was committed atomically:

1. **Task 1: Convert centrifuger task from scalar to array inputs with bash loop** - `537dac67` (feat)
2. **Task 2: Create centrifuger_single.wdl and centrifuger_multi.wdl workflow wrappers** - `a02d75ff` (feat)

**Plan metadata:** (docs commit follows)

## Files Created/Modified
- `pipes/WDL/tasks/tasks_metagenomics.wdl` - centrifuger task: Array[File] reads_bams input, bash loop, glob outputs, updated disk formula and parameter_meta
- `pipes/WDL/workflows/centrifuger_single.wdl` - Single-sample workflow wrapper (scalar BAM → 1-element array call)
- `pipes/WDL/workflows/centrifuger_multi.wdl` - Batch workflow wrapper (Array[File]+ BAMs, no scatter, 480 GB / 16 CPU defaults)

## Decisions Made
- Removed FASTQ inputs (`reads_fastq1`, `reads_fastq2`, `reads_fastq_unpaired`) from the task — BAM-only in array mode matches repo convention; FASTQ path was dead code given both workflow wrappers only accept BAMs
- Sample names derived inside bash loop via `$(basename $bam .bam)` — no `String sample_name` input needed in array mode
- Per-sample FASTQ files cleaned up (`rm -f`) inside loop to free disk during multi-sample runs
- centrifuger_multi exposes both `machine_mem_gb` (default 480) and `cpu` (default 16) at workflow level for runtime override

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All three files pass `miniwdl check`
- Both workflows are ready for test input JSONs and Dockstore registration (Phase 6)
- No blockers

---
*Phase: 05-centrifuger-single-and-centrifuger-multi-workflow-wrappers*
*Completed: 2026-04-02*
