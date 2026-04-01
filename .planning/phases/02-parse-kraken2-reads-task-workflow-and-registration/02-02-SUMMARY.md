---
phase: 02-parse-kraken2-reads-task-workflow-and-registration
plan: "02"
subsystem: wdl
tags: [wdl, miniwdl, dockstore, kraken2, metagenomics, duckdb]

# Dependency graph
requires:
  - phase: 02-parse-kraken2-reads-task-workflow-and-registration
    provides: "parse_kraken2_reads WDL task in tasks_metagenomics.wdl (plan 02-01)"
provides:
  - "Standalone workflow pipes/WDL/workflows/parse_kraken2_reads.wdl wrapping the task"
  - "Placeholder test input JSON test/input/WDL/miniwdl-local/test_inputs-parse_kraken2_reads-local.json"
  - "Dockstore registration entry for parse_kraken2_reads in .dockstore.yml"
affects: [terra, dockstore, ci-validation, miniwdl-tests]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Standalone workflow thin wrapper: import task, mirror inputs (no docker), single call with alias, mirror output"
    - "Call alias required when task name matches workflow name to avoid WDL name conflict"
    - "Dockstore entry: name + subclass + primaryDescriptorPath, no testParameterFiles for placeholder paths"

key-files:
  created:
    - pipes/WDL/workflows/parse_kraken2_reads.wdl
    - test/input/WDL/miniwdl-local/test_inputs-parse_kraken2_reads-local.json
  modified:
    - .dockstore.yml

key-decisions:
  - "Use 'as parse_reads' call alias to avoid WDL name conflict (call name cannot equal containing workflow name)"
  - "Insert Dockstore entry between classify_reads_by_contig and coverage_table (existing file position)"
  - "No testParameterFiles in Dockstore entry (placeholder paths not for CI execution)"

patterns-established:
  - "Call alias pattern: when task and workflow share a name, use `call task as alias`"

requirements-completed: [REQ-08, REQ-09, REQ-10]

# Metrics
duration: 2min
completed: 2026-04-01
---

# Phase 02 Plan 02: parse_kraken2_reads Workflow and Registration Summary

**Standalone WDL workflow parse_kraken2_reads.wdl wrapping the DuckDB taxonomy annotation task, with placeholder test JSON and Dockstore registration entry**

## Performance

- **Duration:** 2 min
- **Started:** 2026-04-01T16:40:10Z
- **Completed:** 2026-04-01T16:42:00Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- Created `pipes/WDL/workflows/parse_kraken2_reads.wdl` — thin wrapper workflow with `allowNestedInputs: true`, mirrors task inputs (no docker), passes miniwdl check
- Created `test/input/WDL/miniwdl-local/test_inputs-parse_kraken2_reads-local.json` — placeholder test input JSON with workflow-scoped keys
- Added Dockstore entry in `.dockstore.yml` with correct `primaryDescriptorPath` and no `testParameterFiles`

## Task Commits

Each task was committed atomically:

1. **Task 1: Create standalone workflow and test input JSON** - `f7efd937` (feat)
2. **Task 2: Add Dockstore registration entry** - `8053d767` (feat)

**Plan metadata:** (docs commit follows)

## Files Created/Modified
- `pipes/WDL/workflows/parse_kraken2_reads.wdl` - Standalone workflow wrapping parse_kraken2_reads task with meta block, 4 inputs, single call alias, 1 output
- `test/input/WDL/miniwdl-local/test_inputs-parse_kraken2_reads-local.json` - Placeholder test JSON with parse_kraken2_reads.kraken2_reads_output and parse_kraken2_reads.taxonomy_db keys
- `.dockstore.yml` - Added parse_kraken2_reads entry after classify_reads_by_contig, before coverage_table

## Decisions Made
- Used `call metagenomics.parse_kraken2_reads as parse_reads` alias because WDL disallows a call name equal to the containing workflow name (miniwdl error: "Call's name may not equal the containing workflow's"). Output reference updated to `parse_reads.read_taxonomy`.
- Dockstore entry inserted between `classify_reads_by_contig` and `coverage_table` entries per plan's specified position (not strictly alphabetical, consistent with file's existing ordering).

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Added call alias to resolve WDL name conflict**
- **Found during:** Task 1 (Create standalone workflow)
- **Issue:** miniwdl check failed: "Call's name may not equal the containing workflow's" — `call metagenomics.parse_kraken2_reads` creates a call named `parse_kraken2_reads` which conflicts with the workflow name `parse_kraken2_reads`
- **Fix:** Changed call to `call metagenomics.parse_kraken2_reads as parse_reads` and updated output reference from `parse_kraken2_reads.read_taxonomy` to `parse_reads.read_taxonomy`
- **Files modified:** pipes/WDL/workflows/parse_kraken2_reads.wdl
- **Verification:** `miniwdl check` exits 0
- **Committed in:** f7efd937 (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Required fix for WDL validity. The call alias pattern is consistent with `classify_reads_by_contig.wdl` which also uses `as classify_reads`. No scope creep.

## Issues Encountered
None beyond the call alias fix documented above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All three deliverables complete: workflow (REQ-08), test JSON (REQ-09), Dockstore entry (REQ-10)
- Phase 02 fully complete — parse_kraken2_reads task, workflow, and registration all done
- Workflow is discoverable on Dockstore and testable with miniwdl (using real inputs)

---
*Phase: 02-parse-kraken2-reads-task-workflow-and-registration*
*Completed: 2026-04-01*

## Self-Check: PASSED

- FOUND: pipes/WDL/workflows/parse_kraken2_reads.wdl (in worktree)
- FOUND: test/input/WDL/miniwdl-local/test_inputs-parse_kraken2_reads-local.json (in worktree)
- FOUND: .dockstore.yml parse_kraken2_reads entry
- FOUND: commit f7efd937 (Task 1)
- FOUND: commit 8053d767 (Task 2)
