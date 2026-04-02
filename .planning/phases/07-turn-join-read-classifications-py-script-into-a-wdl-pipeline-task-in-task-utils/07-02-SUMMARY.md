---
phase: 07-turn-join-read-classifications-py-script-into-a-wdl-pipeline-task-in-task-utils
plan: 02
subsystem: metagenomics
tags: [wdl, dockstore, workflow-wrapper, parquet, kallisto, kraken2, virnucpro, genomad]

# Dependency graph
requires:
  - phase: 07-01
    provides: join_read_classifications WDL task in tasks_metagenomics.wdl
provides:
  - join_read_classifications standalone workflow wrapper
  - Placeholder test input JSON for join_read_classifications
  - Dockstore registration entry for join_read_classifications
affects:
  - Dockstore discoverability of join_read_classifications on Terra/DNAnexus

# Tech tracking
tech-stack:
  added: []
  patterns:
    - Standalone workflow wrapper importing tasks_metagenomics.wdl with task alias (as join_reads)
    - allowNestedInputs: true for Terra-compatible placeholder test JSON

key-files:
  created:
    - pipes/WDL/workflows/join_read_classifications.wdl
    - test/input/WDL/miniwdl-local/test_inputs-join_read_classifications-local.json
  modified:
    - .dockstore.yml

key-decisions:
  - "Call alias as join_reads — WDL disallows call name = containing workflow name (D-10)"
  - "No testParameterFiles in Dockstore entry — placeholder paths not CI-runnable (D-09)"
  - "allowNestedInputs: true — required for Terra-compatible test JSON with workflow-level input keys"

patterns-established:
  - "Workflow wrapper with task alias avoids call-name = workflow-name collision"

requirements-completed:
  - JRC-02
  - JRC-03
  - JRC-04

# Metrics
duration: 2min
completed: 2026-04-02
---

# Phase 7 Plan 02: join_read_classifications Workflow Wrapper Summary

**Standalone workflow wrapper + placeholder test JSON + Dockstore entry for join_read_classifications — completing full task+workflow+test+registration treatment consistent with all prior phases**

## Performance

- **Duration:** 2 min
- **Started:** 2026-04-02T21:37:02Z
- **Completed:** 2026-04-02T21:38:35Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments

- Created `join_read_classifications.wdl` workflow wrapper importing `tasks_metagenomics.wdl` and calling the task with alias `join_reads` (per D-10)
- Added `allowNestedInputs: true` in meta block for Terra-compatible test JSON
- All 4 optional `File?` inputs and required `String sample_id` passed through to task; `classifications_parquet` exposed as workflow output
- `miniwdl check` exits 0 with no errors
- Created placeholder test input JSON with workflow-level keys for all 5 inputs
- Appended `join_read_classifications` entry to `.dockstore.yml` with `subclass: WDL`, no `testParameterFiles`

## Task Commits

Each task was committed atomically:

1. **Task 1: Create standalone workflow wrapper** - `d05801d0` (feat)
2. **Task 2: Create test input JSON and Dockstore entry** - `a812d49b` (feat)

## Files Created/Modified

- `pipes/WDL/workflows/join_read_classifications.wdl` - Standalone workflow wrapper (33 lines)
- `test/input/WDL/miniwdl-local/test_inputs-join_read_classifications-local.json` - Placeholder test JSON (7 lines)
- `.dockstore.yml` - Appended join_read_classifications entry (3 lines added)

## Decisions Made

- Used call alias `as join_reads` — WDL disallows a call with the same name as the enclosing workflow. This follows the established pattern (`parse_reads`, `classify_contigs`, `classify_reads`).
- No `testParameterFiles` in Dockstore entry — placeholder paths are not CI-runnable (consistent with all prior phase entries per STATE.md decision log).
- `allowNestedInputs: true` added to meta — required for Terra to accept the `workflow.input` key format in test JSON.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Phase 07 is complete — `join_read_classifications` task + workflow + test + registration fully implemented
- The workflow is discoverable on Dockstore and callable as a standalone workflow on Terra/DNAnexus

## Self-Check: PASSED

- FOUND: pipes/WDL/workflows/join_read_classifications.wdl
- FOUND: test/input/WDL/miniwdl-local/test_inputs-join_read_classifications-local.json
- FOUND: join_read_classifications entry in .dockstore.yml
- FOUND: commit d05801d0
- FOUND: commit a812d49b

---
*Phase: 07-turn-join-read-classifications-py-script-into-a-wdl-pipeline-task-in-task-utils*
*Completed: 2026-04-02*
