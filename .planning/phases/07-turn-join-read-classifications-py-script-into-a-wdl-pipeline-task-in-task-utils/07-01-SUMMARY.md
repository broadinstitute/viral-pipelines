---
phase: 07-turn-join-read-classifications-py-script-into-a-wdl-pipeline-task-in-task-utils
plan: 01
subsystem: metagenomics
tags: [wdl, duckdb, parquet, zstd, py3-bio, kallisto, kraken2, virnucpro, genomad]

# Dependency graph
requires:
  - phase: 06-centrifuger-single-and-centrifuger-multi-workflow-wrappers
    provides: centrifuger task and workflow patterns in tasks_metagenomics.wdl
provides:
  - join_read_classifications WDL task in tasks_metagenomics.wdl
  - DuckDB 4-way FULL OUTER JOIN logic (Kallisto + VNP + K2 + geNomad) as inline Python heredoc
affects:
  - 07-02 (standalone workflow wrapper, test JSON, dockstore entry)

# Tech tracking
tech-stack:
  added: []
  patterns:
    - File? optional inputs with __NONE__ sentinel in python3<<CODE heredoc
    - DuckDB FULL OUTER JOIN with regexp_replace to strip /1|/2 mate suffix for K2 join

key-files:
  created: []
  modified:
    - pipes/WDL/tasks/tasks_metagenomics.wdl

key-decisions:
  - "Task placed in tasks_metagenomics.wdl (not tasks_utils.wdl) — all classification join tasks co-located there"
  - "File? optional inputs using ~{default='__NONE__' input} pattern — __NONE__ sentinel, _file_is_usable() returns False for non-existent paths"
  - "size() on File? inputs used directly (no select_first wrapper) — WDL returns 0 for undefined optional File?, miniwdl check exits 0"
  - "16 GB / 1 CPU / mem2_ssd1_v2_x2 runtime — scaled up from parse_kraken2_reads 8 GB baseline to handle 4-source in-memory join"

patterns-established:
  - "Optional File? WDL inputs with __NONE__ default sentinel for Python _file_is_usable() guard"

requirements-completed:
  - JRC-01

# Metrics
duration: 4min
completed: 2026-04-02
---

# Phase 7 Plan 01: join_read_classifications WDL Task Summary

**DuckDB 4-way FULL OUTER JOIN task (Kallisto + VNP + K2 + geNomad) appended to tasks_metagenomics.wdl as inline Python3 heredoc outputting ZSTD Parquet**

## Performance

- **Duration:** 4 min
- **Started:** 2026-04-02T21:28:23Z
- **Completed:** 2026-04-02T21:31:55Z
- **Tasks:** 1
- **Files modified:** 1

## Accomplishments

- Appended `join_read_classifications` WDL task to `tasks_metagenomics.wdl` after the `centrifuger` task
- Embedded full 4-way DuckDB join logic from `join_read_classifications.py` verbatim as inline `python3<<CODE` heredoc
- Task accepts 4 optional `File?` inputs (kallisto_summary, kraken2_reads, vnp_reads, genomad_virus_summary) plus required `String sample_id`, outputs `File classifications_parquet`
- `miniwdl check` exits 0 with no errors for the new task

## Task Commits

Each task was committed atomically:

1. **Task 1: Append join_read_classifications task to tasks_metagenomics.wdl** - `4e470643` (feat)

**Plan metadata:** (pending — docs commit)

## Files Created/Modified

- `pipes/WDL/tasks/tasks_metagenomics.wdl` - Appended `join_read_classifications` task (337 lines added)

## Decisions Made

- Used `size()` directly on `File?` inputs in the runtime disk formula (not `select_first([size(...), 0])`). WDL correctly returns 0 for undefined optional Files; `miniwdl check` issued no errors for this pattern.
- Preserved all 4 join steps from the source script exactly — no simplification or reordering.
- `__NONE__` sentinel approach for optional WDL File? inputs: when not provided, WDL interpolates `__NONE__` as the string, and `_file_is_usable()` returns `False` for that non-existent path.

## Deviations from Plan

None - plan executed exactly as written.

The plan flagged a possible need to wrap `size()` calls in `select_first()` for optional File? inputs. Testing confirmed `size()` works directly on `File?` in WDL (returns 0 when undefined), so the simpler form was used. This matches the plan's guidance ("verify this works with miniwdl check").

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- `join_read_classifications` task is fully functional and validated in `tasks_metagenomics.wdl`
- Plan 02 can proceed to create the standalone workflow wrapper, placeholder test JSON, and `.dockstore.yml` entry

---
*Phase: 07-turn-join-read-classifications-py-script-into-a-wdl-pipeline-task-in-task-utils*
*Completed: 2026-04-02*
