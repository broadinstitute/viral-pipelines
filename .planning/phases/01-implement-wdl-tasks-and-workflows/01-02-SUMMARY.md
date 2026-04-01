---
phase: 01-implement-wdl-tasks-and-workflows
plan: 02
subsystem: wdl-tasks
tags: [wdl, duckdb, python, metagenomics, virnucpro, paf, classification]

# Dependency graph
requires:
  - phase: 01-implement-wdl-tasks-and-workflows
    plan: 01
    provides: classify_virnucpro_contigs task that produces contig_classifications TSV consumed here
provides:
  - classify_reads_by_contig task in tasks_metagenomics.wdl with full inline DuckDB Python pipeline
  - classify_reads_by_contig.wdl standalone workflow wrapping the task
affects: [any workflow chaining PAF alignment to read-level viral classification]

# Tech tracking
tech-stack:
  added: [duckdb (runtime pip install)]
  patterns: [python3<<CODE heredoc with pip install before script, WDL call alias to avoid name collision]

key-files:
  created:
    - pipes/WDL/workflows/classify_reads_by_contig.wdl
  modified:
    - pipes/WDL/tasks/tasks_metagenomics.wdl

key-decisions:
  - "WDL call alias 'classify_reads' required — call name cannot equal containing workflow name"
  - "Parquet output path removed entirely — WDL version always outputs TSV via DuckDB COPY"
  - "pip install duckdb --quiet --no-cache-dir placed before python3<<CODE heredoc per REQ-11"

patterns-established:
  - "Call alias pattern: when task name equals workflow name, use 'as <alias>' on the call"
  - "DuckDB tasks: pip install before heredoc, TSV output via COPY with FORMAT CSV DELIMITER tab"

requirements-completed: [REQ-07, REQ-08, REQ-09, REQ-10, REQ-11, REQ-12, REQ-14, REQ-18]

# Metrics
duration: 4min
completed: 2026-04-01
---

# Phase 01 Plan 02: classify_reads_by_contig Task and Workflow Summary

**classify_reads_by_contig WDL task with full inline DuckDB SQL pipeline joining PAF alignments to contig classifications, plus standalone wrapper workflow**

## Performance

- **Duration:** 4 min
- **Started:** 2026-04-01T14:10:12Z
- **Completed:** 2026-04-01T14:14:14Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- Added `classify_reads_by_contig` task to `tasks_metagenomics.wdl` after `classify_virnucpro_contigs`
- Embedded full `prepare_paf_file()` function and DuckDB SQL aggregation pipeline inline
- Created standalone `classify_reads_by_contig.wdl` workflow wrapping the task
- Both files pass `miniwdl check` with no new errors

## Task Commits

Each task was committed atomically:

1. **Task 1: Add classify_reads_by_contig task to tasks_metagenomics.wdl** - `17cb2e34` (feat)
2. **Task 2: Create classify_reads_by_contig standalone workflow** - `111a2ce7` (feat)

**Plan metadata:** (docs commit follows)

## Files Created/Modified
- `pipes/WDL/tasks/tasks_metagenomics.wdl` - Added classify_reads_by_contig task (303 lines) after classify_virnucpro_contigs
- `pipes/WDL/workflows/classify_reads_by_contig.wdl` - New standalone workflow wrapping the task

## Decisions Made
- Used call alias `classify_reads` (not `classify_reads_by_contig`) because WDL disallows a call name equal to the containing workflow name — same pattern already established in Plan 01 for `classify_virnucpro_contigs.wdl`
- Removed parquet output entirely; WDL version always outputs TSV via DuckDB COPY statement
- Kept `prepare_paf_file()` function verbatim — handles both `.gz` and plain PAF transparently

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Used call alias to fix WDL name collision**
- **Found during:** Task 2 (workflow creation + miniwdl check)
- **Issue:** `miniwdl check` returned "Call's name may not equal the containing workflow's" — plan spec showed bare `call metagenomics.classify_reads_by_contig` without alias
- **Fix:** Added `as classify_reads` alias on the call; updated output reference to `classify_reads.read_classifications`
- **Files modified:** pipes/WDL/workflows/classify_reads_by_contig.wdl
- **Verification:** `miniwdl check pipes/WDL/workflows/classify_reads_by_contig.wdl` exits 0
- **Committed in:** 111a2ce7 (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (Rule 1 - Bug, WDL name collision)
**Impact on plan:** Required fix — without alias the workflow is syntactically invalid. Same pattern already documented in STATE.md from Plan 01.

## Issues Encountered
- Worktree `worktree-agent-ade8ca6e` was missing Plan 01 commits (classify_virnucpro_contigs task). Cherry-picked commits `e8ba90a3` and `1c4d2314` from `worktree-agent-a2ead822` before proceeding. No content conflicts.

## Known Stubs
None - task embeds complete DuckDB pipeline logic; no placeholder logic or hardcoded empty returns.

## Next Phase Readiness
- Phase 01 now complete: both classify_virnucpro_contigs and classify_reads_by_contig tasks and workflows implemented
- Ready for test input JSON creation (AGENTS.md requirement) if a Phase 02 exists
- Full VirNucPro analysis chain (raw reads -> viral scoring -> contig classification -> read classification) is now available in WDL

## Self-Check: PASSED

---
*Phase: 01-implement-wdl-tasks-and-workflows*
*Completed: 2026-04-01*
