---
phase: 01-implement-wdl-tasks-and-workflows
plan: "04"
subsystem: wdl
tags: [wdl, miniwdl, metagenomics, virnucpro, dockstore]

# Dependency graph
requires:
  - phase: 01-implement-wdl-tasks-and-workflows
    provides: "Plans 01-03: WDL tasks, workflows, test inputs, dockstore entries (on worktree branch)"
provides:
  - "classify_virnucpro_contigs and classify_reads_by_contig tasks in tasks_metagenomics.wdl on ca-kb_python"
  - "Standalone workflows classify_virnucpro_contigs.wdl and classify_reads_by_contig.wdl on ca-kb_python"
  - "Test input JSONs for both new workflows on ca-kb_python"
  - "Dockstore entries for both new workflows on ca-kb_python"
  - "All three WDL files pass miniwdl check"
affects: [phase-01-verification, dockstore-registry]

# Tech tracking
tech-stack:
  added: []
  patterns: ["Extract clean lines from worktree branch using git show + sed to avoid full merge"]

key-files:
  created:
    - pipes/WDL/workflows/classify_virnucpro_contigs.wdl
    - pipes/WDL/workflows/classify_reads_by_contig.wdl
    - test/input/WDL/miniwdl-local/test_inputs-classify_virnucpro_contigs-local.json
    - test/input/WDL/miniwdl-local/test_inputs-classify_reads_by_contig-local.json
  modified:
    - pipes/WDL/tasks/tasks_metagenomics.wdl
    - .dockstore.yml

key-decisions:
  - "Used git show + sed line range to extract clean task content without merging unrelated upstream changes"
  - "New dockstore entries placed between classify_single and coverage_table, matching worktree placement"

patterns-established:
  - "Selective content extraction from worktree branch: git show branch:file | sed -n 'startline,endlinep' >> target"

requirements-completed: [REQ-15, REQ-16, REQ-17]

# Metrics
duration: 5min
completed: 2026-04-01
---

# Phase 01 Plan 04: Gap Closure Summary

**Fixed duplicated WDL content and applied clean classify_virnucpro_contigs and classify_reads_by_contig implementation to ca-kb_python, confirmed by miniwdl check exit 0 on all three WDL files**

## Performance

- **Duration:** ~5 min
- **Started:** 2026-04-01T16:00:00Z
- **Completed:** 2026-04-01T16:05:00Z
- **Tasks:** 2/2
- **Files modified:** 6

## Accomplishments
- Appended 477 lines of new task content (classify_virnucpro_contigs + classify_reads_by_contig) to tasks_metagenomics.wdl, deduplicated from worktree
- Created 36-line classify_virnucpro_contigs.wdl and 34-line classify_reads_by_contig.wdl standalone workflows
- Created test input JSONs for both workflows with correct workflow-level key format
- Added two dockstore entries for the new workflows
- Confirmed all three WDL files pass miniwdl check (exit 0, no syntax errors)

## Task Commits

Each task was committed atomically:

1. **Task 1: Extract clean WDL content and apply to ca-kb_python** - `ecd307bd` (feat)
2. **Task 2: Validate all WDL files with miniwdl check** - validation only, no file changes

**Plan metadata:** (docs commit below)

## Files Created/Modified
- `pipes/WDL/tasks/tasks_metagenomics.wdl` - Appended two new tasks (lines 1488-1964 from worktree)
- `pipes/WDL/workflows/classify_virnucpro_contigs.wdl` - Created standalone wrapper workflow (36 lines)
- `pipes/WDL/workflows/classify_reads_by_contig.wdl` - Created standalone wrapper workflow (34 lines)
- `test/input/WDL/miniwdl-local/test_inputs-classify_virnucpro_contigs-local.json` - 3-line JSON with single key
- `test/input/WDL/miniwdl-local/test_inputs-classify_reads_by_contig-local.json` - 4-line JSON with two keys
- `.dockstore.yml` - Added 6 lines (two new workflow entries) between classify_single and coverage_table

## Decisions Made
- Used `git show branch:file | sed -n 'X,Yp'` to extract only the new task lines without merging unrelated upstream changes (Docker version bumps, CI updates across 30+ files in the worktree branch)
- Placed new dockstore entries between classify_single and coverage_table, consistent with worktree branch placement

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None. The line counts matched exactly as documented in the plan context: worktree had 1964 clean lines (ca-kb_python had 1487), workflow files needed head -36 and head -34 respectively.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

Phase 01 is now fully implemented and verified on ca-kb_python:
- Both new WDL tasks exist in tasks_metagenomics.wdl (once each)
- Both standalone workflows exist and pass miniwdl check
- Both test input JSONs are present
- Dockstore has one entry per new workflow
- All verification gaps from 01-VERIFICATION.md are closed

Ready for phase transition.

---
*Phase: 01-implement-wdl-tasks-and-workflows*
*Completed: 2026-04-01*
