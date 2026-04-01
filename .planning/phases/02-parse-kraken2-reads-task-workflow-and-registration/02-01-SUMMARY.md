---
phase: 02-parse-kraken2-reads-task-workflow-and-registration
plan: "01"
subsystem: wdl
tags: [wdl, kraken2, taxonomy, duckdb, metagenomics, py3-bio]

# Dependency graph
requires:
  - phase: 01-implement-wdl-tasks-and-workflows
    provides: classify_reads_by_contig task structure (neighbor pattern for new task)
provides:
  - parse_kraken2_reads WDL task in tasks_metagenomics.wdl with inline DuckDB taxonomy annotation Python
affects:
  - 02-02 (workflow wrapping parse_kraken2_reads)
  - Dockstore registration

# Tech tracking
tech-stack:
  added: [duckdb (pip-installed at runtime via py3-bio:0.1.3 image)]
  patterns:
    - python3<<CODE heredoc with pip install preamble (consistent with classify_reads_by_contig)
    - Dynamic disk sizing with ceil(size()) expressions embedded directly in runtime strings (no variable)
    - Boolean WDL-to-Python interpolation via ~{true="True" false="False" var}

key-files:
  created: []
  modified:
    - pipes/WDL/tasks/tasks_metagenomics.wdl

key-decisions:
  - "parse_kraken2_reads task appended after classify_reads_by_contig in tasks_metagenomics.wdl"
  - "DuckDBTaxonomyDatabase class extracted verbatim from parse_kraken2_reads.py (lines 169-238)"
  - "output_format parameter removed from parse_kraken2_output; _write_tsv called unconditionally"
  - "resolve_strains handled at DB load time in DuckDBTaxonomyDatabase.__init__; driver block does not pass it to parse_kraken2_output"
  - "Dynamic disk: ceil(size(kraken2_reads_output)*3 + size(taxonomy_db) + 20) embedded directly (no disk_size variable per Pitfall 3)"

patterns-established:
  - "parse_kraken2_reads: kraken2 reads + duckdb taxonomy -> 6-col TSV (SAMPLE_ID, READ_ID, TAXONOMY_ID, TAX_NAME, KINGDOM, TAX_RANK)"

requirements-completed: [REQ-01, REQ-02, REQ-03, REQ-04, REQ-05, REQ-06, REQ-07, REQ-11]

# Metrics
duration: 15min
completed: 2026-04-01
---

# Phase 02 Plan 01: parse_kraken2_reads WDL Task Summary

**parse_kraken2_reads WDL task with inline DuckDB taxonomy annotation Python, outputting 6-column TSV of per-read NCBI taxonomy (SAMPLE_ID, READ_ID, TAXONOMY_ID, TAX_NAME, KINGDOM, TAX_RANK)**

## Performance

- **Duration:** ~15 min
- **Started:** 2026-04-01T16:20:00Z
- **Completed:** 2026-04-01T16:35:20Z
- **Tasks:** 1
- **Files modified:** 1

## Accomplishments
- Appended `parse_kraken2_reads` task to `tasks_metagenomics.wdl` immediately after `classify_reads_by_contig`
- Python heredoc extracts `DuckDBTaxonomyDatabase`, `parse_kraken2_output`, and `_write_tsv` from source `parse_kraken2_reads.py` with adaptation (output_format removed, TSV-only)
- Task passes `miniwdl check` with no new errors
- Full `parameter_meta` block with patterns for File inputs and description+category for non-File inputs
- Runtime sized at 8 GB / 1 CPU / mem2_ssd1_v2_x2 per requirements

## Task Commits

1. **Task 1: Implement parse_kraken2_reads WDL task with Python heredoc** - `7303803b` (feat)

**Plan metadata:** (see final docs commit)

## Files Created/Modified
- `pipes/WDL/tasks/tasks_metagenomics.wdl` - Added `task parse_kraken2_reads` (232 lines) after `classify_reads_by_contig`

## Decisions Made
- Extracted only `DuckDBTaxonomyDatabase`, `parse_kraken2_output`, `_write_tsv` from source — excluded `TaxonomyDatabase`, `_write_parquet`, argparse, pathlib, pyarrow
- `resolve_strains` handled at `DuckDBTaxonomyDatabase.__init__` time; driver block does not pass it to `parse_kraken2_output` since strain resolution already applied during DB load
- Disk expression embedded directly in runtime strings (no `disk_size` variable) to match plan Pitfall 3 guidance and avoid unused declaration lint warnings
- Merged `ca-kb_python` branch into worktree to obtain the full 1953-line `tasks_metagenomics.wdl` (worktree started from an older base without Phase 01 additions)

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Merged ca-kb_python into worktree before task implementation**
- **Found during:** Setup (before Task 1)
- **Issue:** The worktree branch `worktree-agent-ad0e29cb` was based on an older master commit; `tasks_metagenomics.wdl` had only 943 lines and lacked `classify_reads_by_contig` (added in Phase 01). The plan requires appending after that task.
- **Fix:** Ran `git merge ca-kb_python` to bring in all Phase 01 WDL additions. Resolved a minor `requirements-modules.txt` conflict (kept newer pangolin/nextclade versions from HEAD, added virnucpro-cuda from ca-kb_python).
- **Files modified:** All Phase 01 WDL files + `requirements-modules.txt`
- **Verification:** `tasks_metagenomics.wdl` now 1953 lines; `classify_reads_by_contig` present; merge commit clean
- **Committed in:** `630bdb55` (merge + conflict resolution commit)

---

**Total deviations:** 1 auto-fixed (Rule 3 - blocking)
**Impact on plan:** Merge required to access Phase 01 additions. No scope creep.

## Issues Encountered
None beyond the worktree base mismatch resolved by merging ca-kb_python.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- `parse_kraken2_reads` task is ready for wrapping in a workflow (plan 02-02)
- No blockers

---
*Phase: 02-parse-kraken2-reads-task-workflow-and-registration*
*Completed: 2026-04-01*
