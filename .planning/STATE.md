# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-12)

**Core value:** Reliable viral discovery from metagenomic assemblies must work. If genomad can identify viruses in the input, the pipeline must capture and annotate them correctly.
**Current focus:** Phase 3 - Terra Integration

## Current Position

Phase: 3 of 3 (Terra Integration)
Plan: 2 of 2 complete
Status: Complete - Phase 3 execution
Last activity: 2026-02-12 — Completed plan 03-02 (Gap Closure - True Optional Primitives)

Progress: [██████████] 100%

## Performance Metrics

**Velocity:**
- Total plans completed: 7
- Average duration: 2.6 min
- Total execution time: 0.31 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-core-tasks-and-single-sample-workflow | 3/3 | 7.6 min | 2.5 min |
| 02-multi-sample-workflow | 2/3 | 3.7 min | 1.8 min |
| 03-terra-integration | 2/2 | 8.1 min | 4.0 min |

**Recent Trend:**
- Last 5 plans: 02-01 (1.5 min), 02-02 (2.2 min), 03-01 (7.0 min), 03-02 (1.1 min)
- Trend: Stable (03-02 quick execution due to focused refactoring)

*Updated after each plan completion*

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

**Phase 1 (Complete):**
- Fixed resource allocation (32GB RAM, 8 CPU, 50GB disk) for genomad_end_to_end instead of dynamic sizing
- Use --splits 8 for memory safety with large genomad databases
- **REVISED:** Always-present File outputs instead of File? pattern (testing revealed glob()[0] complexity)
- **REMOVED:** min_length parameter (genomad v1.11.2 doesn't support --min-length option)
- **FIXED:** Database path needs /genomad_db suffix to match tarball structure
- Follow classify_single.wdl structure for consistency
- Use parameter_meta with patterns and category for Terra UI
- Expose both File outputs and primitive statistics for Terra integration
- Use G5012.3.fasta as test assembly (Filovirus test data)
- Created 706MB test database fixture with zstd -19 compression

**Phase 2 (Incomplete - 2/3 plans complete):**
- Single scatter block pattern for genomad_multi (only one task called, not multiple)
- No per-sample summary tasks in multi-sample workflows (matches classify_multi.wdl pattern)
- Maintain cleanup=true default from Phase 1 approved pattern
- Three test samples for multi-sample validation (G5012.3, LASV_NGA_2018_0026, LASV_NGA_2018_0097)
- Separate miniWDL and Cromwell JSON files (not symlinks) per project convention
- Docker authentication gate documented (deferred - not blocking Phase 3)

**Phase 3 (Complete - 2/2 plans complete):**
- WDL 1.0 limitation: Optional Int?/Float? outputs return 0 instead of undefined when no hits
- Removed provirus from report_genomad_summary per user decision
- Implemented correct virus sorting by score+length with tie-breaking
- Score rounding to 2 decimals and family-level taxonomy extraction with Title Case
- Array[File] glob pattern at task level for optional primitives (not Int?/Float? with else-0)
- Conditional extraction at workflow level using if-blocks for true null behavior
- Verified pattern produces null (not 0) in Terra data tables when zero hits found

### Pending Todos

None.

### Blockers/Concerns

None. Phase 1 tested and validated.

## Session Continuity

Last session: 2026-02-12 (Phase 3 execution - plan 03-02)
Stopped at: Completed 03-02-PLAN.md (Gap Closure - True Optional Primitives)
Resume file: None

## Phase 1 Completion Summary

**Implementation:** 3 plans executed (01-01, 01-02, 01-03)
**Testing:** Runtime validation with miniwdl + local Docker image
**Bugs Fixed:** 3 critical bugs discovered and fixed during testing
- Removed invalid --min-length parameter
- Fixed database path (/genomad_db suffix required)
- Simplified output handling (always-present files vs File? pattern)

**Test Results:**
- ✓ Classified G5012.3 as Filovirus (Filoviridae family)
- ✓ Memory efficient: 4GB used vs 32GB allocated
- ✓ Empty result handling verified (0 plasmids, 0 proviruses)
- ✓ Summary statistics correct (1 virus, taxonomy populated)

**Commits:** 8 total (implementation + bug fixes)
- Plans 01-01 through 01-03: Core implementation
- Commit 15173df6: Bug fixes from testing
