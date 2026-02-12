# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-12)

**Core value:** Reliable viral discovery from metagenomic assemblies must work. If genomad can identify viruses in the input, the pipeline must capture and annotate them correctly.
**Current focus:** Phase 2 - Multi-Sample Workflow

## Current Position

Phase: 2 of 3 (Multi-Sample Workflow)
Plan: None (ready to plan)
Status: Ready to plan Phase 2
Last activity: 2026-02-12 — Phase 1 complete, tested, and approved

Progress: [███░░░░░░░] 33%

## Performance Metrics

**Velocity:**
- Total plans completed: 3
- Average duration: 2.5 min
- Total execution time: 0.13 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-core-tasks-and-single-sample-workflow | 3/3 | 7.6 min | 2.5 min |

**Recent Trend:**
- Last 5 plans: 01-01 (3.3 min), 01-02 (2.1 min), 01-03 (2.2 min)
- Trend: Improving

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

### Pending Todos

None.

### Blockers/Concerns

None. Phase 1 tested and validated.

## Session Continuity

Last session: 2026-02-12 (Phase 1 execution, testing, and approval)
Stopped at: Phase 1 approved after successful testing - ready for Phase 2 planning
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
