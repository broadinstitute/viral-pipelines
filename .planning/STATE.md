# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-12)

**Core value:** Reliable viral discovery from metagenomic assemblies must work. If genomad can identify viruses in the input, the pipeline must capture and annotate them correctly.
**Current focus:** Phase 1 - Core Tasks and Single-Sample Workflow

## Current Position

Phase: 1 of 3 (Core Tasks and Single-Sample Workflow)
Plan: 3 of 3 (completed)
Status: Complete
Last activity: 2026-02-12 — Completed plan 01-03 (testing infrastructure and platform registration)

Progress: [██████████] 100%

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

- Fixed resource allocation (32GB RAM, 8 CPU, 50GB disk) for genomad_end_to_end instead of dynamic sizing (Plan 01-01)
- Use --splits 8 for memory safety with large genomad databases (Plan 01-01)
- Optional File? outputs with glob()[0] pattern to handle empty results gracefully (Plan 01-01)
- Use default='' pattern for optional File? inputs in command blocks to handle WDL substitution (Plan 01-01)
- [Phase 01-02]: Follow classify_single.wdl structure for consistency with existing workflows
- [Phase 01-02]: Use parameter_meta with patterns and category fields for Terra UI integration
- [Phase 01-02]: Expose both optional File? outputs and primitive statistics for flexible downstream usage
- [Phase 01-03]: Use G5012.3.fasta as test assembly input (existing viral test data)
- [Phase 01-03]: Reference genomad_db-tinytest.tar.zst for test database (to be created separately)
- [Phase 01-03]: Follow existing pattern: JSON in parent WDL/ directory with symlink in miniwdl-local/

### Pending Todos

None yet.

### Blockers/Concerns

None yet.

## Session Continuity

Last session: 2026-02-12 (plan 01-03 execution)
Stopped at: Completed 01-03-PLAN.md (testing infrastructure and platform registration) - Phase 1 complete
Resume file: None
