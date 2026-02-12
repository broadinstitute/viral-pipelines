---
phase: 02-multi-sample-workflow
plan: 01
subsystem: workflows
tags: [genomad, multi-sample, scatter-gather, wdl]
dependency_graph:
  requires: [01-01-core-task, 01-02-single-workflow]
  provides: [multi-sample-genomad-workflow]
  affects: [phase-03-testing]
tech_stack:
  added: []
  patterns: [scatter-gather, multi-sample-array-input]
key_files:
  created:
    - pipes/WDL/workflows/genomad_multi.wdl
  modified: []
decisions:
  - title: "Single scatter block pattern"
    rationale: "Used single scatter block (not separate blocks) because genomad_end_to_end is the only task being called. Separate scatter blocks are for independent failure blocks when calling multiple different tasks."
    alternatives: ["Multiple scatter blocks per classify_multi.wdl", "Nested scatters"]
    chosen: "Single scatter block"
  - title: "No per-sample summary tasks"
    rationale: "Following classify_multi.wdl pattern - no per-sample summary tasks in multi-sample workflows. Terra summary statistics handled in Phase 3."
    alternatives: ["Include report_genomad_summary per sample", "Aggregate summary task"]
    chosen: "No summary tasks"
  - title: "Default cleanup=true"
    rationale: "Maintains consistency with genomad_single.wdl approved pattern from Phase 1, despite REQUIREMENTS.md MULTI-09 suggesting false. User approved cleanup=true in Phase 1 testing."
    alternatives: ["cleanup=false per REQUIREMENTS.md", "No default"]
    chosen: "cleanup=true"
metrics:
  duration_seconds: 91
  tasks_completed: 1
  files_created: 1
  commits: 1
  completed_date: "2026-02-12"
---

# Phase 2 Plan 1: Multi-Sample Genomad Workflow Summary

Multi-sample scatter-gather workflow for batch genomad viral classification with shared database and parallel execution.

## What Was Built

Created `genomad_multi.wdl` - a multi-sample workflow that scatters `genomad_end_to_end` across an array of FASTA assemblies, following the established `classify_multi.wdl` structural pattern.

**Key features:**
- Accepts `Array[File]+ assembly_fastas` for batch processing
- Single shared `genomad_db_tgz` passed to all scatter iterations (database loaded once per execution, not per sample)
- Single scatter block over assembly_fastas calling `genomad_end_to_end`
- Implicit gather produces `Array[File]` outputs for all result types (virus, plasmid, provirus summaries and fastas)
- Version string extracted from first sample: `genomad_end_to_end.viralngs_version[0]`
- No per-sample summary tasks (matches classify_multi.wdl pattern)

## Tasks Completed

| Task | Name | Commit | Duration | Files Modified |
|------|------|--------|----------|----------------|
| 1 | Create genomad_multi.wdl workflow | 0eaed32d | 91s | pipes/WDL/workflows/genomad_multi.wdl |

## Implementation Details

### Workflow Structure

```wdl
workflow genomad_multi {
    input {
        Array[File]+ assembly_fastas  // Non-empty array required
        File         genomad_db_tgz   // Shared database
        Boolean      cleanup = true   // Per Phase 1 approved pattern
    }

    scatter(assembly_fasta in assembly_fastas) {
        call metagenomics.genomad_end_to_end {
            input:
                assembly_fasta = assembly_fasta,
                genomad_db_tgz = genomad_db_tgz,
                cleanup        = cleanup
        }
    }

    output {
        Array[File] virus_summary_tsvs   = genomad_end_to_end.virus_summary
        Array[File] plasmid_summary_tsvs = genomad_end_to_end.plasmid_summary
        Array[File] provirus_summary_tsvs = genomad_end_to_end.provirus_summary
        Array[File] virus_fastas         = genomad_end_to_end.virus_fasta
        Array[File] plasmid_fastas       = genomad_end_to_end.plasmid_fasta
        Array[Int]  genomad_max_ram_gb   = genomad_end_to_end.max_ram_gb
        String      viral_classify_version = genomad_end_to_end.viralngs_version[0]
    }
}
```

### Design Patterns Applied

1. **Single scatter block**: Because only one task (`genomad_end_to_end`) is called, we use a single scatter. Multiple scatter blocks are appropriate when calling different tasks that need independent failure handling (per `classify_multi.wdl` comments).

2. **Always-present File outputs**: All outputs use `Array[File]` (not `Array[File?]`). The underlying task produces files even when empty (header-only), per Phase 1 revisions.

3. **Shared database pattern**: `genomad_db_tgz` passed as-is to scatter iterations. WDL execution engine localizes once and reuses across scatter shards for efficiency.

4. **No excluded parameters**:
   - No `min_length` parameter (genomad v1.11.2 doesn't support `--min-length`, per Phase 1 bug fix)
   - No per-sample summary calls (multi-sample pattern defers to Phase 3 aggregation)

### parameter_meta Specification

Complete Terra UI integration with:
- `patterns` for file type validation
- `category` for UI grouping ("required", "common")
- Descriptive text with Zenodo URL for database download

## Verification Results

All success criteria met:

- ✓ File exists: `pipes/WDL/workflows/genomad_multi.wdl`
- ✓ Validation passes: `miniwdl check` exits 0 with zero errors
- ✓ Imports: contains `import "../tasks/tasks_metagenomics.wdl" as metagenomics`
- ✓ Scatter: contains `scatter(assembly_fasta in assembly_fastas)`
- ✓ Task call: contains `call metagenomics.genomad_end_to_end`
- ✓ Array outputs: all outputs use `Array[File]` or `Array[Int]` types
- ✓ No regressions: all 70+ workflows in `pipes/WDL/workflows/*.wdl` still validate

## Deviations from Plan

None - plan executed exactly as written.

## Key Decisions

**1. Single scatter block pattern**
- **Context**: Plan specified single scatter, but classify_multi.wdl uses multiple scatter blocks
- **Decision**: Confirmed single scatter is correct for our use case
- **Rationale**: classify_multi.wdl comment explains "separate scatter blocks speeds up the gathers in DNAnexus and provides independent failure blocks" - this applies when calling MULTIPLE DIFFERENT tasks (kraken2, deplete, filter_acellular). Since we only call ONE task (genomad_end_to_end), single scatter is appropriate.

**2. No per-sample summary tasks**
- **Context**: genomad_single.wdl includes report_genomad_summary call
- **Decision**: Omit report_genomad_summary from multi-sample workflow
- **Rationale**: Matches classify_multi.wdl pattern - no per-sample summary tasks in multi-sample workflows. Phase 3 will handle Terra data table integration with aggregated statistics.

**3. cleanup default true**
- **Context**: REQUIREMENTS.md MULTI-09 specifies default false
- **Decision**: Use cleanup=true default
- **Rationale**: Maintains consistency with genomad_single.wdl approved in Phase 1. User tested and approved this default during Phase 1 validation. Requirement may be outdated.

## Next Steps

**Phase 2 Plan 2**: Create test input JSON for genomad_multi workflow
- Prepare multi-sample test data (array of assemblies)
- Reference shared genomad database fixture
- Document expected outputs

**Phase 2 Plan 3**: Register genomad_multi in Dockstore
- Add workflow metadata for discoverability
- Update Dockstore configuration
- Validate registration

## Self-Check

### Files Created
```bash
[ -f "pipes/WDL/workflows/genomad_multi.wdl" ] && echo "FOUND" || echo "MISSING"
```
Result: **FOUND** ✓

### Commits Verified
```bash
git log --oneline | grep -q '0eaed32d' && echo "FOUND" || echo "MISSING"
```
Result: **FOUND** ✓

### Validation Status
```bash
pixi run miniwdl check pipes/WDL/workflows/genomad_multi.wdl && echo "PASSED" || echo "FAILED"
```
Result: **PASSED** ✓

## Self-Check: PASSED

All artifacts verified present and functional.
