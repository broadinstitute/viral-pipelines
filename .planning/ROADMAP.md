# Roadmap: Genomad Pipeline Integration

## Overview

This roadmap delivers genomad viral discovery capability to viral-pipelines through three phases: Core Tasks and Single-Sample Workflow establishes the foundation with task module, single-sample processing, and cross-platform validation; Multi-Sample Workflow extends to batch processing with scatter-gather pattern; and Terra Integration adds workspace-friendly summary statistics. Each phase delivers working, testable functionality that builds on the previous phase.

## Phases

**Phase Numbering:**
- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

- [x] **Phase 1: Core Tasks and Single-Sample Workflow** - Task module, single-sample WDL, cross-platform validation
- [x] **Phase 2: Multi-Sample Workflow** - Batch processing with scatter-gather pattern
- [ ] **Phase 3: Terra Integration** - Workspace summary statistics and reporting

## Phase Details

### Phase 1: Core Tasks and Single-Sample Workflow
**Goal**: User can process individual metagenomic assemblies through genomad to identify viruses, plasmids, and proviruses with complete annotations and taxonomy
**Depends on**: Nothing (first phase)
**Requirements**: TASK-01 through TASK-19, SINGLE-01 through SINGLE-11, TEST-01, TEST-03, TEST-04, TEST-05, TEST-06, TEST-08, PLATFORM-01 through PLATFORM-08, DOC-01 through DOC-05
**Success Criteria** (what must be TRUE):
  1. User can submit single FASTA assembly to genomad_single.wdl and receive virus_summary.tsv, plasmid_summary.tsv, provirus.tsv, and annotated sequences
  2. Empty FASTA inputs produce valid empty outputs without task failure
  3. Workflow executes successfully on miniWDL local, Cromwell cloud, Terra workspace, and DNAnexus platform
  4. User-provided 5GB genomad database decompresses and executes without disk space errors
  5. Complex sample names with dots and dashes (e.g., "S20.l1.xxxx.2024-01-15") produce correctly named outputs
**Plans:** 3 plans

Plans:
- [x] 01-01-PLAN.md — Genomad task module (genomad_end_to_end + report_genomad_summary in tasks_metagenomics.wdl)
- [x] 01-02-PLAN.md — Single-sample workflow (genomad_single.wdl) + WDL validation
- [x] 01-03-PLAN.md — Test infrastructure (test input JSON, Dockstore registration, documentation)

### Phase 2: Multi-Sample Workflow
**Goal**: User can process batches of metagenomic assemblies in parallel, receiving per-sample outputs as arrays
**Depends on**: Phase 1
**Requirements**: MULTI-01 through MULTI-12, TEST-02, TEST-07, TEST-09, TEST-10
**Success Criteria** (what must be TRUE):
  1. User can submit array of 10+ FASTA assemblies to genomad_multi.wdl and receive arrays of outputs (one per sample)
  2. Scatter pattern parallelizes genomad execution across samples efficiently
  3. Partial failures (some samples empty, others succeed) produce valid outputs for successful samples
  4. Multi-sample test completes successfully with miniWDL and passes CI validation
**Plans:** 2 plans

Plans:
- [x] 02-01-PLAN.md — Multi-sample workflow (genomad_multi.wdl) with scatter-gather pattern
- [x] 02-02-PLAN.md — Test infrastructure (test input JSON, Dockstore registration, validation)

### Phase 3: Terra Integration
**Goal**: User can view genomad summary statistics directly in Terra workspace data tables without downloading TSV files
**Depends on**: Phase 2
**Requirements**: TERRA-01 through TERRA-16
**Success Criteria** (what must be TRUE):
  1. User sees total_viruses, total_plasmids as optional numeric columns in Terra data table (undefined when zero)
  2. User sees top_virus_name, top_virus_score, top_virus_taxonomy as columns in Terra data table
  3. Empty genomad results produce undefined counts and empty string text without task failure
  4. Summary task follows report_primary_kraken_taxa pattern for consistency with existing workflows
**Plans:** 2 plans

Plans:
- [x] 03-01-PLAN.md — Rewrite report_genomad_summary task + update workflow outputs (sorting, optional types, family taxonomy, provirus removal)
- [ ] 03-02-PLAN.md — Gap closure: refactor optional outputs for true null (not 0) using Array[File] + workflow-level conditional extraction

## Progress

**Execution Order:**
Phases execute in numeric order: 1 → 2 → 3

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Core Tasks and Single-Sample Workflow | 3/3 | ✓ Complete | 2026-02-12 |
| 2. Multi-Sample Workflow | 2/2 | ✓ Complete | 2026-02-12 |
| 3. Terra Integration | 1/2 | In progress | - |

---
*Roadmap created: 2026-02-12*
*Last updated: 2026-02-12 — Phase 3 gap closure plan added (03-02)*
