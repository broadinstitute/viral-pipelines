# Project: VirNucPro Contig & Read Classification WDL Pipelines

## What This Is

Two production-ready WDL task wrappers that extend the VirNucPro pipeline in `viral-pipelines` with post-processing steps: contig-level viral classification (from VirNucPro chunk scores) and read-level classification (by joining PAF alignments to contig classifications). Both scripts originate from `respiratory-surveillance/bin/` and are promoted into this repo as first-class WDL tasks, validated with miniwdl, and registered in Dockstore.

## Core Value

Enables VirNucPro users to complete the full analysis chain — from raw reads through viral scoring to per-read classification — entirely within the existing WDL infrastructure on Terra/DNAnexus.

## Context

- **Codebase:** `viral-pipelines` — WDL workflows/tasks for viral NGS data analysis
- **Shipped:** v1.0 (2026-04-01) — Phase 1 complete, 4 plans, all tasks verified
- **Shipped:** v1.1 (2026-04-01) — Phase 2 complete, 2 plans, parse_kraken2_reads task verified
- **Current state:** `tasks_metagenomics.wdl` contains all three tasks; three standalone workflows exist; all three registered in `.dockstore.yml`; all files pass `miniwdl check`
- **Tech stack:** WDL 1.0, Python 3 inline heredocs, DuckDB (runtime pip install), `quay.io/broadinstitute/py3-bio:0.1.3`
- **Working branch:** `ca-kb_python` — implementation lives here; not yet merged to master

## Requirements

### Validated

- ✓ WDL task/workflow pattern established (`tasks_metagenomics.wdl`, `classify_virnucpro_single.wdl`) — existing
- ✓ Docker image pattern: `quay.io/broadinstitute/py3-bio:0.1.3` for Python tasks — existing
- ✓ Inline `python3<<CODE` heredoc pattern for embedding scripts — v1.0
- ✓ Test input JSONs in `test/input/WDL/miniwdl-local/` — v1.0
- ✓ `classify_virnucpro_contigs` WDL task in `tasks_metagenomics.wdl` — v1.0
- ✓ `classify_reads_by_contig` WDL task in `tasks_metagenomics.wdl` — v1.0
- ✓ Standalone workflow `classify_virnucpro_contigs.wdl` — v1.0
- ✓ Standalone workflow `classify_reads_by_contig.wdl` — v1.0
- ✓ Both tasks pass `miniwdl check` validation — v1.0
- ✓ Test input JSONs for both workflows — v1.0 (placeholder paths; real test data deferred)
- ✓ Docker image strategy: py3-bio + runtime pip install for duckdb — v1.0
- ✓ `.dockstore.yml` entries for both workflows — v1.0
- ✓ `parse_kraken2_reads` WDL task in `tasks_metagenomics.wdl` — v1.1
- ✓ Standalone workflow `parse_kraken2_reads.wdl` — v1.1
- ✓ Test input JSON for `parse_kraken2_reads` — v1.1
- ✓ `.dockstore.yml` entry for `parse_kraken2_reads` — v1.1

## Completed Milestone: v1.1 Kraken2 Read Taxonomy Annotation WDL Task

**Shipped:** 2026-04-01 — Phase 2 complete (2/2 plans), verification passed (5/5 must-haves)

**Delivered:**
- ✓ `parse_kraken2_reads` task in `tasks_metagenomics.wdl` (DuckDB path only; 6-column TSV output; 8 GB / 1 CPU; py3-bio + pip install duckdb runtime)
- ✓ Standalone workflow `parse_kraken2_reads.wdl` (uses `as parse_reads` call alias per WDL constraint)
- ✓ Test input JSON in `test/input/WDL/miniwdl-local/`
- ✓ `.dockstore.yml` entry for the new workflow

### Active

### Out of Scope

- Combined end-to-end workflow chaining all three stages — user chose standalone only
- Custom Docker image build — pip install at runtime is sufficient for now
- Parquet output for `classify_reads_by_contig` — TSV default matches repo conventions
- Real test data files for `miniwdl run` execution — test JSONs use placeholders; deferred

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Tasks in `tasks_metagenomics.wdl` | Keep all virnucpro classification co-located | ✓ Good — clean grouping |
| py3-bio + pip install duckdb | Avoids new image build; py3-bio already has pandas | ✓ Good — no infra change needed |
| Scripts as inline `python3<<CODE` heredocs | Consistent with existing repo pattern | ✓ Good — passes miniwdl check |
| Standalone workflows only | User preference — no combined pipeline needed yet | ✓ Good — simpler, composable |
| WDL call alias required (`classify_contigs`, `classify_reads`) | WDL disallows call name = containing workflow name | ✓ Good — discovered during execution, now documented |
| Selective `git show + sed` for content extraction | Avoided merging 30+ unrelated upstream changes from worktree branch | ✓ Good — clean cherry-pick of only new content |
| `allowNestedInputs: true` + workflow-level JSON keys | Enables direct workflow input in test JSONs without call-alias indirection | ✓ Good — matches miniwdl local run pattern |
| Dockstore entries without `testParameterFiles` | Placeholder JSON paths not for CI execution | — Pending — add when real test data exists |

## Evolution

This document evolves at phase transitions and milestone boundaries.

---
*Last updated: 2026-04-01 — Milestone v1.1 complete (Phase 2)*
