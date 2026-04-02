# Project: VirNucPro Contig & Read Classification WDL Pipelines

## What This Is

Three production-ready WDL task wrappers that extend the VirNucPro pipeline in `viral-pipelines` with post-processing steps: contig-level viral classification (from VirNucPro chunk scores), read-level classification (by joining PAF alignments to contig classifications), and read-level Kraken2 taxonomy annotation (by joining Kraken2 per-read output to a pre-built DuckDB taxonomy database). All three scripts originate from `respiratory-surveillance/bin/` and are promoted into this repo as first-class WDL tasks, validated with miniwdl, and registered in Dockstore.

## Core Value

Enables VirNucPro users to complete the full analysis chain — from raw reads through viral scoring to per-read classification — entirely within the existing WDL infrastructure on Terra/DNAnexus.

## Context

- **Codebase:** `viral-pipelines` — WDL workflows/tasks for viral NGS data analysis
- **Shipped:** v1.0 (2026-04-01) — Phase 1 complete, 4 plans, all tasks verified
- **Shipped:** v1.1 (2026-04-01) — Phase 2 complete, 2 plans, parse_kraken2_reads task verified
- **Shipped:** v2.0 (2026-04-01) — Phase 3 complete, 3 plans, summarize_kb_extract_reads task verified
- **Current state:** `tasks_metagenomics.wdl` contains all three tasks; three standalone workflows exist; all three registered in `.dockstore.yml`; all files pass `miniwdl check`
- **Tech stack:** WDL 1.0, Python 3 inline heredocs, DuckDB (runtime pip install), zstandard compression, `quay.io/broadinstitute/py3-bio:0.1.3`
- **Working branch:** `ca-kb_python` — implementation lives here; not yet merged to master
- **Next milestone:** Planning phase — new requirements to be defined

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
- ✓ `summarize_kb_extract_reads` WDL task in `tasks_metagenomics.wdl` — v2.0
- ✓ Standalone workflow `summarize_kb_extract_reads.wdl` — v2.0
- ✓ Test input JSON for `summarize_kb_extract_reads` — v2.0
- ✓ `.dockstore.yml` entry for `summarize_kb_extract_reads` — v2.0

## Current Milestone: v3.0 Centrifuger Taxonomic Classification WDL

**Goal:** Add Centrifuger as a first-class taxonomic classifier in viral-pipelines, with single-sample and multi-sample workflows modeled on the existing `classify_single`/`classify_multi` pattern — multi-sample mode designed to amortize the 200+ GB database load across all samples on a node.

**Target features:**
- `centrifuger` WDL task in `tasks_metagenomics.wdl` (wraps Centrifuger binary; accepts pre-built index, FASTQ inputs, outputs classification report and per-read results)
- `centrifuger_single.wdl` — standalone workflow for one sample per call (Terra scatter-friendly)
- `centrifuger_multi.wdl` — batch workflow that loads the database once and scatters classification across multiple samples on the same node (1–4 nodes for a full cohort)
- Test input JSONs for both workflows
- Dockstore registration entries for both workflows

### Validated

- ✓ `centrifuger` WDL task in `tasks_metagenomics.wdl` — v3.0 Phase 4 (input: `File centrifuger_db_tgz` + `String db_name`; outputs: classification_tsv, kreport, centrifuger_log)
- ✓ BAM→FASTQ conversion via picard SamToFastq with `/1`/`/2` suffixes — v3.0 Phase 4

### Active

(CFGR-03, CFGR-04 — centrifuger_single and centrifuger_multi workflows)

### Out of Scope

- Combined end-to-end workflow chaining Centrifuger with downstream steps — standalone only for now
- Custom Docker image build — will use existing tool image or runtime install
- Real test data files for `miniwdl run` execution — test JSONs use placeholders; deferred
- Parquet output for `classify_reads_by_contig` — TSV default matches repo conventions

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Tasks in `tasks_metagenomics.wdl` | Keep all virnucpro classification co-located | ✓ Good — clean grouping |
| py3-bio + pip install duckdb | Avoids new image build; py3-bio already has pandas | ✓ Good — no infra change needed |
| Scripts as inline `python3<<CODE` heredocs | Consistent with existing repo pattern | ✓ Good — passes miniwdl check |
| Standalone workflows only | User preference — no combined pipeline needed yet | ✓ Good — simpler, composable |
| WDL call alias required (`classify_contigs`, `classify_reads`, `parse_reads`) | WDL disallows call name = containing workflow name | ✓ Good — pattern documented, applied consistently |
| Selective `git show + sed` for content extraction | Avoided merging 30+ unrelated upstream changes from worktree branch | ✓ Good — clean cherry-pick of only new content |
| `allowNestedInputs: true` + workflow-level JSON keys | Enables direct workflow input in test JSONs without call-alias indirection | ✓ Good — matches miniwdl local run pattern |
| Dockstore entries without `testParameterFiles` | Placeholder JSON paths not for CI execution | — Pending — add when real test data exists |
| DuckDB-only heredoc extraction for `parse_kraken2_reads` | Exclude `TaxonomyDatabase`, `_write_parquet`, `argparse` — only DuckDB path needed | ✓ Good — leaner task, no unused dependencies |
| 8 GB / 1 CPU / `mem2_ssd1_v2_x2` for `parse_kraken2_reads` | DuckDB in-memory join is memory-bound, not CPU-bound | ✓ Good — locked in CONTEXT.md from research phase |
| py3-bio + runtime zstandard for `summarize_kb_extract_reads` | Consistent with v1.0/v1.1 pattern, no custom image needed | ✓ Good — no infrastructure changes |
| Inline Python heredoc for kb extract task | Matches existing repo conventions | ✓ Good — passes miniwdl check |
| Zstd compression for TSV output | Space-efficient, consistent with repo conventions | ✓ Good — smaller outputs |

## Evolution

This document evolves at phase transitions and milestone boundaries.

---
*Last updated: 2026-04-02 — Phase 4 complete: centrifuger WDL task added to tasks_metagenomics.wdl*
