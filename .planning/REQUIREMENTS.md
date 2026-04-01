# Requirements: v1.1 Kraken2 Read Taxonomy Annotation WDL Task

## Overview

Add a WDL task and standalone workflow to wrap `parse_kraken2_reads.py` from `respiratory-surveillance/bin/`. The task reads Kraken2 per-read output and annotates each read with NCBI taxonomy (name, kingdom, rank) using a pre-built DuckDB taxonomy database. Follows the same inline-heredoc and runtime patterns established in v1.0.

## Functional Requirements

### Task: parse_kraken2_reads

- [ ] **REQ-01** — Task `parse_kraken2_reads` added to `pipes/WDL/tasks/tasks_metagenomics.wdl` after the existing `classify_reads_by_contig` task
- [ ] **REQ-02** — Inputs: `File kraken2_reads_output` (Kraken2 per-read output, may be gzipped), `File taxonomy_db` (pre-built DuckDB file from build_taxonomy_db.py), `String sample_id` (default derived from filename stem), `Boolean resolve_strains` (default false)
- [ ] **REQ-03** — Output: `File read_taxonomy` — TSV with columns: SAMPLE_ID, READ_ID, TAXONOMY_ID, TAX_NAME, KINGDOM, TAX_RANK
- [ ] **REQ-04** — Python logic embedded inline as `python3<<CODE` heredoc; DuckDB taxonomy path only — no nodes.dmp/names.dmp fallback
- [ ] **REQ-05** — Docker: `quay.io/broadinstitute/py3-bio:0.1.3`; command preamble runs `pip install duckdb --quiet --no-cache-dir` before the Python block; `set -e` guards the install
- [ ] **REQ-06** — Runtime: 8 GB memory, 1 CPU, dynamic disk size (`ceil(size(kraken2_reads_output)*3 + size(taxonomy_db) + 20)`); `disks`/`disk` dual annotation, `dx_instance_type: "mem2_ssd1_v2_x2"`, `preemptible: 2`, `maxRetries: 2`
- [ ] **REQ-07** — Task passes `miniwdl check` with no errors

### Standalone Workflow

- [ ] **REQ-08** — `pipes/WDL/workflows/parse_kraken2_reads.wdl` — thin workflow wrapping the `parse_kraken2_reads` task; `meta` block includes `description`, `author`, `email`, and `allowNestedInputs: true`

### Validation & Registration

- [ ] **REQ-09** — Test input JSON `test/input/WDL/miniwdl-local/test_inputs-parse_kraken2_reads-local.json` with placeholder file paths
- [ ] **REQ-10** — `.dockstore.yml` entry added for `parse_kraken2_reads` workflow (no `testParameterFiles` — placeholder paths not for CI execution)

## Non-Functional Requirements

- [ ] **REQ-11** — Task follows `tasks_metagenomics.wdl` conventions: `parameter_meta` with `description`/`patterns`/`category` for each input, task structure matches immediate neighbours (`classify_virnucpro_contigs`, `classify_reads_by_contig`)

## Out of Scope

- `build_taxonomy_db.py` WDL wrapper — DuckDB assumed pre-built and passed as input; wrapping the build step deferred
- Parquet output format — TSV matches repo conventions; pyarrow not in py3-bio image
- Multi-sample scatter workflow — not requested
- Combined end-to-end workflow chaining Kraken2 → parse → downstream — not requested

## Traceability

| REQ-ID | Phase | Plan |
|--------|-------|------|
| REQ-01 | Phase 2 | TBD |
| REQ-02 | Phase 2 | TBD |
| REQ-03 | Phase 2 | TBD |
| REQ-04 | Phase 2 | TBD |
| REQ-05 | Phase 2 | TBD |
| REQ-06 | Phase 2 | TBD |
| REQ-07 | Phase 2 | TBD |
| REQ-08 | Phase 2 | TBD |
| REQ-09 | Phase 2 | TBD |
| REQ-10 | Phase 2 | TBD |
| REQ-11 | Phase 2 | TBD |

---
*Last updated: 2026-04-01 — traceability updated by roadmap creation*
