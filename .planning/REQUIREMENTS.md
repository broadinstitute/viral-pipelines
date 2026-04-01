# Requirements

## Overview

Add two WDL tasks and standalone workflows to wrap the `classify_virnucpro_contigs.py` and `classify_reads_by_contig.py` scripts from `respiratory-surveillance/bin/`. Tasks go in `tasks_metagenomics.wdl`. Both tasks embed their Python scripts as inline `python3<<CODE` heredocs, consistent with the repo-wide pattern.

## Functional Requirements

### Task: classify_virnucpro_contigs

- [x] **REQ-01** — Task `classify_virnucpro_contigs` added to `pipes/WDL/tasks/tasks_metagenomics.wdl` after the existing `classify_virnucpro` task
- [x] **REQ-02** — Inputs: `File virnucpro_scores_tsv` (TSV of chunk scores from VirNucPro), optional tuning params (`min_viral_prop`, `min_nonviral_prop`, `min_chunks`, `id_col`, `id_pattern`)
- [x] **REQ-03** — Output: `File contig_classifications` (TSV with columns: ID, call, tier, weighted_delta, n_chunks, etc.)
- [x] **REQ-04** — Script embedded inline as `python3<<CODE` heredoc — no external file dependency
- [x] **REQ-05** — Docker: `quay.io/broadinstitute/py3-bio:0.1.3` (pandas is pre-installed)
- [x] **REQ-06** — Runtime: 2 GB memory, 1 CPU, 50 GB HDD; includes `# TES` disk annotation

### Task: classify_reads_by_contig

- [x] **REQ-07** — Task `classify_reads_by_contig` added to `pipes/WDL/tasks/tasks_metagenomics.wdl` after `classify_virnucpro_contigs`
- [x] **REQ-08** — Inputs: `File paf_file` (gzipped PAF from BamToPAF), `File contig_classifications` (from `classify_virnucpro_contigs`), optional filter params (`min_mapq`, `min_identity`, `min_query_cov`)
- [x] **REQ-09** — Output: `File read_classifications` (TSV with read_id, call, tier, mapping metrics)
- [x] **REQ-10** — Script embedded inline as `python3<<CODE` heredoc
- [x] **REQ-11** — Docker: `quay.io/broadinstitute/py3-bio:0.1.3`; command preamble runs `pip install duckdb --quiet --no-cache-dir` before the Python block; `set -e` guards the install
- [x] **REQ-12** — Runtime: 4 GB memory, 2 CPU, 100 GB HDD; includes `# TES` disk annotation

### Standalone Workflows

- [x] **REQ-13** — `pipes/WDL/workflows/classify_virnucpro_contigs.wdl` — thin workflow wrapping `classify_virnucpro_contigs` task; `meta` block with description/author/email/allowNestedInputs
- [x] **REQ-14** — `pipes/WDL/workflows/classify_reads_by_contig.wdl` — thin workflow wrapping `classify_reads_by_contig` task; same meta block pattern

### Validation

- [x] **REQ-15** — Both workflows pass `miniwdl check`
- [x] **REQ-16** — Test input JSON files in `test/input/WDL/miniwdl-local/` for both workflows
- [x] **REQ-17** — `.dockstore.yml` entries added for both new workflows

## Non-Functional Requirements

- [x] **REQ-18** — Tasks follow existing `tasks_metagenomics.wdl` conventions: `parameter_meta` with `description`/`patterns`/`category`, `disks`/`disk` dual runtime annotation, `preemptible: 2` default
- [x] **REQ-19** — `requirements-modules.txt` does not need a new entry (py3-bio already listed)

## Out of Scope

- Combined end-to-end workflow chaining all three steps — user chose standalone only
- Custom Docker image build — pip install at runtime is sufficient
- Parquet output format — TSV default matches repo conventions
- Multi-sample scatter workflows — not requested

---
*Last updated: 2026-03-31 after initialization*
