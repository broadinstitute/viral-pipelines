# Requirements: v3.0 Centrifuger Taxonomic Classification WDL

**Defined:** 2026-04-02
**Core Value:** Adds Centrifuger as a first-class taxonomic classifier in viral-pipelines, with single-sample and multi-sample workflows that amortize the 200+ GB database load across all samples on a node.

## v1 Requirements

### Centrifuger Task

- [x] **CFGR-01**: WDL task `centrifuger` in `tasks_metagenomics.wdl` that wraps the `centrifuger` binary
  - Accepts pre-built index via `Directory centrifuger_db` and `String db_name` (index prefix)
  - Accepts optional FASTQ inputs: `File? reads_fastq1`, `File? reads_fastq2` (paired-end), `File? reads_fastq_unpaired` (single-end/long)
  - Accepts optional BAM input: `File? reads_bam` — converted to FASTQ internally via `picard SamToFastq`
  - Accepts `String sample_name` for output file naming
  - Outputs `File classification_tsv` — per-read classification results
  - Outputs `File kreport` — Kraken-style summary report (via `centrifuger-kreport`)
  - ~~Outputs `File krona_html`~~ — **deferred**: KronaTools absent from `ghcr.io/broadinstitute/docker-centrifuger:1.0.0`; will be addressed in a future phase when Docker image is updated
  - Outputs `File centrifuger_log` — captured stderr/runtime log
  - Runtime sized for 200+ GB index: `cpu: 8`, `memory: "240 GB"`, `disks: "local-disk 400 SSD"`
  - Docker: `ghcr.io/broadinstitute/docker-centrifuger:1.0.0`

- [x] **CFGR-02**: BAM→FASTQ conversion via `picard SamToFastq`
  - Paired-end detection: run `picard SamToFastq` with both `FASTQ=R1.fq` and `SECOND_END_FASTQ=R2.fq`; check if R2 is non-empty to determine paired vs. single-end
  - Paired-end reads: pass `-1 R1.fq -2 R2.fq` to centrifuger; read IDs must carry `/1` and `/2` suffixes (use Picard `READ1_SUFFIX=/1 READ2_SUFFIX=/2`)
  - Single-end reads: pass `-u reads.fq` to centrifuger (note: `-r` is a centrifuger-build flag for reference sequences, not for classification reads)
  - At least one of `reads_bam`, `reads_fastq1`, or `reads_fastq_unpaired` must be provided

### Single-Sample Workflow

- [x] **CFGR-03**: Standalone workflow `pipes/WDL/workflows/centrifuger_single.wdl`
  - Scalar inputs (one sample per workflow call — Terra scatter-friendly)
  - Imports `tasks_metagenomics.wdl` and calls task with alias `as run_centrifuger`
  - `meta { allowNestedInputs: true }` for Terra-compatible test JSON
  - Passes `miniwdl check` validation

### Multi-Sample Workflow

- [x] **CFGR-04**: Standalone workflow `pipes/WDL/workflows/centrifuger_multi.wdl`
  - Accepts `Array[File]+ reads_bams` OR `Array[File]+ reads_fastqs` as batch input
  - Uses a single WDL task call — sample loop lives inside the task command block (krakenuniq bash-loop pattern; NOT WDL scatter)
  - Database loaded once per task call; all samples classified in sequence on same node
  - Outputs parallel arrays for all four file types (one entry per sample)
  - Runtime sized for full cohort on 1–4 nodes: `cpu: 16`, `memory: "480 GB"`, `disks: "local-disk 600 SSD"`
  - Passes `miniwdl check` validation

### Infrastructure

- [x] **CFGR-05**: Test input JSON `test/input/WDL/miniwdl-local/test_inputs-centrifuger_single-local.json`
  - Placeholder paths for index directory and FASTQ/BAM inputs
  - Workflow-level input keys; follows v1.0–v2.0 test file conventions

- [x] **CFGR-06**: Test input JSON `test/input/WDL/miniwdl-local/test_inputs-centrifuger_multi-local.json`
  - Placeholder paths for Array inputs
  - Workflow-level input keys; follows v1.0–v2.0 test file conventions

- [x] **CFGR-07**: Dockstore registration entries in `.dockstore.yml` for both `centrifuger_single.wdl` and `centrifuger_multi.wdl`
  - Subclass: WDL
  - No `testParameterFiles` (placeholder paths not CI-runnable)

## v2 Requirements — join_read_classifications

### WDL Task

- [ ] **JRC-01**: WDL task `join_read_classifications` in `tasks_metagenomics.wdl`
  - Accepts 4 optional File inputs: `File? kallisto_summary` (Parquet), `File? kraken2_reads` (Parquet), `File? vnp_reads` (Parquet), `File? genomad_virus_summary` (TSV)
  - Accepts required `String sample_id` — filters Kallisto/K2 tables, stamps SAMPLE_ID on output
  - Embeds full join_read_classifications.py logic as inline `python3<<CODE` heredoc
  - 4-way DuckDB join: Kallisto FULL OUTER VNP on READ_ID, then FULL OUTER K2 with /1|/2 suffix strip, then LEFT JOIN geNomad via VNP_CONTIG_ID
  - Outputs `File classifications_parquet` — ZSTD-compressed Parquet
  - Docker: `quay.io/broadinstitute/py3-bio:0.1.3` with `pip install duckdb`
  - Runtime: 16 GB memory, 1 CPU, disk autoscaled from input sizes

### Standalone Workflow

- [ ] **JRC-02**: Standalone workflow `pipes/WDL/workflows/join_read_classifications.wdl`
  - Imports `tasks_metagenomics.wdl`, calls task with alias `as join_reads`
  - `meta { allowNestedInputs: true }` for Terra-compatible test JSON
  - Passes `miniwdl check` validation

### Infrastructure

- [ ] **JRC-03**: Test input JSON `test/input/WDL/miniwdl-local/test_inputs-join_read_classifications-local.json`
  - Placeholder paths for all 4 optional File inputs plus sample_id
  - Workflow-level input keys; follows existing test file conventions

- [ ] **JRC-04**: Dockstore registration entry in `.dockstore.yml` for `join_read_classifications.wdl`
  - Subclass: WDL
  - No `testParameterFiles` (placeholder paths not CI-runnable)

## Out of Scope

| Feature | Reason |
|---------|--------|
| Centrifuger index build task | Use pre-built index; build is a separate concern |
| Combined end-to-end pipeline chaining Centrifuger with downstream steps | Standalone-first per project pattern |
| Custom Docker image build | Using `ghcr.io/broadinstitute/docker-centrifuger:1.0.0` (Broad-maintained) |
| Real test data for `miniwdl run` execution | Placeholder JSONs acceptable; real data deferred |
| WDL scatter for multi-sample | Defeats the purpose of amortizing database load; bash-loop is the correct pattern |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| CFGR-01 | Phase 4 | Complete |
| CFGR-02 | Phase 4 | Complete |
| CFGR-03 | Phase 5 | Complete |
| CFGR-04 | Phase 5 | Complete |
| CFGR-05 | Phase 6 | Complete |
| CFGR-06 | Phase 6 | Complete |
| CFGR-07 | Phase 6 | Complete |
| JRC-01 | Phase 7 | Planned |
| JRC-02 | Phase 7 | Planned |
| JRC-03 | Phase 7 | Planned |
| JRC-04 | Phase 7 | Planned |

**Coverage:**
- v1 requirements: 7 total, 7 complete
- v2 requirements: 4 total, 0 complete
- Unmapped: 0

---
*Requirements defined: 2026-04-02*
*Last updated: 2026-04-02 — Phase 7 requirements added*
