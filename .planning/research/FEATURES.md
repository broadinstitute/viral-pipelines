# Feature Research

**Domain:** Centrifuger WDL integration — taxonomic classifier wrapping for viral-pipelines
**Researched:** 2026-04-01
**Confidence:** HIGH (CLI verified against official GitHub README + Genome Biology paper; WDL patterns verified against repo source)

---

## Centrifuger CLI and I/O Reference

This section is the authoritative reference for requirements definition.

### Centrifuger Binary Inputs

#### Required

| Input | Flag | Type | Notes |
|-------|------|------|-------|
| Index prefix | `-x FILE` | String (prefix path) | References `[prefix].1.cfr`, `.2.cfr`, `.3.cfr`, `.4.cfr` — four files |
| Reads (paired R1) | `-1 FILE` | FASTQ | Paired-end mode |
| Reads (paired R2) | `-2 FILE` | FASTQ | Paired-end mode |
| Reads (single-end) | `-u FILE` | FASTQ | Single-end mode; mutually exclusive with `-1/-2` |
| Reads (interleaved) | `-i FILE` | FASTQ | Interleaved paired-end; mutually exclusive with `-1/-2` and `-u` |

**Confidence: HIGH** — Verified from official GitHub README (mourisl/centrifuger).

The repo currently passes BAM files to Kraken2 (`reads_bam: File`) and calls `metagenomics kraken2` (viral-ngs Python wrapper), which handles BAM-to-FASTQ conversion internally. Centrifuger takes FASTQ directly; BAM-to-FASTQ conversion must be handled in the WDL command block explicitly, or by using the viral-ngs BAM conversion utilities already present in the repo.

#### Optional Classification Flags

| Flag | Type | Default | Purpose | Expose in WDL? |
|------|------|---------|---------|----------------|
| `-t INT` | Int | 1 | Number of threads | YES — tie to `cpu` in runtime |
| `-k INT` | Int | 1 | Report up to k distinct primary assignments per read pair | YES — optional input |
| `--min-hitlen INT` | Int | auto | Minimum hit length filter | YES — optional input |
| `--hitk-factor INT` | Int | 40 | Hit resolution factor | NO — advanced tuning, not needed for MVP |
| `--un STR` | String (prefix) | none | Output unclassified reads to `<prefix>_1/2.fq.gz` | NO — optional, adds output management complexity |
| `--cl STR` | String (prefix) | none | Output classified reads to `<prefix>_1/2.fq.gz` | NO — same reason |
| `--merge-readpair` | Boolean | false | Overlap merge and adapter trim for paired reads | NO — pre-trimming happens upstream |

**Not applicable to this use case:** `--barcode`, `--UMI`, `--read-format`, `--barcode-whitelist` — single-cell features, not relevant for bulk viral metagenomics.

**Confidence: HIGH** — Verified from official GitHub README.

### Centrifuger Index File Structure

The index is a set of four files: `[prefix].1.cfr`, `[prefix].2.cfr`, `[prefix].3.cfr`, `[prefix].4.cfr`. These are run-block compressed BWT indices (FM-index). For WDL purposes, these must be delivered as a tarball identical in convention to `kraken2_db_tgz` — a compressed archive decompressed in the command block.

**Pre-built index sizes (reference, for disk sizing):**
- RefSeq human + bacteria + archaea + virus: ~41 GB index, ~43 GB peak RAM
- GTDB r226: ~164 GB index
- NCBI core nt: ~212 GB index

These are 2–5x smaller than equivalent Kraken2 databases for the same reference set.

**Confidence: HIGH** — From official README and Genome Biology paper (PMC11046777).

### Centrifuger Output Files

#### Primary Classification Output (stdout → TSV)

Centrifuger writes classification results to **stdout by default**. The command must redirect to a file (e.g., `centrifuger ... > sample.centrifuger.tsv`).

| Column | Name | Description |
|--------|------|-------------|
| 1 | readID | Read identifier |
| 2 | seqID | Sequence ID where read classified |
| 3 | taxID | NCBI taxonomy ID |
| 4 | score | Weighted hit sum (higher = more specific hit) |
| 5 | 2ndBestScore | Score for next-best classification |
| 6 | hitLength | Base pairs matching genomic sequence |
| 7 | queryLength | Length of read (or combined mate pair length) |
| 8 | numMatches | Number of matching positions |

This file is analogous to Kraken2's `kraken2.reads.txt` (per-read output). In the existing repo, this is compressed with `pigz` and output as `.txt.gz`.

**Confidence: HIGH** — Verified from official README.

#### Kraken-Style Report (centrifuger-kreport)

`centrifuger-kreport` is a script bundled with Centrifuger that converts the per-read TSV to a Kraken-style summary report.

**Invocation:**
```
centrifuger-kreport -x <index_prefix> [OPTIONS] <centrifuger_output.tsv> > report.txt
```

| Option | Purpose |
|--------|---------|
| `-x INDEX` | Required — index prefix (used to pull taxonomy data) |
| `--no-lca` | Report count fractions at each taxon instead of LCA |
| `--show-zeros` | Include clades with zero reads |
| `--min-score INT` | Filter reads below alignment score |
| `--min-length INT` | Filter alignments below length threshold |
| `--report-score-data` | Add extra column with best classification score |

Output columns match the Kraken report format: `%reads | clade_reads | taxon_reads | rank_code | taxID | name`

This report is analogous to Kraken2's `kraken2.report.txt`. The existing `krona` task in `tasks_metagenomics.wdl` accepts kreport-format files via `--inputType kraken2`. Centrifuger kreport output is compatible with this same Krona task.

**Confidence: HIGH** — Verified from centrifuge-kreport source code review (format is identical to centrifuge kreport, which is the same Perl script with `--report-score-data` added).

#### Quantification Output (centrifuger-quant, optional)

`centrifuger-quant` produces an abundance/profiling table from the per-read TSV.

```
centrifuger-quant -x <index_prefix> -c <classification_tsv>
```

Output columns: `name | taxID | taxRank | genomeSize | numReads | numUniqueReads | abundance`

This is analogous to Bracken or Kraken2's summary report but includes genome-size-normalized abundance. Optional for the MVP — not needed to match the existing `classify_kraken2.wdl` feature set.

**Confidence: HIGH** — Verified from official README.

---

## Feature Landscape

### Table Stakes (Users Expect These)

Features that users expect if Centrifuger is positioned as a drop-in complement to Kraken2.

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| Per-read classification TSV | Core output of any classifier | LOW | centrifuger stdout → redirect to file; compress with pigz (same as kraken2 task) |
| Kraken-style summary report | Downstream tools (Bracken, Krona, aggregate reports) consume kreport format | LOW | `centrifuger-kreport -x index classification.tsv > report.txt`; same format as kraken2 report |
| Krona HTML visualization | Standard output of every classifier in this repo | LOW | Existing `krona` task in tasks_metagenomics.wdl accepts kreport format; `--inputType kraken2` works for centrifuger kreport |
| Paired-end FASTQ input | Most viral sequencing is paired-end | LOW | `-1 R1.fq -2 R2.fq`; if input is BAM, need BAM-to-FASTQ step in command block |
| Single-end FASTQ input | Some protocols produce SE reads | LOW | `-u reads.fq` |
| Tarball index as WDL File | Consistent with `kraken2_db_tgz` pattern; cloud storage requires single-file inputs | MEDIUM | Four `.cfr` files must be packed into a tarball; decompressed in command block |
| centrifuger_single.wdl standalone workflow | Users expect one-sample-at-a-time workflow for Terra scatter | LOW | Direct model of `classify_kraken2.wdl` |
| centrifuger_multi.wdl batch workflow | Amortize 40+ GB database load across samples | MEDIUM | Database extracted once outside scatter block; samples scatter over in-memory DB |
| Test input JSONs | Required for miniwdl check pass + Dockstore registration | LOW | Placeholder paths following established test/input/ pattern |
| Dockstore registration | All workflows in this repo are registered | LOW | `.dockstore.yml` entries following existing pattern |

### Differentiators (Centrifuger-Specific Capabilities)

Features that Centrifuger enables but Kraken2 does not, worth exposing in the WDL wrapper.

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Database-load-once multi-sample mode | Centrifuger's 43 GB index (vs Kraken2's ~100+ GB) makes in-memory DB feasible across more samples per node; reduces I/O cost | MEDIUM | `centrifuger_multi.wdl` design: extract index once, scatter classification over Array[File] samples within same VM |
| `centrifuger-quant` abundance table | Genome-size-normalized abundance output beyond what kreport provides; useful for comparing samples | MEDIUM | Optional second output from `centrifuger-quant -x index -c reads_tsv`; not required for MVP but clean to expose |
| Per-read sequence ID assignment | Centrifuger reports the specific sequence (seqID) that classified a read, not just taxID; enables strain-level resolution | LOW | Already in default output TSV column 2; no extra flags needed |

### Anti-Features (Commonly Requested, Often Problematic)

| Feature | Why Requested | Why Problematic | Alternative |
|---------|---------------|-----------------|-------------|
| BAM input as primary interface | Existing repo uses BAM for Kraken2 via viral-ngs wrapper | Centrifuger only accepts FASTQ; adding BAM support requires either a new viral-ngs wrapper (not available) or explicit BAM-to-FASTQ in the command block using samtools or similar — adds dependency and complexity | Accept FASTQ inputs directly; document BAM-to-FASTQ as an upstream step; the existing `tasks_read_utils.wdl` has conversion tools that can be called before the centrifuger task |
| Direct Krona generation in centrifuger task | Kraken2 task generates Krona in the same task | Centrifuger-kreport requires the index to be present (reads taxonomy from index files); if the index is already decompressed in the same task, this is feasible — but it adds significant disk pressure when the index is 40+ GB and the krona taxonomy DB is separate | Either: (a) generate kreport + krona in the same task while the index is decompressed, OR (b) output the kreport and let the caller chain to the existing standalone `krona` task. Option (a) matches kraken2 pattern but requires careful disk sizing. Recommend option (a) to match repo conventions. |
| `--un` / `--cl` unclassified/classified read outputs | Useful for human depletion pipeline | Adds two more output files with file-prefix naming that doesn't fit clean WDL output declaration; read filtering in this repo uses `filter_bam_to_taxa` (Kraken2-based), not centrifuger | Keep centrifuger task focused on classification output only; route read filtering through existing infrastructure |
| centrifuger-build index construction task | Building custom databases | Index build requires nodes.dmp, names.dmp, reference FASTAs, seqID-to-taxID tables — a separate workflow with very different resource requirements (RAM-intensive build); outside the scope of classification wrapping | Provide pre-built index tarballs as WDL File inputs; defer `centrifuger-build` to a future separate workflow |

---

## Feature Dependencies

```
[BAM input]
    └──requires──> [BAM-to-FASTQ conversion] (upstream step, not in centrifuger task)
                       └──requires──> [centrifuger task]

[centrifuger task]
    └──produces──> [per-read TSV]
                       └──requires──> [centrifuger-kreport] (needs index + TSV)
                                          └──produces──> [kreport.txt]
                                                             └──requires──> [krona task] (existing)
                                                                                └──produces──> [krona.html]

[centrifuger task]
    └──produces──> [per-read TSV]
                       └──requires──> [centrifuger-quant] (needs index + TSV)
                                          └──produces──> [abundance.tsv] (optional/differentiator)

[centrifuger_multi.wdl]
    └──requires──> [centrifuger task] (with Array[File] scatter pattern)
    └──depends-on──> [index extracted once before scatter]
```

### Dependency Notes

- **kreport requires index:** `centrifuger-kreport -x <index>` reads taxonomy data from the index files. This means the kreport generation step needs the index to still be decompressed, creating a strong argument for generating kreport + krona within the same task (while index is on-disk), not as a separate downstream task.
- **krona task is reusable:** The existing `krona` task in `tasks_metagenomics.wdl` already accepts kreport format files (`--inputType kraken2`). If kreport is output from the centrifuger task, callers can chain to the existing `krona` task without modification.
- **centrifuger_multi.wdl requires scatter-after-db-extract:** The multi-sample workflow cannot use a simple `scatter(sample in samples) { call centrifuger }` pattern because that would decompress the 40+ GB index separately per sample. The pattern used in this repo's Kallisto workflows (see `classify_kallisto_multi.wdl`) should be studied as a model for in-task batching.

---

## WDL Input/Output Table (Required vs Optional)

This is the primary deliverable for requirements definition.

### centrifuger WDL Task Inputs

#### Required Inputs

| WDL Name | Type | Centrifuger Flag | Notes |
|----------|------|-----------------|-------|
| `reads_fq1` | `File` | `-1` | Paired-end R1 FASTQ (may be .gz) |
| `reads_fq2` | `File?` | `-2` | Paired-end R2 FASTQ; optional (absent = single-end mode) |
| `reads_fq_se` | `File?` | `-u` | Single-end FASTQ; mutually exclusive with paired inputs |
| `centrifuger_db_tgz` | `File` | `-x` (after decompression) | Tarball containing `[prefix].1.cfr` through `[prefix].4.cfr` |
| `krona_taxonomy_db_tgz` | `File` | (krona step) | Krona taxonomy.tab — same as kraken2 task |

**Note on input design:** The WDL task should accept either (paired R1 + optional R2) OR (SE reads). A common pattern is `File reads_fq1, File? reads_fq2` where the command block detects whether R2 is defined to switch between `-1/-2` and `-u` mode. This is simpler than a conditional input with separate modes.

#### Optional Inputs with Defaults

| WDL Name | Type | Default | Centrifuger Flag | Notes |
|----------|------|---------|-----------------|-------|
| `min_hit_len` | `Int?` | (centrifuger auto) | `--min-hitlen` | Leave unset to use centrifuger's heuristic |
| `report_top_k` | `Int?` | 1 | `-k` | Expose for users who want multiple assignments |
| `threads` | `Int` | 16 | `-t` | Match CPU allocation in runtime block |
| `machine_mem_gb` | `Int` | 50 | (runtime memory) | Must exceed peak RAM; 43 GB for RefSeq index + overhead |
| `docker` | `String` | `quay.io/biocontainers/centrifuger:1.1.0--h4349b82_0` | (runtime docker) | BioContainers tag — must be verified before use |

#### Runtime Block Target

| Resource | Value | Rationale |
|----------|-------|-----------|
| `memory` | `50 GB` | 43 GB peak RAM for RefSeq index; 50 GB provides headroom |
| `cpu` | 16 | Centrifuger scales with threads; match kraken2 task |
| `disks` | auto-calculated | Index tarball + decompressed index (2x) + reads + output; floor at 375 GB for SSD allocation |
| `dx_instance_type` | `mem3_ssd1_v2_x8` | Matches kraken2 task; memory-bound workload |
| `preemptible` | 2 | Consistent with kraken2 task |

### centrifuger WDL Task Outputs

| WDL Name | Type | Source | Description |
|----------|------|--------|-------------|
| `centrifuger_reads_tsv_gz` | `File` | stdout redirect + pigz | Per-read classification TSV (8 columns) compressed; analogous to `kraken2_reads_report` |
| `centrifuger_summary_report` | `File` | `centrifuger-kreport` stdout redirect | Kraken-style summary report; analogous to `kraken2_summary_report` |
| `krona_report_html` | `File` | existing `metagenomics krona` call | Krona HTML visualization; analogous to `krona_report_html` in kraken2 task |
| `max_ram_gb` | `Int` | `/sys/fs/cgroup/memory.peak` | Peak RAM consumed; standard monitoring output in this repo |
| `viralngs_version` | `String` | `metagenomics --version` | Tool version tracking; standard in this repo |

**Note:** `centrifuger_summary_report` and `krona_report_html` are generated within the same task while the index is still decompressed (kreport requires index on-disk). This matches the Kraken2 task pattern precisely.

---

## Krona Integration Path

The existing `krona` task in `tasks_metagenomics.wdl` accepts kreport-format files and takes `--inputType kraken2`. Since `centrifuger-kreport` produces output in the identical Kraken-style format, the integration path is:

1. Run `centrifuger-kreport -x $DB_DIR/centrifuger ... > sample.centrifuger.report.txt` inside the centrifuger task (while the decompressed index is still present at `$DB_DIR/centrifuger`)
2. Call `metagenomics krona sample.centrifuger.report.txt ... --inputType kraken2` — the same viral-ngs Krona wrapper already used for Kraken2, KrakenUniq, etc.

**This is a direct reuse of existing infrastructure.** No new Krona conversion tooling is needed. The `centrifuger_summary_report` output type is compatible with `reports.aggregate_metagenomics_reports` in `tasks_reports.wdl` for cross-sample summary, same as Kraken2.

**Confidence: HIGH** — The `krona` task's `input_type` parameter passes `--inputType` to the viral-ngs `metagenomics krona` command; kreport format is the same for Kraken2 and Centrifuger.

---

## Complexity Comparison: Centrifuger vs Kraken2 WDL Wrapping

| Dimension | Kraken2 WDL Task | Centrifuger WDL Task | Delta |
|-----------|-----------------|---------------------|-------|
| Input read format | BAM (via viral-ngs wrapper) | FASTQ (direct binary call) | BAM-to-FASTQ required upstream or in task |
| Index format | `hash.k2d + opts.k2d + taxo.k2d` tarball | `[prefix].{1,2,3,4}.cfr` tarball | Same tarball pattern, different internal files |
| Index RAM at peak | ~100+ GB (standard RefSeq) | ~43 GB (same RefSeq) | **Centrifuger is 2–3x more memory-efficient** |
| Classification output | Per-read TXT (redirected to .txt.gz) | Per-read TSV (redirected to .tsv.gz) | Functionally identical pattern |
| Summary report | `kraken2 --report` flag during classification | Separate `centrifuger-kreport -x index` step | **Extra step requiring index to still be on-disk** |
| Krona generation | Calls `metagenomics krona --inputType kraken2` | Same call, same format | No difference |
| Docker image | `quay.io/broadinstitute/viral-ngs:3.0.10-classify` | `quay.io/biocontainers/centrifuger:1.1.0--...` | Different image; no viral-ngs wrapper available |
| Viral-ngs wrapper available | YES — `metagenomics kraken2` handles BAM conversion | NO — must call `centrifuger` binary directly | **Must write raw bash, not use viral-ngs wrapper** |
| Multi-sample disk amortization | Kraken2 index too large for easy in-memory multi-sample | Centrifuger index small enough to hold in RAM across samples | **Centrifuger genuinely benefits from multi-sample mode** |
| WDL wrapping complexity | LOW — existing working task as reference | MEDIUM — no viral-ngs wrapper; FASTQ input; kreport as separate step | Slightly more complex command block |

**Overall assessment:** Centrifuger WDL wrapping is **marginally more complex** than Kraken2 due to: (1) no viral-ngs Python wrapper to lean on (must call binary directly and handle BAM→FASTQ if BAM inputs are required), (2) kreport generation requiring the index to remain decompressed. The fundamental WDL patterns (tarball input, runtime sizing, output naming) are identical. Estimated effort is 1.2–1.5x the Kraken2 task development effort.

---

## MVP Definition

### v3.0 Launch With

- [x] `centrifuger` WDL task in `tasks_metagenomics.wdl` — core classifier wrapping; paired-end FASTQ input, tarball index, classification TSV + kreport + Krona outputs
- [x] `centrifuger_single.wdl` — standalone workflow for one sample; direct model of `classify_kraken2.wdl`
- [x] `centrifuger_multi.wdl` — batch workflow; Array[File] samples, database extracted once, scatter classification
- [x] Test input JSONs for both workflows in `test/input/WDL/miniwdl-local/`
- [x] `.dockstore.yml` entries for both workflows
- [x] All files pass `miniwdl check`

### Add After Validation (v3.x)

- [ ] `centrifuger-quant` output — abundance table with genome-size normalization; add as optional second output from the task when users request it
- [ ] Real test data files for `miniwdl run` execution — deferred per established repo pattern
- [ ] BAM input mode — if users require BAM as primary input, add BAM-to-FASTQ conversion step within the task using `samtools fastq`

### Future Consideration (v4+)

- [ ] `centrifuger-build` index construction workflow — separate resource profile, out of scope for classification milestone
- [ ] Centrifuger integration with `classify_single.wdl` / `classify_multi.wdl` — once single-tool workflows validate, consider adding Centrifuger as a parallel classification arm in the multi-tool workflows

---

## Feature Prioritization Matrix

| Feature | User Value | Implementation Cost | Priority |
|---------|------------|---------------------|----------|
| centrifuger task (classification + kreport + krona) | HIGH | MEDIUM | P1 |
| centrifuger_single.wdl | HIGH | LOW | P1 |
| centrifuger_multi.wdl (db-load-once) | HIGH | MEDIUM | P1 |
| Test input JSONs | HIGH (CI/Dockstore) | LOW | P1 |
| Dockstore entries | HIGH (discoverability) | LOW | P1 |
| centrifuger-quant output | MEDIUM | LOW | P2 |
| BAM input mode | MEDIUM | MEDIUM | P2 |
| centrifuger-build workflow | LOW | HIGH | P3 |

---

## Sources

- [Centrifuger GitHub README — mourisl/centrifuger](https://github.com/mourisl/centrifuger) — HIGH confidence; CLI flags, index format, output format
- [Centrifuger Genome Biology paper (PMC11046777)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11046777/) — HIGH confidence; memory benchmarks, sensitivity comparison with Kraken2
- [centrifuge-kreport source (DaehwanKimLab/centrifuge)](https://github.com/DaehwanKimLab/centrifuge/blob/master/centrifuge-kreport) — HIGH confidence; kreport output column format (centrifuger-kreport is derived from same script)
- [Bioconda centrifuger recipe](https://bioconda.github.io/recipes/centrifuger/README.html) — MEDIUM confidence; version 1.1.0 availability, quay.io/biocontainers/centrifuger tag pattern
- viral-pipelines `tasks_metagenomics.wdl` (local) — HIGH confidence; established patterns for kraken2 task, krona task, runtime blocks
- viral-pipelines `classify_kraken2.wdl`, `classify_multi.wdl` (local) — HIGH confidence; workflow patterns to follow

---

*Feature research for: Centrifuger WDL integration (viral-pipelines v3.0 milestone)*
*Researched: 2026-04-01*
