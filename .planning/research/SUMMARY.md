# Project Research Summary

**Project:** Centrifuger WDL Integration — viral-pipelines v3.0
**Domain:** WDL task/workflow wrapping of a large-database FM-index taxonomic classifier (Centrifuger) on Terra/GCP/DNAnexus
**Researched:** 2026-04-01
**Confidence:** HIGH

## Executive Summary

This project adds Centrifuger taxonomic classification to the viral-pipelines WDL library as the v3.0 milestone. Centrifuger is a successor to Centrifuge (not Kraken2) that uses an FM-index rather than a full k-mer hash, achieving 2-3x lower RAM requirements for the same reference set. The recommended approach is a single WDL task in `tasks_metagenomics.wdl` that accepts `Array[File]+` BAM inputs, converts them to FASTQ via samtools, runs classification in a bash loop against a single loaded index, and emits per-read TSV and Kraken-style summary report outputs. Two thin wrapper workflows (`centrifuger_single.wdl`, `centrifuger_multi.wdl`) follow the exact pattern already established by the VirNucPro and KrakenUniq milestones.

The correct Docker image is `quay.io/biocontainers/centrifuger:1.1.0--hf426362_0` from BioContainers — no Broad image includes the `centrifuger` binary, and no custom image build is needed. All outputs (per-read TSV, Kraken-style kreport, Krona HTML) are compatible with the existing metagenomics aggregation and visualization tooling. The kreport generation step (`centrifuger-kreport`) requires the index to remain decompressed on disk, making in-task generation the only practical approach.

The single critical architectural decision is that `centrifuger_multi.wdl` must NOT scatter over samples. Scattering would localize the 200+ GB index independently per scatter shard, eliminating the amortization goal that justifies a multi workflow. Instead, the task must accept the full `Array[File]+` and iterate samples in a bash loop — the same pattern already proven in the `krakenuniq` task at lines 83-84 of `tasks_metagenomics.wdl`. Pitfall risk is moderate: index format confusion (`.cfr` vs `.cf`), disk allocation math, and RAM sizing are all addressable with explicit guards and overridable runtime inputs.

## Key Findings

### Recommended Stack

The stack is fully determined by direct codebase inspection and external registry confirmation. All tasks that run Python scripts use `quay.io/broadinstitute/py3-bio:0.1.3` with the `python3<<CODE ... CODE` inline heredoc pattern. The Centrifuger task is a binary-invocation task (not a Python task), so it uses the BioContainers image instead. The `samtools` binary needed for BAM-to-FASTQ conversion is available in the BioContainers Centrifuger image via conda dependencies. No custom Docker builds are required for any deliverable in v3.0.

**Core technologies:**
- `quay.io/biocontainers/centrifuger:1.1.0--hf426362_0`: Centrifuger classifier — only image with the `centrifuger` binary; BioContainers auto-published from Bioconda, no custom build needed
- `centrifuger` + `centrifuger-kreport`: Classification and report generation — bundled together in the BioContainers image; kreport requires index still on-disk so must be called within the same task
- `samtools fastq`: BAM-to-FASTQ conversion — required because Centrifuger has no native BAM support; available in the BioContainers image via conda
- `quay.io/broadinstitute/py3-bio:0.1.3` + `pip install duckdb`: Python analysis tasks (v1.0-v2.0) — established codebase pattern; duckdb installed at runtime since it is absent from py3-bio

### Expected Features

Full feature research is in `.planning/research/FEATURES.md`.

**Must have (table stakes):**
- Per-read classification TSV compressed with pigz — core output; analogous to `kraken2_reads_report`
- Kraken-style summary report via `centrifuger-kreport` — required for downstream aggregation (`aggregate_metagenomics_reports`) and Krona
- Krona HTML visualization — all classifiers in this repo emit Krona; reuses existing `metagenomics krona --inputType kraken2` call
- Paired-end FASTQ input (via BAM-to-FASTQ conversion) — dominant input type in the repo
- Tarball index as WDL `File` — consistent with `kraken2_db_tgz` pattern; cloud storage requires single-file inputs
- `centrifuger_single.wdl` standalone workflow — one sample per Terra submission
- `centrifuger_multi.wdl` batch workflow — amortize 200+ GB database load across samples on one node
- Test input JSONs and Dockstore registration entries — required for `miniwdl check` pass and workflow discoverability

**Should have (differentiators):**
- `centrifuger-quant` abundance table — genome-size-normalized abundance output; adds value over kreport alone; optional second output
- Per-read seqID assignment — Centrifuger reports the specific sequence a read classified to (not just taxID), enabling strain-level resolution; already in default output, no extra flags

**Defer (v2+):**
- `centrifuger-build` index construction workflow — separate resource profile, high complexity, not a classification concern
- BAM input as primary interface (without samtools conversion) — requires explicit samtools step; acceptable to document this as an upstream caller responsibility
- Integration with `classify_single.wdl` / `classify_multi.wdl` multi-tool workflows — once single-tool workflows validate

### Architecture Approach

The architecture is additive: one new task appended to `tasks_metagenomics.wdl` (after line 2353), plus two new workflow files and two new test JSONs, plus two new `.dockstore.yml` entries. The `centrifuger` task uses the `krakenuniq` bash-loop pattern: `Array[File]+` input, single DB decompression before the loop, per-sample classification inside the loop, `glob()` output collection. The `centrifuger_single.wdl` wrapper passes `reads_bams = [reads_bam]` (single-element array) to the same task — no separate task definition is needed. The `centrifuger_multi.wdl` passes the full array directly with no scatter.

**Major components:**
1. `task centrifuger` (new, in `tasks_metagenomics.wdl`) — decompress index once, BAM-to-FASTQ loop, classify each sample, generate kreport and krona per-sample, emit Array outputs via glob
2. `centrifuger_single.wdl` (new) — thin wrapper; scalar `File reads_bam` wrapped as `[reads_bam]`; passes through task outputs
3. `centrifuger_multi.wdl` (new) — thin wrapper; `Array[File]+` passed directly; no scatter; documents batching rationale in meta block
4. Test JSONs and `.dockstore.yml` entries — required scaffolding; modeled on `classify_virnucpro_single` and `classify_virnucpro_multi` precedents

### Critical Pitfalls

Full pitfall research is in `.planning/research/PITFALLS.md`.

1. **Scatter in `centrifuger_multi.wdl`** — each Terra scatter shard is an independent VM; using `scatter(bam in reads_bams) { call centrifuger }` downloads and decompresses the 200+ GB index N times. Prevention: implement as single-task bash loop (krakenuniq pattern), never scatter.

2. **Index format confusion (`.cfr` vs `.cf`)** — Centrifuge and Centrifuger use incompatible binary index formats despite near-identical names. Prevention: add `parameter_meta` `patterns: ["*.cfr"]` and a runtime extension check guard in the command block (`ls $DB_DIR/*.cfr || exit 1`).

3. **BAM passed directly to Centrifuger** — Centrifuger accepts only FASTQ; no native BAM support. Prevention: explicit `samtools fastq` (or `samtools collate | samtools fastq`) step in the command block before calling `centrifuger`.

4. **Disk under-allocation for 200+ GB tarball** — tarball + decompressed index + BAM expansion can exceed 1500 GB for large databases. Prevention: use the `ceil((2 * size(db, "GB") + 10 * length(bams) + 50) / 750.0) * 750` formula with 750 GB floor, matching the `kraken2` task pattern exactly.

5. **RAM ceiling mismatch** — Centrifuger loads the FM-index entirely into RAM (no partial loading); a 200+ GB NT-scale database requires 240+ GB RAM. Prevention: expose `machine_mem_gb` as an overridable input with a documented default; do not copy `machine_mem_gb = 90` from the Kraken2 task without adjustment.

6. **Executor content duplication (process pitfall)** — appending to large WDL files without a pre-check has caused duplicate task bodies in prior milestones. Prevention: executor plan must include `grep -c "task centrifuger" tasks_metagenomics.wdl` check returning 0 before writing.

## Implications for Roadmap

Based on research, suggested phase structure:

### Phase 1: Core Task Implementation
**Rationale:** The task is the dependency for both workflow wrappers; it must pass `miniwdl check` before any workflow file is written. BAM conversion, kreport generation, Krona integration, and resource sizing are all resolved at this layer. This is the highest-complexity deliverable.
**Delivers:** `task centrifuger` in `tasks_metagenomics.wdl` — classification, BAM-to-FASTQ conversion, kreport, Krona HTML, per-read TSV output, max_ram_gb capture, disk autoscaling formula.
**Addresses:** Per-read TSV, kreport, Krona (all table-stakes features); `machine_mem_gb` and disk autoscaling (resource pitfalls).
**Avoids:** BAM passthrough bug, RAM OOM, disk OOM, output format mismatch.

### Phase 2: Single and Multi Workflow Wrappers
**Rationale:** Purely structural once the task exists; both wrappers are thin pass-throughs with no logic. The architectural decision (no scatter in multi) is already captured in the task design. Risk is low — this step is mechanical.
**Delivers:** `centrifuger_single.wdl` and `centrifuger_multi.wdl` — both pass `miniwdl check`; both include `allowNestedInputs: true` and correct call alias where needed.
**Uses:** BioContainers Centrifuger image (from STACK.md); krakenuniq-model array input pattern (from ARCHITECTURE.md).
**Avoids:** Scatter anti-pattern in multi workflow; call alias collision; missing `allowNestedInputs`.

### Phase 3: Test JSONs and Dockstore Registration
**Rationale:** Scaffolding that unlocks CI validation and Dockstore discoverability; low complexity, no implementation decisions. Follows the established test input pattern exactly.
**Delivers:** `centrifuger_single.json`, `centrifuger_multi.json` in `test/input/WDL/miniwdl-local/`; two new entries in `.dockstore.yml`. All placeholder paths, no real test data (consistent with v1.0-v2.0 precedent).
**Avoids:** Adding `testParameterFiles` to `.dockstore.yml` (placeholder JSONs cause Dockstore CI failures).

### Phase Ordering Rationale

- Task must precede workflows because workflow import resolves at `miniwdl check` time; a failing task import obscures the root error.
- Dockstore entries should come last to avoid registering workflows before they pass validation.
- BAM-to-FASTQ conversion, kreport, and Krona generation are all concentrated in Phase 1 because all three require the index to be present on disk in the same execution context — separating them across phases is architecturally unsound.
- Resource sizing (disk formula, RAM default) must be finalized in Phase 1 before the task can be tested; they are not scaffolding.

### Research Flags

Phases with well-documented patterns (no additional research needed):
- **Phase 2 (workflow wrappers):** Direct models exist in `classify_virnucpro_single.wdl` and `classify_virnucpro_multi.wdl`. No unknowns.
- **Phase 3 (test JSONs / Dockstore):** Pattern is fully established; no research required.

Phases that may benefit from spot-checks before execution:
- **Phase 1 (task implementation):** The BioContainers image tag `1.1.0--hf426362_0` should be confirmed against the live registry before writing the `docker` declaration. STACK.md marks this HIGH confidence but notes it requires registry lookup. Also confirm `samtools` is available in the BioContainers image before committing BAM conversion logic.

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | All images and patterns verified via direct codebase inspection; BioContainers tag confirmed via TRS API |
| Features | HIGH | CLI flags verified against official GitHub README and Genome Biology paper; WDL I/O table derived from codebase patterns |
| Architecture | HIGH | All structural patterns confirmed from direct codebase inspection (krakenuniq, classify_virnucpro, kraken2 tasks); insertion point confirmed (line 2353) |
| Pitfalls | HIGH | Critical pitfalls grounded in codebase precedents (RETROSPECTIVE.md lessons), Terra/Cromwell scatter semantics (documented), and Centrifuger paper memory numbers |

**Overall confidence:** HIGH

### Gaps to Address

- **Exact BioContainers image tag:** STACK.md reports `1.1.0--hf426362_0` as confirmed via TRS API. Executor should verify this tag is still pullable before writing the `docker` declaration. Fallback: check `quay.io/repository/biocontainers/centrifuger?tab=tags` directly.
- **`samtools` availability in BioContainers Centrifuger image:** Assumed present via conda dependencies; not independently verified. If absent, BAM-to-FASTQ must be handled upstream or a different image must be selected. Executor should confirm with `docker run quay.io/biocontainers/centrifuger:1.1.0--hf426362_0 which samtools` before committing conversion logic.
- **RAM default for production database scope:** ARCHITECTURE.md recommends 240 GB default for NT-scale; FEATURES.md suggests 50 GB for RefSeq-scope. The correct default depends on which pre-built index users will be directed to. This should be resolved with the team before Phase 1 is finalized.
- **`duckdb` absence from py3-bio image:** STACK.md marks this MEDIUM confidence — absence confirmed by lack of any import, but py3-bio 0.1.3 was not independently inspected. If `duckdb` is actually present, the `pip install duckdb` step can be removed from the `classify_reads_by_contig` task.

## Sources

### Primary (HIGH confidence)
- `pipes/WDL/tasks/tasks_metagenomics.wdl` (direct inspection) — krakenuniq/kraken2 task patterns, disk formula, runtime block structure, glob output pattern
- `pipes/WDL/workflows/classify_virnucpro_single.wdl`, `classify_virnucpro_multi.wdl` (direct inspection) — thin wrapper pattern, `allowNestedInputs`, single-element array wrap
- `.planning/PROJECT.md` (direct inspection) — key decisions, v3.0 milestone requirements, out-of-scope items
- `.planning/RETROSPECTIVE.md` (direct inspection) — process pitfalls from prior milestones
- [mourisl/centrifuger GitHub README](https://github.com/mourisl/centrifuger) — CLI flags, index format, output columns, binary names
- [Centrifuger Genome Biology 2024 (PMC11046777)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11046777/) — memory benchmarks, sensitivity comparison
- [BioContainers TRS API](https://api.biocontainers.pro/ga4gh/trs/v2/tools/centrifuger/versions) — exact image tag `1.1.0--hf426362_0`, 133.9 MB compressed

### Secondary (MEDIUM confidence)
- [Bioconda centrifuger recipe](https://bioconda.github.io/recipes/centrifuger/README.html) — version availability, quay.io/biocontainers tag pattern
- [WARP WDL cost optimization / disk autosizing patterns](https://broadinstitute.github.io/warp/docs/Best_practices/) — `ceil(size(...))` formula precedent
- [Terra scatter-gather documentation](https://support.terra.bio/hc/en-us/articles/360037128572-Scatter-gather-parallelism) — per-shard VM semantics confirming scatter anti-pattern

### Tertiary (LOW confidence)
- `duckdb` absent from py3-bio:0.1.3 — inferred from absence of any `import duckdb` in codebase; image not directly inspected

---
*Research completed: 2026-04-01*
*Ready for roadmap: yes*
