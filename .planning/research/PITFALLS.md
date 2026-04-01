# Pitfalls Research

**Domain:** WDL pipeline wrapping a large-database taxonomic classifier (Centrifuger) on Terra/DNAnexus
**Researched:** 2026-04-01
**Confidence:** HIGH (codebase inspection + official docs) / MEDIUM (platform behavior, disk math)

---

## Critical Pitfalls

### Pitfall 1: Index Format Confusion — Centrifuger `.cfr` vs Centrifuge `.cf`

**What goes wrong:**
A pre-built Centrifuge index (files ending in `.1.cf`, `.2.cf`, etc.) is passed to `centrifuger`. The tool silently fails or produces nonsense output. The two tools use incompatible binary index formats despite having near-identical names.

**Why it happens:**
Centrifuger is the direct successor of Centrifuge and shares most CLI flags, so it is easy to assume the databases are interchangeable. Public databases are often distributed without explicit tool attribution, and older community resources (databases labeled "centrifuge-compatible") predate Centrifuger.

**How to avoid:**
- Declare `patterns: ["*.cfr"]` in the WDL `parameter_meta` for the index input to make the distinction explicit to Terra users.
- Add a runtime guard in the task command block that checks for the `.cfr` extension before calling `centrifuger`: e.g., `ls ${DB_DIR}/*.cfr || { echo "ERROR: index files must end in .cfr (Centrifuger format, not Centrifuge .cf)"; exit 1; }`.
- Document the distinction prominently in task `meta.description`.
- The Centrifuger README explicitly states the index prefix targets `.1.cfr`, `.2.cfr`, `.3.cfr`, `.4.cfr` files; treat this as the authoritative naming contract.

**Warning signs:**
- User provides a tarball with `.cf` extension files and no `.cfr` files.
- `centrifuger` exits with a cryptic error about index loading, not a classification failure.
- The task produces an empty output TSV with no per-read assignments.

**Phase to address:**
Requirements definition phase — specify the exact index format in the input contract. Executor phase — add the extension check guard in the task command block.

---

### Pitfall 2: Disk Size Under-Allocation for 200+ GB Index Tarballs

**What goes wrong:**
The WDL task allocates insufficient disk and the task fails mid-execution when extracting the index tarball. On GCP/Terra with local SSDs (`LOCAL` disk type), GCP allocates SSDs in 375 GB pairs. A disk formula that gives 500 GB will round up to 750 GB — but if the formula gives 749 GB it stays at 749 GB and the task fails. Mismatches between the tarball size and the decompressed index size cause the same failure.

**Why it happens:**
The existing Kraken2 task disk formula (`ceil((8 * size(reads_bam, "GB") + 3 * size(kraken2_db_tgz, "GB") + 50) / 750.0) * 750`) takes the localized tarball plus the decompressed files into account. For Centrifuger, the math is different: Centrifuger's compressed `.cfr` index for 140 Gbp RefSeq prokaryotes is ~41 GB on disk, but it decompresses to ~43 GB in RAM — not to disk. However a 200+ GB tarball (e.g., a larger-scope database) still needs at minimum: tarball (200 GB) + extracted files (200 GB) + FASTQ conversion buffer (3–8x BAM size) + outputs + overhead. This can exceed 1500 GB.

**How to avoid:**
- Model the disk formula explicitly as: `tarball_size + extracted_size + (BAM_to_FASTQ_expansion * bam_size) + overhead`. Expose `machine_disk_gb` as an overridable input with a computed default, using the same `ceil(.../ 750.0) * 750` rounding pattern already established in the `kraken2` task.
- Use `LOCAL` (not `HDD`) disk type for the Centrifuger task — consistent with every other heavy-classification task in this codebase (`kraken2`, `krakenuniq`, `blast`).
- Set a hard floor of 750 GB (one GCP SSD pair): `Int disk_size = if disk_size_auto < 750 then 750 else disk_size_auto`.
- Verify against the actual pre-built NCBI RefSeq full database before committing the formula.

**Warning signs:**
- Terra job fails with "No space left on device" or "disk quota exceeded" error in task log.
- Task fails during `tar -xf` or `tar -xzf` extraction, not during classification itself.
- `du -hs $DB_DIR` line (if present) shows file size close to allocated disk limit.

**Phase to address:**
Requirements definition — document database disk footprint. WDL task implementation phase — implement and validate disk formula with realistic sizes before shipping.

---

### Pitfall 3: Multi-Sample Workflow Localizes the 200 GB Database Once Per Scatter Shard

**What goes wrong:**
`centrifuger_multi.wdl` is designed to amortize the database load. But if the task is implemented the same way as `kraken2` (database as a `File` input, extracted inside the task `command` block), then a simple `scatter` over samples will localize and decompress the full database independently for every scatter shard. On Terra/GCP, each shard runs on a separate VM; the database is downloaded and extracted N times. For 20 samples, this is 20 × 200 GB = 4 TB of egress and 20 serial extraction operations, eliminating the amortization goal.

**Why it happens:**
WDL scatter semantics: every scatter shard is an independent task call on an independent VM. There is no "shared disk" between shards in Cromwell/Terra. The `classify_virnucpro_multi.wdl` pattern (scatter over `reads_input_set`) works cheaply because VirNucPro has no large shared database — each shard is independent. That pattern does not transfer to a 200+ GB database scenario without architectural changes.

**How to avoid:**
Choose one of two valid strategies and commit to it in the requirements phase — do not mix them:

**Option A (Terra scatter, per-shard localization):** Accept repeated localization. Each `centrifuger_single.wdl` call localizes the database independently. This is the simpler, Terra-idiomatic pattern. For cohort runs, amortize at the submission level by calling the single workflow across samples in parallel (Terra data table + workflow submission) rather than in a multi workflow. The database is still downloaded once per VM but that is the expected cost.

**Option B (True multi-sample amortization):** Implement `centrifuger_multi.wdl` as a single task (not a scatter) that takes `Array[File]+ reads_bams` and loops over samples inside the `command` block. The database loads once; classification loops. This is architecturally different from the existing `classify_virnucpro_multi.wdl` scatter pattern and requires the task to produce `Array[File]` outputs via file naming conventions, not WDL scatter gather. This approach has precedent in the codebase (e.g., `classify_kallisto_multi.wdl` handles multi-sample classification in a single command).

**Warning signs:**
- `centrifuger_multi.wdl` is written as `scatter(reads_bam in reads_bams) { call metagenomics.centrifuger { ... } }` — this is Option A behavior, not Option B.
- Terra job shows N identical database localization log entries for N samples.
- `centrifuger_multi.wdl` costs the same as N `centrifuger_single.wdl` submissions.

**Phase to address:**
Requirements definition phase — explicitly decide between Option A and Option B before writing WDL. This decision determines whether `centrifuger_multi.wdl` is a scatter-wrapper or a single batching task.

---

### Pitfall 4: BAM Input — Centrifuger Requires FASTQ; No Native BAM Support

**What goes wrong:**
The WDL task passes a BAM file directly to `centrifuger` via `-u` or `-1 -2`, expecting the tool to handle it as the `kraken2` WDL task does (where the viral-ngs Python wrapper handles BAM-to-FASTQ conversion internally via `metagenomics.kraken2`). Centrifuger has no BAM support — it accepts only FASTQ (gzipped or plain) and FASTA. Passing a BAM path results in a parse error or garbled per-read IDs.

**Why it happens:**
The existing `kraken2` task in this codebase delegates BAM handling to the `metagenomics` viral-ngs Python wrapper (`metagenomics kraken2 $DB_DIR/kraken2 ~{reads_bam} ...`), which handles BAM internally. Centrifuger is invoked as a direct binary (`centrifuger -x ...`), not through the viral-ngs wrapper, so the BAM abstraction is absent.

**How to avoid:**
Add an explicit BAM-to-FASTQ conversion step in the task command block before calling `centrifuger`. The established pattern in this codebase uses `samtools`:
- For paired-end BAM: `samtools sort -n ~{reads_bam} | samtools fastq -1 reads_1.fq.gz -2 reads_2.fq.gz -0 /dev/null -s /dev/null`
- Then call `centrifuger ... -1 reads_1.fq.gz -2 reads_2.fq.gz`
- For single-end: `samtools fastq ~{reads_bam} | gzip > reads.fq.gz && centrifuger ... -u reads.fq.gz`

The `samtools` binary is available in the `viral-ngs` Docker image. If using `py3-bio`, confirm samtools availability or switch to the viral-ngs image for this task (consistent with how `kraken2` uses `quay.io/broadinstitute/viral-ngs:3.0.10-classify`).

Expose `paired_end` as a WDL Boolean input to select the conversion branch.

**Warning signs:**
- `centrifuger` exits with `Error: Reads file is not in FASTQ format`.
- Per-read output TSV has all `taxID = 0` (unclassified) entries — indicates reads were not parsed.
- Task command block uses `~{reads_bam}` as a direct argument to `centrifuger` without a prior `samtools fastq` step.

**Phase to address:**
WDL task implementation phase — include BAM conversion in the command block specification. Research phase — confirm the Docker image has `samtools` available.

---

### Pitfall 5: Output Format Mismatch — Per-Read TSV vs Kraken-Style Report

**What goes wrong:**
The task outputs only the native Centrifuger per-read TSV (8 columns: readID, seqID, taxID, score, 2ndBestScore, hitLength, queryLength, numMatches). Downstream tools in the pipeline (`aggregate_metagenomics_reports`, `krona`) expect a Kraken-style summary report format. No Krona HTML is generated. The workflow integrates superficially with existing metagenomics outputs but cannot be aggregated or visualized using the existing tooling.

**Why it happens:**
Centrifuger's primary output is the per-read TSV (equivalent to Centrifuge's `-S` output). The Kraken-style summary report requires a second step: `centrifuger-kreport` (or the `centrifuger-quant` tool with `--output-format centrifuge`). This two-step pattern is not obvious from the tool name or basic documentation.

**How to avoid:**
Plan for two distinct output files from the task:
1. Per-read classification TSV (native Centrifuger format): `~{out_basename}.centrifuger.reads.tsv`
2. Kraken-style summary report: `~{out_basename}.centrifuger.report.txt` (generated via `centrifuger-kreport` from the per-read output)

If Krona visualization is needed, add a `centrifuger-quant` or `centrifuger-kreport` call after classification. The `centrifuger-kreport` binary is bundled in the same release as `centrifuger`.

Declare both output files in the WDL `output {}` block.

**Warning signs:**
- Task `output {}` block declares only one output file.
- Task command block calls `centrifuger` but does not call `centrifuger-kreport` or `centrifuger-quant`.
- `aggregate_metagenomics_reports` downstream receives Centrifuger per-read TSV instead of summary report.

**Phase to address:**
Requirements definition — decide what outputs the task must produce (per-read only, or per-read + summary). WDL task implementation — implement both outputs in the command block.

---

### Pitfall 6: Memory Ceiling Mismatch — RAM Required vs Machine Class

**What goes wrong:**
Centrifuger loads the entire FM-index into RAM at startup. For a RefSeq prokaryote + viral database (140 Gbp scope), this is ~43 GB. For a full NCBI NT-scope database (200+ GB scope), RAM requirements are proportionally higher and poorly documented. A task requesting 90 GB RAM (matching the `kraken2` default of `machine_mem_gb = 90`) may OOM if the Centrifuger database scope is larger than the Kraken2 database in use.

**Why it happens:**
Centrifuger memory is fully determined by the index size. Unlike Kraken2, which allows partial loading, Centrifuger loads the FM-index completely before processing any reads. The `--build-mem` parameter affects index construction only, not classification. There is no partial-index mode for classification.

**How to avoid:**
- Make `machine_mem_gb` an overridable workflow input (not hardcoded) with a default sized to the expected database. For the NCBI RefSeq full prokaryote database (~43 GB), 64 GB RAM is the safe minimum; 90 GB provides headroom for BAM conversion buffers. For a larger virus+prokaryote+human scope database, 128–256 GB may be required.
- Document the RAM requirement alongside the `dx_instance_type` selection in the task `runtime` block comment.
- Add a cgroup memory peak capture (`/sys/fs/cgroup/memory.peak`) to the task outputs, as the existing `kraken2` task does, to establish empirical baselines.

**Warning signs:**
- Task fails with "Killed" signal (OOM killer).
- `dmesg` or VM logs show out-of-memory events.
- The database tarball size is >60 GB compressed (decompressed index will require >90 GB RAM).

**Phase to address:**
Requirements definition — document database size and expected RAM. WDL task implementation — set `machine_mem_gb` default appropriately; expose as overridable input.

---

## Technical Debt Patterns

| Shortcut | Immediate Benefit | Long-term Cost | When Acceptable |
|----------|-------------------|----------------|-----------------|
| Hardcode `machine_mem_gb = 90` matching Kraken2 task | No decision needed now | Task silently OOMs on larger-scope databases | Never — expose as input with a documented default |
| Omit `centrifuger-kreport` step | Simpler command block | Output incompatible with existing aggregation/Krona tooling | Only if task is explicitly "per-read only" with no aggregation use case |
| Skip BAM extension guard; require user to pre-convert | Smaller task command block | User confusion; silent garbled reads | Only if the workflow explicitly mandates FASTQ-only input and enforces it via `patterns:` |
| Use `scatter` in multi-sample workflow without documenting repeated localization | Matches existing `classify_virnucpro_multi.wdl` pattern exactly | Users expect database amortization per PROJECT.md v3.0 goal; architecture misleads | Never — must document the behavior explicitly or implement true batching |
| Accept Centrifuger 1-file index prefix (no extension check) | No guard code needed | Centrifuge `.cf` indexes passed silently; classification produces wrong results | Never for production — always validate index extension |

---

## Integration Gotchas

| Integration | Common Mistake | Correct Approach |
|-------------|----------------|------------------|
| Terra scatter + large DB | Scatter over samples expecting shared DB localization | Each Terra scatter shard is an independent VM; DB downloaded per shard. Use single-task batching loop for true amortization. |
| DNAnexus `dx_instance_type` | Selecting `mem3_ssd1_v2_x32` (same as Kraken2) without verifying RAM for Centrifuger DB size | Verify that instance RAM exceeds the decompressed FM-index size; use `mem3_ssd1_v2_x32` (244 GB RAM) as minimum for large databases |
| Dockstore test parameter files | Adding `testParameterFiles` pointing to placeholder JSONs | Only add `testParameterFiles` when the JSON contains real, runnable test data. Placeholder JSONs cause Dockstore CI failures. |
| `miniwdl check` validation | Task passes `miniwdl check` but fails at runtime due to BAM input type | `miniwdl check` validates syntax only; it cannot validate that the tool accepts the input file format. Runtime behavior must be tested separately. |
| Centrifuger index as tarball | Tarballing the index the same way as Kraken2 (which uses `read_utils.extract_tarball` helper) | Centrifuger index files are `.cfr` numbered files; they can be tarballed but must preserve the prefix naming. The extraction path must match the `-x prefix` argument. |
| `viral-ngs` wrapper vs direct binary | Expecting `metagenomics centrifuger` CLI to exist in viral-ngs Docker image | Centrifuger is not (as of this writing) wrapped in the viral-ngs Python toolchain. The task must invoke the `centrifuger` binary directly, requiring the correct Docker image. |

---

## Performance Traps

| Trap | Symptoms | Prevention | When It Breaks |
|------|----------|------------|----------------|
| Per-shard DB localization in scatter multi | 20-sample run costs same as 20 single-sample runs; wall time determined by download speed | Implement true batching in `centrifuger_multi` command block; scatter only sample-level work, not DB loading | Every multi-sample run with N > 1 |
| Disk over-provisioning for SSD on GCP | `ceil(...)` formula gives 1501 GB → rounded to 2250 GB (6 × 375 GB) — correct but expensive | Round to nearest 750 GB (2-SSD pair) ceiling, not 375 GB (1-SSD) ceiling | When disk formula result falls between 750 GB multiples |
| BAM-to-FASTQ expansion underestimated | Disk full during FASTQ conversion before classification even starts | Use 8x BAM size for FASTQ expansion budget in disk formula (consistent with `kraken2` task comment) | BAM files with high complexity or large read counts |
| `preemptible: 2` on a 6-hour classification | GCP preemptible VMs are reclaimed after 24h max but jobs can be interrupted during the long DB load/extraction phase, wasting hours of compute | Set `preemptible: 0` for Centrifuger tasks exceeding 2 hours, or `preemptible: 1` with task-level restart (no partial output recovery) | Any run where DB extraction + classification exceeds 2 hours |

---

## "Looks Done But Isn't" Checklist

- [ ] **Index format guard:** Task declares `centrifuger -x $DB_DIR/centrifuger` — verify the tarball actually contains `.cfr` files and the extraction path matches the `-x` prefix argument.
- [ ] **BAM-to-FASTQ step present:** Task command block contains `samtools fastq` or equivalent before `centrifuger` invocation — not just `centrifuger ... ~{reads_bam}` directly.
- [ ] **Both output files declared:** Task `output {}` block has both per-read TSV and summary report — not just one.
- [ ] **Disk formula tested against actual database:** Disk `ceil(...)` formula evaluated with the real database file size; not copy-pasted from Kraken2 task without adjustment.
- [ ] **RAM default verified:** `machine_mem_gb` default is >= decompressed index size + 20 GB buffer — not blindly copied from Kraken2's `90`.
- [ ] **Docker image has `centrifuger` binary:** The selected Docker image actually contains the `centrifuger` binary — not just the old `centrifuge` binary. These are different packages.
- [ ] **Multi-sample architecture decision documented:** `centrifuger_multi.wdl` has a comment or meta block stating whether it uses per-shard localization (scatter) or single-task batching (loop), and why.
- [ ] **`miniwdl check` passes after BAM handling logic added:** Shell branching (`if [[ "~{reads_bam}" == *.bam ]]`) inside the heredoc is valid bash but verify miniwdl does not flag it.
- [ ] **`.dockstore.yml` entries without `testParameterFiles`:** Placeholder JSONs are not added as test parameter files.
- [ ] **`dx_instance_type` RAM verified:** DNAnexus instance type provides RAM >= `machine_mem_gb` value.

---

## Recovery Strategies

| Pitfall | Recovery Cost | Recovery Steps |
|---------|---------------|----------------|
| Index format confusion at runtime | LOW | Update test input JSON to point to `.cfr` tarball; re-run. No code change needed if guard is in place. |
| Disk OOM during extraction | LOW | Add `Int? machine_disk_gb` override to task; re-submit with larger disk. |
| RAM OOM during classification | LOW–MEDIUM | Re-submit with higher `machine_mem_gb` override. If no override input exists, requires a code edit and re-miniwdl-check. |
| Multi-sample workflow localizes DB N times | MEDIUM | Requires architectural rewrite of `centrifuger_multi.wdl` from scatter to single batching task. Cannot be patched at submission time. |
| BAM passed directly to centrifuger | LOW | Add `samtools fastq` step to command block; re-miniwdl-check; re-submit. |
| Wrong output format (no summary report) | LOW | Add `centrifuger-kreport` call to command block; add output declaration; re-miniwdl-check. |

---

## Pitfall-to-Phase Mapping

| Pitfall | Prevention Phase | Verification |
|---------|------------------|--------------|
| Index format confusion (`.cf` vs `.cfr`) | Requirements definition: document exact format in input contract | `parameter_meta.patterns` lists `*.cfr`; command block has extension check |
| Disk under-allocation for 200+ GB DB | Requirements: document DB disk footprint; WDL task implementation: validate formula | Disk formula evaluated with actual database size; formula matches or exceeds test calculation |
| Multi-sample localization anti-pattern | Requirements definition: decide scatter vs batching before writing WDL | `centrifuger_multi.wdl` meta block documents strategy; terra cost estimate matches expectation |
| BAM input (no native BAM support) | WDL task implementation: include BAM conversion step | Command block contains `samtools fastq` before `centrifuger`; task tested with BAM input |
| Output format mismatch (per-read only) | Requirements: enumerate required output files | Task `output {}` block has both TSV and summary report; downstream tools receive correct format |
| Memory ceiling for large index | Requirements: document RAM requirement; WDL task implementation: expose `machine_mem_gb` | `machine_mem_gb` default >= index size + 20 GB; value documented with rationale |
| Preemptible VM on long-running job | WDL task implementation | `preemptible` set to 0 or 1 for tasks expected to exceed 2 hours |
| Docker image lacks `centrifuger` binary | Requirements: confirm image availability; WDL task implementation | `docker pull` test before committing the image reference |

---

## Process Pitfalls (From Prior Milestones in This Repo)

These recur across milestones and apply directly to this work.

### P1: Executor Content Duplication When Appending to Large Files

**What goes wrong:**
Executor agent appends the new `centrifuger` task block to `tasks_metagenomics.wdl` without checking whether the content already exists. The task body is written twice. (Happened in v1.0; required a dedicated gap-closure plan.)

**How to avoid:**
Executor plan must include an explicit pre-append content check: `grep -c "task centrifuger" pipes/WDL/tasks/tasks_metagenomics.wdl` — must return 0 before writing. If it returns 1, skip the write.

**Phase to address:** Every executor wave that appends to an existing WDL task file.

---

### P2: Uncommitted Planning Files Cause Merge Conflicts on Worktree Branches

**What goes wrong:**
In-session planning changes (ROADMAP.md, STATE.md) are not committed before executor worktree branches are created. When the worktree branch is merged back, `STATE.md` conflicts require repeated stash/pop cycles. The `phase complete` CLI update to ROADMAP.md can be silently overwritten. (Happened in v1.1.)

**How to avoid:**
Before spawning each executor wave, commit all in-session planning changes: `git add -f .planning/ && git commit -m "docs: update planning state before executor wave"`. After merging each worktree, re-verify ROADMAP.md and STATE.md reflect the expected state.

**Phase to address:** Between every executor wave during implementation phases.

---

### P3: WDL Call Alias Required When Task Name Equals Workflow Name

**What goes wrong:**
Standalone workflow `centrifuger_single.wdl` (workflow name: `centrifuger_single`) calls `metagenomics.centrifuger`. If the workflow is named `centrifuger` (matching the task name), WDL 1.0 disallows the call without an alias and `miniwdl check` fails with a name collision error.

**How to avoid:**
Name standalone workflows `centrifuger_single` and `centrifuger_multi` (not `centrifuger`). If a plain `centrifuger.wdl` workflow is needed, use `call metagenomics.centrifuger as run_centrifuger`. Apply the alias pattern as a rule for all new standalone workflows (documented in PROJECT.md Key Decisions).

**Phase to address:** WDL workflow implementation — verify `miniwdl check` after writing each standalone workflow file.

---

## Sources

- [Centrifuger GitHub README (mourisl/centrifuger)](https://github.com/mourisl/centrifuger) — index format, CLI flags, output columns, memory requirements
- [Centrifuger Genome Biology 2024 paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03244-4) — 43 GB RAM for 140 Gbp RefSeq database
- [WARP WDL cost optimization](https://broadinstitute.github.io/warp/docs/Best_practices/GC_cost_optimization) — localization timing, reference file amortization
- [WARP disk autosizing patterns](https://broadinstitute.github.io/warp/docs/Best_practices/autosize) — `ceil(size(...))` pattern
- [Terra scatter-gather documentation](https://support.terra.bio/hc/en-us/articles/360037128572-Scatter-gather-parallelism) — per-shard VM semantics
- [Cromwell directory localization issue #5737](https://github.com/broadinstitute/cromwell/issues/5737) — large file repeated copy anti-pattern
- [nf-core/taxprofiler Centrifuge output documentation](https://nf-co.re/taxprofiler/1.1.8/docs/output) — expected output file names
- Direct codebase inspection:
  - `pipes/WDL/tasks/tasks_metagenomics.wdl` lines 206–341 — `kraken2` task (disk formula, decompress pattern, runtime attributes)
  - `pipes/WDL/workflows/classify_multi.wdl` — multi-sample scatter pattern
  - `pipes/WDL/workflows/classify_virnucpro_multi.wdl` — scatter-based multi-sample pattern
  - `.planning/RETROSPECTIVE.md` — v1.0 and v1.1 lessons (duplication guard, commit-before-wave, call alias)
  - `.planning/codebase/CONCERNS.md` — known issues (OOM, Docker sync, fragile areas)

---
*Pitfalls research for: Centrifuger WDL integration with 200+ GB database*
*Researched: 2026-04-01*
