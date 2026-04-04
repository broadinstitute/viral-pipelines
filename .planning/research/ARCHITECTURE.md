# Architecture Research: Centrifuger WDL Integration

**Domain:** WDL task/workflow integration — taxonomic classifier with large on-disk index
**Researched:** 2026-04-01
**Confidence:** HIGH — all structural findings from direct codebase inspection; docker image tag
requires registry confirmation (MEDIUM)

---

## Standard Architecture

### System Overview

```
                        centrifuger_multi.wdl
 ┌──────────────────────────────────────────────────────────────┐
 │  input: Array[File]+ reads_bams, File centrifuger_db_tgz     │
 │                                                              │
 │  call metagenomics.centrifuger as centrifuger_batch {        │
 │    input:                                                    │
 │      reads_bams         = reads_bams   (whole array)         │
 │      centrifuger_db_tgz = centrifuger_db_tgz                 │
 │  }                                                           │
 │     ^^^ single task call — ONE node, ONE db load             │
 │         task iterates all BAMs in bash loop                  │
 └──────────────────────────────────────────────────────────────┘

                        centrifuger_single.wdl
 ┌──────────────────────────────────────────────────────────────┐
 │  input: File reads_bam, File centrifuger_db_tgz              │
 │                                                              │
 │  call metagenomics.centrifuger {                             │
 │    input:                                                    │
 │      reads_bams         = [reads_bam]   (wrapped in array)   │
 │      centrifuger_db_tgz = centrifuger_db_tgz                 │
 │  }                                                           │
 └──────────────────────────────────────────────────────────────┘

                        tasks_metagenomics.wdl
 ┌──────────────────────────────────────────────────────────────┐
 │  task centrifuger {                                          │
 │    input: Array[File]+ reads_bams, File centrifuger_db_tgz   │
 │    command <<<                                               │
 │      # 1. decompress centrifuger_db_tgz once to $DB_DIR      │
 │      # 2. for bam in reads_bams: bam2fq, run centrifuger -x  │
 │      #    $DB_DIR/cfr_idx, run centrifuger-quant, gzip out   │
 │      # 3. emit Array[File] outputs via glob()                │
 │    >>>                                                       │
 │  }                                                           │
 └──────────────────────────────────────────────────────────────┘
```

### Why centrifuger_multi does NOT scatter

The critical design difference from `classify_multi.wdl` (which scatters individual `kraken2`
calls) is that `centrifuger_multi.wdl` does not scatter. Instead the task accepts
`Array[File]+` and iterates samples in a bash loop. The 200+ GB index is decompressed once
and held in RAM for all samples on a single node.

The authoritative precedent in this codebase is the `krakenuniq` task
(`tasks_metagenomics.wdl`, lines 83-84), which explicitly comments:

> "execute on all inputs and outputs serially, but with a single database load into ram"

Scattering `centrifuger` the way `classify_multi.wdl` scatters `kraken2` would force each
shard to stage and decompress the 200+ GB index independently. For a 20-sample batch that
becomes 20 separate node allocations and 4+ TB of data movement purely to load the database.

`centrifuger_single.wdl` passes a single-element array literal (`[reads_bam]`) to the same
task — no separate task definition is needed.

---

## Recommended Project Structure

```
pipes/WDL/
├── tasks/
│   └── tasks_metagenomics.wdl          MODIFY — ADD task centrifuger at end of file
└── workflows/
    ├── centrifuger_single.wdl           NEW FILE
    └── centrifuger_multi.wdl            NEW FILE

test/input/WDL/miniwdl-local/
    ├── centrifuger_single.json          NEW FILE
    └── centrifuger_multi.json           NEW FILE

.dockstore.yml                           MODIFY — ADD 2 entries
```

No other files require modification. All changes are additive except the two file appends
(tasks_metagenomics.wdl and .dockstore.yml).

---

## Architectural Patterns

### Pattern 1: Array-input task with bash loop (krakenuniq model)

**What:** The task declares `Array[File]+ reads_bams` and iterates samples in a bash loop
inside `command <<<`. The database is decompressed before the loop. Per-sample outputs use a
derived basename. The output block uses `glob()` to collect all per-sample files into arrays.

**When to use:** Any tool with a large in-memory database where database load cost dominates.
Krakenuniq (320 GB RAM, 750 GB disk, 32 CPU) uses this exact pattern. Centrifuger's
decompressed index is 200+ GB, so this is the correct model.

**Trade-offs:**
- Pro: Database loaded once; all samples run on one node; compute cost proportional to actual
  classification work, not to index I/O repeated N times.
- Con: One sample failure in the bash loop (if not trapped) can abort the whole batch. Mitigate
  with per-iteration error handling or accept all-or-nothing semantics for simplicity.

**Bash loop skeleton (adapted from krakenuniq, lines 77-90):**

```bash
for bam in ~{sep=' ' reads_bams}; do
  base=$(basename "$bam" .bam)
  samtools collate -u "$bam" \
    | samtools fastq \
        -1 "${base}.R1.fastq" \
        -2 "${base}.R2.fastq" \
        -0 "${base}.unpaired.fastq" \
        -s /dev/null -n
  centrifuger \
    -x "$DB_DIR/cfr_idx" \
    -1 "${base}.R1.fastq" -2 "${base}.R2.fastq" \
    -t ~{cpu} \
    > "${base}.centrifuger.tsv"
  centrifuger-quant \
    -x "$DB_DIR/cfr_idx" \
    -c "${base}.centrifuger.tsv" \
    > "${base}.centrifuger.report.txt"
  pigz "${base}.centrifuger.tsv"
done
```

**Output block (glob collection):**

```wdl
output {
  Array[File] centrifuger_hits_tsv_gz = glob("*.centrifuger.tsv.gz")
  Array[File] centrifuger_report_txt  = glob("*.centrifuger.report.txt")
  Int         max_ram_gb              = ceil(read_float("MEM_BYTES") / 1000000000)
  String      centrifuger_version     = read_string("VERSION")
}
```

The `max_ram_gb` cgroup memory capture idiom is present in krakenuniq (line 112) and kraken2
(line 318). Include it in the centrifuger task identically.

### Pattern 2: centrifuger_single.wdl thin wrapper

**What:** A minimal workflow that wraps a scalar `File reads_bam` input in a single-element
array and passes it to the `centrifuger` task. No `parameter_meta`, no scatter, no business
logic. Modeled on `classify_virnucpro_single.wdl` (28 lines total).

**When to use:** One sample per Terra submission, or when an outer workflow scatters and calls
`centrifuger_single` as its inner workflow.

**Exact structure to follow:**

```wdl
version 1.0
import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow centrifuger_single {
    meta {
        description: "Taxonomic classification of a single sample via Centrifuger."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }
    input {
        File  reads_bam
        File  centrifuger_db_tgz
    }
    call metagenomics.centrifuger {
        input:
            reads_bams         = [reads_bam],
            centrifuger_db_tgz = centrifuger_db_tgz
    }
    output {
        Array[File] centrifuger_hits_tsv_gz = centrifuger.centrifuger_hits_tsv_gz
        Array[File] centrifuger_report_txt  = centrifuger.centrifuger_report_txt
        Int         max_ram_gb              = centrifuger.max_ram_gb
        String      centrifuger_version     = centrifuger.centrifuger_version
    }
}
```

Note: the call name (`centrifuger`) will match the task name — no alias needed because the
workflow name (`centrifuger_single`) differs from the call name (`centrifuger`). This is
consistent with the existing WDL call-alias rule documented in PROJECT.md.

### Pattern 3: centrifuger_multi.wdl direct array pass-through

**What:** Passes `Array[File]+` directly to the task. No scatter. This is the same topology as
`classify_krakenuniq.wdl` (which calls the array-input `krakenuniq` task directly) and
contrasts with `classify_multi.wdl` (which scatters because kraken2 is per-sample).

**When to use:** Batch of 2–50 samples on one node. User accepts all-or-nothing execution and
wants to amortize the single large database load.

**Exact structure to follow:**

```wdl
version 1.0
import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow centrifuger_multi {
    meta {
        description: "Taxonomic classification of multiple samples via Centrifuger, amortizing the large index across all samples on a single node."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }
    input {
        Array[File]+ reads_bams
        File         centrifuger_db_tgz
    }
    parameter_meta {
        reads_bams: {
          description: "Set of BAM files to classify. May be unmapped or mapped, paired-end or single-end.",
          patterns: ["*.bam"]
        }
        centrifuger_db_tgz: {
          description: "Pre-built Centrifuger index tarball. Index files use the .cfr extension.",
          patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
        }
    }
    call metagenomics.centrifuger as centrifuger_batch {
        input:
            reads_bams         = reads_bams,
            centrifuger_db_tgz = centrifuger_db_tgz
    }
    output {
        Array[File] centrifuger_hits_tsv_gz = centrifuger_batch.centrifuger_hits_tsv_gz
        Array[File] centrifuger_report_txt  = centrifuger_batch.centrifuger_report_txt
        Int         max_ram_gb              = centrifuger_batch.max_ram_gb
        String      centrifuger_version     = centrifuger_batch.centrifuger_version
    }
}
```

The call alias `centrifuger_batch` is used to make the call name distinct from the workflow
name when needed; here the workflow is named `centrifuger_multi` so an alias is optional but
included for clarity and forward-compatibility.

---

## Data Flow

### Classification flow

```
reads_bams (Array[File]+) + centrifuger_db_tgz (File)
    |
    v
[task centrifuger — single container on one node]
    |-- 1. decompress centrifuger_db_tgz  -->  $DB_DIR/cfr_idx.*.cfr  (once, ~200 GB)
    |-- 2. capture centrifuger version    -->  VERSION file
    |-- 3. for each BAM:
    |       samtools bam2fq  -->  R1/R2 FASTQ
    |       centrifuger -x $DB_DIR/cfr_idx -1 R1 -2 R2  -->  basename.centrifuger.tsv
    |       centrifuger-quant -x $DB_DIR/cfr_idx -c tsv  -->  basename.centrifuger.report.txt
    |       pigz basename.centrifuger.tsv
    |-- 4. capture peak RAM               -->  MEM_BYTES
    v
centrifuger_hits_tsv_gz: Array[File]   ("*.centrifuger.tsv.gz")
centrifuger_report_txt:  Array[File]   ("*.centrifuger.report.txt")
max_ram_gb:              Int
centrifuger_version:     String
```

### Key data flows

1. **single.wdl path:** Scalar `reads_bam: File` is wrapped: `reads_bams = [reads_bam]`. Task
   outputs are `Array[File]` of length 1. The workflow output block re-exports them identically
   to the multi path — callers need no special-casing.

2. **multi.wdl path:** `reads_bams: Array[File]+` passed directly to the task. No wrapping.
   Task processes all samples serially under the same loaded index.

---

## Integration Points

### tasks_metagenomics.wdl insertion point

The new `task centrifuger` appends to the end of `tasks_metagenomics.wdl`. The file currently
ends at line 2353 (end of `summarize_kb_extract_reads` task). No existing content is modified.

| File | Action | Location |
|------|--------|----------|
| `pipes/WDL/tasks/tasks_metagenomics.wdl` | ADD `task centrifuger` | After line 2353 (end of file) |
| `pipes/WDL/workflows/centrifuger_single.wdl` | CREATE | `pipes/WDL/workflows/` |
| `pipes/WDL/workflows/centrifuger_multi.wdl` | CREATE | `pipes/WDL/workflows/` |
| `test/input/WDL/miniwdl-local/centrifuger_single.json` | CREATE | test input dir |
| `test/input/WDL/miniwdl-local/centrifuger_multi.json` | CREATE | test input dir |
| `.dockstore.yml` | ADD 2 entries | Under WDL workflows section |

### .dockstore.yml entries to add

The pattern used for `classify_virnucpro_single` and `classify_virnucpro_multi` (.dockstore.yml
lines 322-327) is the canonical form. Apply it identically:

```yaml
  - name: centrifuger_single
    subclass: WDL
    primaryDescriptorPath: /pipes/WDL/workflows/centrifuger_single.wdl
  - name: centrifuger_multi
    subclass: WDL
    primaryDescriptorPath: /pipes/WDL/workflows/centrifuger_multi.wdl
```

No `testParameterFiles` entries — consistent with the v1.0–v2.0 precedent of omitting test
parameter paths until real test data exists.

### Docker image

No `quay.io/broadinstitute/viral-classify` tag currently ships centrifuger binary. The correct
image is from Bioconda's BioContainers automated builds:

```
quay.io/biocontainers/centrifuger:<tag>
```

Latest upstream release: v1.1.0 (February 2026). The exact `<tag>` string must be confirmed
against the live registry at `https://quay.io/repository/biocontainers/centrifuger?tab=tags`
before writing the task. This is a blocking prerequisite for Phase 1.

**Confidence: MEDIUM** — The Bioconda package for centrifuger exists and BioContainers
auto-publishes it to quay.io; exact versioned tag requires registry lookup.

---

## Runtime Specification

### 200+ GB index: resource sizing rationale

Published evidence from the centrifuger paper and GitHub:
- RefSeq prokaryotic index (~140G nt): ~43 GB RAM at classification time
- Full NT-scale build requires `--build-mem 240G` on a 300 GB server
- A 200 GB on-disk tarball decompresses to a comparable working-set footprint

The closest analog in this codebase is the `krakenuniq` task
(`tasks_metagenomics.wdl`, lines 126-134): 320 GB RAM, 32 CPU, 750 GB disk.

Centrifuger can be run with fewer threads and the index is more compact (FM-index vs full k-mer
hash), so the RAM floor can be set lower while remaining user-overridable:

| Runtime parameter | Recommended value | Rationale |
|-------------------|-------------------|-----------|
| `machine_mem_gb` | `240` (default) | Covers NT-scale index; user-overridable with `Int machine_mem_gb = 240` in input block |
| `cpu` | `16` (default) | Centrifuger BWT traversal is CPU-bound; matches kraken2 task (line 334); user-overridable |
| `disk_size` floor | `750` | Minimum for any 200 GB index (1x compressed + 1x decompressed + overhead) |
| `disk_size` autoscale | formula below | Scales with actual db size and BAM count |
| `dx_instance_type` | `mem3_ssd1_v2_x16` | Closest DNAnexus instance to 240 GB RAM / 16 CPU |
| `preemptible` | `0` | Long db-load + multi-sample serial run; preemption wastes hours of work |
| disk type | `LOCAL` | Heavy I/O, large database; same as krakenuniq and kraken2 tasks |

**Disk autoscaling formula (WDL, matches kraken2 pattern at lines 248-250):**

```wdl
Int disk_size_raw = ceil((2 * size(centrifuger_db_tgz, "GB") + 10 * length(reads_bams) + 50) / 750.0) * 750
Int disk_size     = if disk_size_raw < 750 then 750 else disk_size_raw
```

The `2 *` factor covers compressed + decompressed index. The `10 * length(reads_bams)` factor
covers BAM-to-FASTQ expansion (~8x) plus classification output (~1x), rounded up. `50` is
fixed overhead for temp files and OS. GCP local SSDs are allocated in 375 GB increments (so
750 GB multiples); the formula rounds up to the next 750 GB boundary, matching the established
kraken2 pattern exactly.

**Full runtime block:**

```wdl
runtime {
  docker:           docker
  memory:           "~{machine_mem_gb} GB"
  cpu:              cpu
  disks:            "local-disk ~{disk_size} LOCAL"
  disk:             "~{disk_size} GB" # TES
  dx_instance_type: "mem3_ssd1_v2_x16"
  preemptible:      0
}
```

---

## Build Order

Dependencies are strictly ordered:

```
Step 1: task centrifuger in tasks_metagenomics.wdl
         → passes: miniwdl check pipes/WDL/tasks/tasks_metagenomics.wdl
              ↓
Step 2: centrifuger_single.wdl  (imports tasks_metagenomics.wdl)
        centrifuger_multi.wdl   (imports tasks_metagenomics.wdl)
         → passes: miniwdl check pipes/WDL/workflows/centrifuger_single.wdl
                   miniwdl check pipes/WDL/workflows/centrifuger_multi.wdl
              ↓
Step 3: test JSONs for both workflows
         → passes: miniwdl check --path ... centrifuger_single.wdl
              ↓
Step 4: .dockstore.yml entries for both workflows
```

Do not write the workflow files before the task passes `miniwdl check`. A task-level syntax
error will cause both workflow imports to fail, obscuring the root error.

---

## Anti-Patterns

### Anti-Pattern 1: Scatter over centrifuger calls in centrifuger_multi

**What people do:** Mirror the `classify_multi.wdl` scatter pattern, creating one call per BAM,
each call receiving the same `centrifuger_db_tgz` file.

**Why it's wrong:** WDL/Cromwell stages file inputs per scatter shard. A 20-sample batch
downloads and decompresses the 200+ GB index 20 times, across 20 separate node allocations.
This is the anti-pattern `classify_multi.wdl` explicitly avoids for krakenuniq by delegating
to the krakenuniq task (which handles the loop internally).

**Do this instead:** Pass the whole `Array[File]+` to a single task call. Bash loop inside
the task command iterates samples against the already-loaded index. This is documented
explicitly in the `krakenuniq` command block comment at line 83-84.

### Anti-Pattern 2: Separate task definition for centrifuger_single

**What people do:** Define a second task `centrifuger_single_task` with scalar `File reads_bam`
input to serve the single-sample workflow, duplicating the command block.

**Why it's wrong:** Identical command block maintained in two places. Any bug fix or flag change
must be applied twice.

**Do this instead:** One task, `Array[File]+` input. `centrifuger_single.wdl` wraps the scalar
with `reads_bams = [reads_bam]`. This pattern is already established in this repo — see how
`classify_virnucpro_single.wdl` and `classify_virnucpro_multi.wdl` both call the same
`classify_virnucpro` task.

### Anti-Pattern 3: Hardcoded disk size without autoscaling

**What people do:** Set `Int disk_size = 750` as a constant.

**Why it's wrong:** A full NT-scale index tarball is 200+ GB compressed; decompressed working
set doubles that. Larger batches add more BAM expansion overhead. A fixed 750 GB will be
insufficient for large batches or future larger indexes.

**Do this instead:** Use the `ceil((2 * size(db, "GB") + 10 * length(bams) + 50) / 750.0) * 750`
autoscaling formula with a 750 GB floor, matching the kraken2 pattern exactly.

### Anti-Pattern 4: Omitting `allowNestedInputs: true`

**What people do:** Forget the `allowNestedInputs` metadata key.

**Why it's wrong:** Without it, `miniwdl run` with a flat JSON (keys like
`centrifuger_single.reads_bam`) silently ignores inputs rather than routing them to the task.
All four workflows added in v1.0–v2.0 of this project include this key.

**Do this instead:** Add `allowNestedInputs: true` to the `meta {}` block of both
`centrifuger_single.wdl` and `centrifuger_multi.wdl`.

### Anti-Pattern 5: Per-sample `centrifuger-quant` separate from classification

**What people do:** Run `centrifuger-quant` in a separate downstream task or as a separate
WDL call after each classification.

**Why it's wrong:** `centrifuger-quant` requires the index (`-x`) and the per-read hits TSV
(`-c`). Running it inside the same task reuses the already-loaded index from `$DB_DIR`
without additional staging.

**Do this instead:** Call `centrifuger-quant` immediately after `centrifuger` for each sample
inside the bash loop, before moving to the next sample. Include both outputs in `glob()`.

---

## Component Boundaries (v3.0 context)

| Component | Responsibility | Input | Output |
|-----------|---------------|-------|--------|
| `centrifuger` task (new) | Taxonomic classification via FM-index; loop over samples, single db load | `Array[File]+ reads_bams`, `File centrifuger_db_tgz` | per-sample hits TSV.gz, report TXT |
| `centrifuger_single.wdl` (new) | Wrap one sample for Terra single-sample workflow | `File reads_bam`, `File centrifuger_db_tgz` | pass-through of task outputs |
| `centrifuger_multi.wdl` (new) | Batch wrapper, amortize db load | `Array[File]+ reads_bams`, `File centrifuger_db_tgz` | pass-through of task outputs |
| `tasks_metagenomics.wdl` | Task library (existing) | imports by workflows | n/a |
| `.dockstore.yml` | Workflow registry (existing) | modified to add 2 entries | n/a |

---

## Sources

- `pipes/WDL/tasks/tasks_metagenomics.wdl` lines 1-135: `krakenuniq` task — array-input,
  bash loop, single db load, glob output, max_ram_gb capture pattern
- `pipes/WDL/tasks/tasks_metagenomics.wdl` lines 206-341: `kraken2` task — disk autoscaling
  formula, runtime block structure, version capture pattern
- `pipes/WDL/workflows/classify_virnucpro_single.wdl` and `classify_virnucpro_multi.wdl`:
  single/multi thin-wrapper pattern; `allowNestedInputs: true`; single-element array wrap
- `pipes/WDL/workflows/classify_krakenuniq.wdl`: direct call to array-input task, no scatter
- `.dockstore.yml` lines 322-327: single/multi registration pattern without testParameterFiles
- `.planning/PROJECT.md` v3.0 milestone: requirements, out-of-scope items, key decisions
- [Centrifuger GitHub](https://github.com/mourisl/centrifuger): CLI flags (`-x`, `-1`, `-2`,
  `-u`, `-t`), per-read TSV output (8-column), centrifuger-quant for reports, `.cfr` index
  format, v1.1.0 latest release (February 2026)
- [Centrifuger Genome Biology 2024](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03244-4):
  memory requirements (43 GB for 140G nt; 240 GB build for NT scale)
- [Bioconda centrifuger recipe](https://bioconda.github.io/recipes/centrifuger/README.html):
  `quay.io/biocontainers/centrifuger` image confirmed; exact tag requires registry lookup

---
*Architecture research for: Centrifuger WDL task and single/multi workflows (v3.0 milestone)*
*Researched: 2026-04-01*
