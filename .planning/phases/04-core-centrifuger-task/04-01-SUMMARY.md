---
phase: 04-core-centrifuger-task
plan: 01
subsystem: wdl-tasks
tags: [centrifuger, wdl, metagenomics, classification, picard, bam-to-fastq]
dependency_graph:
  requires: []
  provides: [centrifuger WDL task in tasks_metagenomics.wdl]
  affects: [pipes/WDL/workflows/ (Phase 5 callers)]
tech_stack:
  added: [centrifuger 1.1.0 (ghcr.io/broadinstitute/docker-centrifuger:1.0.0), picard SamToFastq]
  patterns: [tarball index extraction, optional-File BAM/FASTQ dispatch, ceil disk autoscaling]
key_files:
  modified: [pipes/WDL/tasks/tasks_metagenomics.wdl]
decisions:
  - "Use File centrifuger_db_tgz (tarball) instead of Directory centrifuger_db — WDL 1.0 does not support Directory type; tarball is consistent with kraken2 and krakenuniq index patterns"
  - "Single-end reads use centrifuger -u flag (not -r); -r is a centrifuger-build flag for reference sequences"
  - "kreport generation co-located in same task — centrifuger-kreport requires the index to remain on disk"
  - "Disk formula includes 3x tarball size: ceil((8 * size(reads_fastq1) + 3 * size(centrifuger_db_tgz) + 400) / 400.0) * 400"
metrics:
  duration: ~6 minutes
  completed: 2026-04-02T13:57:03Z
  tasks_completed: 2
  files_modified: 1
---

# Phase 04 Plan 01: centrifuger WDL Task Summary

**One-liner:** centrifuger WDL task appended to tasks_metagenomics.wdl with picard SamToFastq BAM conversion, paired/unpaired FASTQ support, tarball index extraction, and kreport generation — passes miniwdl check.

## What Was Built

Appended the `centrifuger` task to `pipes/WDL/tasks/tasks_metagenomics.wdl` (lines 944–1111). The task:

- Accepts a compressed Centrifuger index tarball (`File centrifuger_db_tgz`) and index prefix (`String db_name`)
- Accepts optional paired-end FASTQ (`reads_fastq1` / `reads_fastq2`), single-end FASTQ (`reads_fastq_unpaired`), or BAM (`reads_bam`)
- Converts BAM inputs to FASTQ at runtime via `picard SamToFastq` with exact CFGR-02 flags (READ1_SUFFIX=/1, READ2_SUFFIX=/2)
- Detects paired vs single-end from non-empty R2.fq after BAM conversion
- Runs `centrifuger` with `-1/-2` (paired) or `-u` (single-end/unpaired)
- Runs `centrifuger-kreport` in the same task (index must remain on disk)
- Emits three outputs: `classification_tsv`, `kreport`, `centrifuger_log`
- Runtime: 240 GB memory, cpu 8, ceil disk autoscaling with 400 GB floor, Docker `ghcr.io/broadinstitute/docker-centrifuger:1.0.0`

## Tasks Completed

| Task | Name | Commit | Files |
|------|------|--------|-------|
| 1 | Pre-flight validation and append centrifuger task | 79de5db4 | pipes/WDL/tasks/tasks_metagenomics.wdl |
| 1 (fix) | Replace Directory with File tarball (auto-fix) | f83e30c3 | pipes/WDL/tasks/tasks_metagenomics.wdl |
| 2 | Validate with miniwdl check — all 6 checks pass | (no new commit needed) | — |

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] WDL 1.0 does not support the Directory type**

- **Found during:** Task 2 (miniwdl check)
- **Issue:** Plan specified `Directory centrifuger_db` as the index input type, but `Directory` is a WDL 2.0 feature. miniwdl 1.13.1 running against `version 1.0` files rejects it with `Unknown type Directory`. The plan's own success criterion #4 (`miniwdl check` must exit 0) conflicts with using `Directory` in a `version 1.0` file.
- **Fix:** Replaced `Directory centrifuger_db` with `File centrifuger_db_tgz` (tarball archive), consistent with the existing `kraken2_db_tgz` and `krakenuniq_db_tar_lz4` pattern. Added runtime tarball extraction logic supporting `.tar.gz`, `.tar.lz4`, `.tar.zst`, and `.tar.bz2`. Updated disk formula to include `3 * size(centrifuger_db_tgz, "GB")` in addition to the 400 GB floor. `INDEX_PREFIX` bash variable replaces the WDL path expression `~{centrifuger_db}/~{db_name}`.
- **Files modified:** `pipes/WDL/tasks/tasks_metagenomics.wdl`
- **Commit:** f83e30c3
- **Impact on REQUIREMENTS.md:** CFGR-01 specifies `Directory centrifuger_db` — this must be updated to reflect `File centrifuger_db_tgz`. Phase 5 callers must use the tarball input pattern.

## Decisions Made

| Decision | Rationale |
|----------|-----------|
| `File centrifuger_db_tgz` over `Directory centrifuger_db` | WDL 1.0 compatibility; consistent with existing kraken2/krakenuniq tarball pattern |
| Single-end flag `-u` (not `-r`) | Centrifuger docs use `-u`; `-r` is centrifuger-build's reference input flag |
| kreport in same task as classification | centrifuger-kreport requires the index to remain on disk; architecturally cannot split |
| 400 GB disk floor | Covers NT-scale centrifuger index; plus 3x tarball size for extraction overhead |
| `preemptible: 2` | Follows newer task pattern in this file (not the older krakenuniq `preemptible: 0`) |

## Success Criteria Verification

| Criterion | Result |
|-----------|--------|
| Exactly 1 `task centrifuger` | PASS — grep -c returns 1 |
| picard SamToFastq with READ1_SUFFIX=/1, READ2_SUFFIX=/2 | PASS |
| Three outputs: classification_tsv, kreport, centrifuger_log | PASS |
| No krona_html output (deferred) | PASS — grep -c returns 0 |
| miniwdl check exits 0 | PASS |
| memory 240 GB, cpu 8, ceil disk formula | PASS |
| Single-end uses `-u` flag (not `-r`) | PASS |

## Self-Check: PASSED

- File exists: pipes/WDL/tasks/tasks_metagenomics.wdl — confirmed (1111 lines)
- Commit 79de5db4 exists — confirmed
- Commit f83e30c3 exists — confirmed
- miniwdl check exits 0 — confirmed
- task centrifuger at line 944 — confirmed
