---
gsd_state_version: 1.0
milestone: v3.0
milestone_name: Centrifuger Taxonomic Classification WDL
current_phase: 04
status: verifying
stopped_at: Completed 04-01-PLAN.md
last_updated: "2026-04-02T13:58:21.413Z"
last_activity: 2026-04-02
progress:
  total_phases: 3
  completed_phases: 1
  total_plans: 1
  completed_plans: 1
  percent: 0
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-04-02)

**Core value:** Enables VirNucPro users to complete the full analysis chain within existing WDL infrastructure on Terra/DNAnexus
**Current focus:** Phase 04 — core-centrifuger-task

## Current Position

Phase: 04 (core-centrifuger-task) — EXECUTING
Plan: 1 of 1
Milestone: v3.0 Centrifuger Taxonomic Classification WDL
Current phase: 04
Status: Phase complete — ready for verification
Last activity: 2026-04-02

Progress: [░░░░░░░░░░] 0% (0/3 phases complete)

## Accumulated Context

### Previous Milestones

- **v1.0** (Shipped 2026-04-01): classify_virnucpro_contigs, classify_reads_by_contig WDL tasks + workflows
- **v1.1** (Shipped 2026-04-01): parse_kraken2_reads WDL task with DuckDB taxonomy annotation
- **v2.0** (Shipped 2026-04-01): summarize_kb_extract_reads WDL task with zstandard compression

### Decisions

Recent decisions affecting current work:

- [v1.0] Tasks in tasks_metagenomics.wdl — keep virnucpro classification co-located
- [v1.0] py3-bio + pip install duckdb — avoids custom image build
- [v1.0] python3<<CODE inline heredoc — consistent repo-wide Python-in-WDL pattern
- [v1.0] Dockstore entries without testParameterFiles — placeholder paths not for CI execution
- [v3.0] Use BioContainers centrifuger image (`quay.io/biocontainers/centrifuger:1.1.0--hf426362_0`) — only image with the centrifuger binary; samtools available via conda deps
- [v3.0] centrifuger_multi.wdl uses bash-loop (krakenuniq pattern), NOT WDL scatter — scatter would re-download 200+ GB index per shard, defeating amortization goal
- [v3.0] BAM inputs converted via samtools fastq inside task — Centrifuger has no native BAM support
- [v3.0] kreport generated inside same task via centrifuger-kreport — index must remain on disk; cannot split into separate task
- [v3.0] Disk autoscaling with ceil formula — same pattern as kraken2 task; do not hardcode disk size
- [v3.0] RAM default 240 GB — sized for NT-scale index; must not copy Kraken2's 90 GB default
- [v3.0] No testParameterFiles in .dockstore.yml — placeholder paths are not CI-runnable
- [Phase 04-core-centrifuger-task]: File centrifuger_db_tgz (tarball) used instead of Directory — WDL 1.0 does not support Directory type; matches kraken2/krakenuniq tarball pattern
- [Phase 04-core-centrifuger-task]: Single-end reads use centrifuger -u flag (not -r); -r is centrifuger-build reference input flag

### Pending Todos

- Verify BioContainers image tag `1.1.0--hf426362_0` is still pullable before writing docker declaration
- Confirm samtools availability in BioContainers centrifuger image before committing BAM conversion logic
- Run `grep -c "task centrifuger" tasks_metagenomics.wdl` before appending task (must return 0)
- Resolve RAM default (240 GB for NT-scale vs 50 GB for RefSeq-scope) with team before Phase 4 is finalized

### Blockers/Concerns

None.

## Session Continuity

Last session: 2026-04-02T13:58:21.411Z
Stopped at: Completed 04-01-PLAN.md
Resume file: None
Next action: `/gsd:plan-phase 4`
