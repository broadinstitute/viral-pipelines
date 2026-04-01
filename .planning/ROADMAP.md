# Roadmap: VirNucPro Contig & Read Classification WDL Pipelines

## Milestones

- ✅ **v1.0 VirNucPro Classification WDL** — Phase 1 (shipped 2026-04-01)
- ✅ **v1.1 Kraken2 Read Taxonomy Annotation WDL Task** — Phase 2 (shipped 2026-04-01)
- 🚧 **v2.0 KB Extract Read Summarization WDL** — Phase 3 (in progress)

## Phases

<details>
<summary>✅ v1.0 VirNucPro Classification WDL (Phase 1) — SHIPPED 2026-04-01</summary>

- [x] Phase 1: Implement WDL Tasks and Workflows (4/4 plans) — completed 2026-04-01

See archive: `.planning/milestones/v1.0-ROADMAP.md`

</details>

<details>
<summary>✅ v1.1 Kraken2 Read Taxonomy Annotation WDL Task (Phase 2) — SHIPPED 2026-04-01</summary>

- [x] Phase 2: parse_kraken2_reads Task, Workflow, and Registration (2/2 plans) — completed 2026-04-01

See archive: `.planning/milestones/v1.1-ROADMAP.md`

</details>

<details>
<summary>🚧 v2.0 KB Extract Read Summarization WDL (Phase 3) — PLANNED</summary>

- [ ] Phase 3: summarize_kb_extract_reads Task, Workflow, and Registration (0/3 plans) — planned 2026-04-01

**Phase Goal:** Convert summarize_extracted_reads.py into production WDL task and workflow

**Requirements Mapped:**
- SUMM-01: WDL task implementation in tasks_metagenomics.wdl
- SUMM-02: Standalone workflow wrapper
- SUMM-03: Test input JSON for miniwdl-local
- SUMM-04: Dockstore registration entry

**Plans Created:**
- [ ] `03-01-PLAN.md` — WDL task implementation (SUMM-01)
- [ ] `03-02-PLAN.md` — Standalone workflow wrapper (SUMM-02)
- [ ] `03-03-PLAN.md` — Test JSON + Dockstore registration (SUMM-03, SUMM-04)

**Success Criteria:**
1. Task passes `miniwdl check` validation
2. Workflow passes `miniwdl check` validation
3. Test JSON created with placeholder paths
4. Dockstore entry added to `.dockstore.yml`
5. All files committed to working branch

</details>

## Progress

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1. Implement WDL Tasks and Workflows | v1.0 | 4/4 | Complete | 2026-04-01 |
| 2. parse_kraken2_reads Task, Workflow, and Registration | v1.1 | 2/2 | Complete | 2026-04-01 |
| 3. summarize_kb_extract_reads Task, Workflow, and Registration | v2.0 | 3/3 | Complete   | 2026-04-01 |

---
*Last updated: 2026-04-01 — Phase 3 plans created (3 plans in 2 waves)*
