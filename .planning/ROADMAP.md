# Roadmap

## Milestone 1: VirNucPro Contig & Read Classification WDL Pipelines

**Goal:** Two production-ready WDL tasks and standalone workflows that wrap the contig and read classification scripts, integrated into the existing `tasks_metagenomics.wdl` and validated with miniwdl.

---

## Phase 1: Implement WDL Tasks and Workflows

**Goal:** Add `classify_virnucpro_contigs` and `classify_reads_by_contig` tasks to `tasks_metagenomics.wdl`, create two standalone wrapper workflows, add test inputs, and update `.dockstore.yml`.

**Success criteria:**
- Both tasks present in `tasks_metagenomics.wdl` after the existing `classify_virnucpro` task
- Both workflows exist in `pipes/WDL/workflows/`
- `miniwdl check` passes on all four new/modified files
- Test input JSONs present in `test/input/WDL/miniwdl-local/`
- `.dockstore.yml` entries added for both workflows

**Requirements covered:** REQ-01 through REQ-19

**Plans:** 2/3 plans executed

Plans:
- [x] 01-01-PLAN.md -- classify_virnucpro_contigs task and standalone workflow
- [x] 01-02-PLAN.md -- classify_reads_by_contig task and standalone workflow
- [ ] 01-03-PLAN.md -- Test input JSONs, dockstore entries, and final validation

---

## Phase Order

```
Phase 1: Wave 1 (Plan 01) -> Wave 2 (Plan 02) -> Wave 3 (Plan 03)
```

Plans 01 and 02 are sequential because both modify `tasks_metagenomics.wdl`. Plan 03 depends on both.

---
*Last updated: 2026-03-31 after plan-phase*
