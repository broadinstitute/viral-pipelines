# Roadmap: VirNucPro Contig & Read Classification WDL Pipelines

## Milestones

- ✅ **v1.0 VirNucPro Classification WDL** — Phase 1 (shipped 2026-04-01)
- ✅ **v1.1 Kraken2 Read Taxonomy Annotation WDL Task** — Phase 2 (shipped 2026-04-01)

## Phases

<details>
<summary>✅ v1.0 VirNucPro Classification WDL (Phase 1) — SHIPPED 2026-04-01</summary>

- [x] Phase 1: Implement WDL Tasks and Workflows (4/4 plans) — completed 2026-04-01

See archive: `.planning/milestones/v1.0-ROADMAP.md`

</details>

### ✅ v1.1 Kraken2 Read Taxonomy Annotation WDL Task (SHIPPED 2026-04-01)

**Milestone Goal:** Promote `parse_kraken2_reads.py` into viral-pipelines as a production WDL task with standalone workflow, test inputs, and Dockstore registration — following the same inline-heredoc and runtime patterns established in v1.0.

- [x] **Phase 2: parse_kraken2_reads Task, Workflow, and Registration** - Implement and formally verify the Kraken2 taxonomy annotation task end-to-end

## Phase Details

### Phase 2: parse_kraken2_reads Task, Workflow, and Registration
**Goal**: The parse_kraken2_reads task and standalone workflow are production-ready, validated, and registered in Dockstore
**Depends on**: Phase 1
**Requirements**: REQ-01, REQ-02, REQ-03, REQ-04, REQ-05, REQ-06, REQ-07, REQ-08, REQ-09, REQ-10, REQ-11
**Success Criteria** (what must be TRUE):
  1. `parse_kraken2_reads` task exists in `tasks_metagenomics.wdl` with correct inputs, outputs, runtime, and `parameter_meta` matching neighbour task conventions — and passes `miniwdl check` with no errors
  2. Standalone workflow `pipes/WDL/workflows/parse_kraken2_reads.wdl` exists with a `meta` block (`description`, `author`, `email`, `allowNestedInputs: true`) and passes `miniwdl check`
  3. Test input JSON `test/input/WDL/miniwdl-local/test_inputs-parse_kraken2_reads-local.json` exists with workflow-level input keys and placeholder file paths
  4. `.dockstore.yml` contains an entry for `parse_kraken2_reads` workflow (no `testParameterFiles`)
  5. Task `parameter_meta` block includes `description`, `patterns`, and `category` for every input, matching the structure of `classify_virnucpro_contigs` and `classify_reads_by_contig`
**Plans:** 2/2 plans executed
Plans:
- [x] 02-01-PLAN.md — Implement parse_kraken2_reads WDL task with Python heredoc in tasks_metagenomics.wdl
- [x] 02-02-PLAN.md — Standalone workflow, test input JSON, and Dockstore registration

## Progress

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1. Implement WDL Tasks and Workflows | v1.0 | 4/4 | Complete | 2026-04-01 |
| 2. parse_kraken2_reads Task, Workflow, and Registration | v1.1 | 2/2 | Complete | 2026-04-01 |

---
*Last updated: 2026-04-01 — Phase 2 complete*
