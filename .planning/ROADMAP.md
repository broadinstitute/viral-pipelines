# Roadmap: VirNucPro Contig & Read Classification WDL Pipelines

## Milestones

- ✅ **v1.0 VirNucPro Classification WDL** — Phase 1 (shipped 2026-04-01)
- ✅ **v1.1 Kraken2 Read Taxonomy Annotation WDL Task** — Phase 2 (shipped 2026-04-01)
- ✅ **v2.0 KB Extract Read Summarization WDL** — Phase 3 (shipped 2026-04-01)
- ✅ **v3.0 Centrifuger Taxonomic Classification WDL** — Phases 4-6 (shipped 2026-04-02)

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
<summary>✅ v2.0 KB Extract Read Summarization WDL (Phase 3) — SHIPPED 2026-04-01</summary>

- [x] Phase 3: summarize_kb_extract_reads Task, Workflow, and Registration (3/3 plans) — completed 2026-04-01

See archive: `.planning/milestones/v2.0-ROADMAP.md`

</details>

<details>
<summary>✅ v3.0 Centrifuger Taxonomic Classification WDL (Phases 4-6) — SHIPPED 2026-04-02</summary>

- [x] **Phase 4: Core centrifuger Task** - WDL task in tasks_metagenomics.wdl with BAM/FASTQ inputs, 4 outputs, input validation, and runtime sizing for 200+ GB index (completed 2026-04-02)
- [x] **Phase 5: centrifuger_single and centrifuger_multi Workflow Wrappers** - Standalone WDL workflows for single-sample and batch classification; multi uses bash-loop (not scatter) (completed 2026-04-02)
- [x] **Phase 6: Test Input JSONs and Dockstore Registration** - Placeholder test JSONs for both workflows and .dockstore.yml entries for Dockstore discoverability (completed 2026-04-02)

</details>

## Phase Details

### Phase 4: Core centrifuger Task
**Goal**: The centrifuger WDL task exists in tasks_metagenomics.wdl, handles all input types, emits all required outputs, and passes miniwdl check
**Depends on**: Nothing (additive to existing file)
**Requirements**: CFGR-01, CFGR-02
**Success Criteria** (what must be TRUE):
  1. `grep -c "task centrifuger" tasks_metagenomics.wdl` returns exactly 1 (no duplicate task body)
  2. The task accepts optional BAM input and converts it to FASTQ via `picard SamToFastq`; paired-end BAMs produce R1/R2 with `/1`/`/2` suffixes; single-end BAMs produce one FASTQ
  3. The task emits three named outputs: `classification_tsv`, `kreport`, `centrifuger_log` (`krona_html` deferred — KronaTools absent from Docker image)
  4. `miniwdl check pipes/WDL/tasks/tasks_metagenomics.wdl` exits 0 with no errors
  5. Runtime block specifies `cpu: 8`, `memory: "240 GB"`, and disk sized using the ceil autoscaling formula (not a hardcoded value)
**Plans**: 1 plan
Plans:
- [x] 04-01-PLAN.md — Append centrifuger task to tasks_metagenomics.wdl

### Phase 5: centrifuger_single and centrifuger_multi Workflow Wrappers
**Goal**: Both standalone workflow WDL files exist, import the task correctly, and pass miniwdl check; centrifuger_multi uses a bash-loop not scatter
**Depends on**: Phase 4
**Requirements**: CFGR-03, CFGR-04
**Success Criteria** (what must be TRUE):
  1. `miniwdl check pipes/WDL/workflows/centrifuger_single.wdl` exits 0
  2. `miniwdl check pipes/WDL/workflows/centrifuger_multi.wdl` exits 0
  3. `centrifuger_multi.wdl` contains no `scatter` block — classification loop lives in the task command block
  4. Both workflows include `meta { allowNestedInputs: true }` enabling Terra-compatible test JSON usage
  5. Both workflows call the task with alias `as run_centrifuger` (no call-name collision with workflow name)
**Plans**: 1 plan
Plans:
- [x] 05-01-PLAN.md — Modify centrifuger task to array inputs + create centrifuger_single and centrifuger_multi workflow wrappers

### Phase 6: Test Input JSONs and Dockstore Registration
**Goal**: Both workflows have placeholder test input JSONs and are registered in .dockstore.yml, making them discoverable on Dockstore
**Depends on**: Phase 5
**Requirements**: CFGR-05, CFGR-06, CFGR-07
**Success Criteria** (what must be TRUE):
  1. `test/input/WDL/miniwdl-local/test_inputs-centrifuger_single-local.json` exists with workflow-level input keys
  2. `test/input/WDL/miniwdl-local/test_inputs-centrifuger_multi-local.json` exists with workflow-level input keys and array-typed BAM/FASTQ placeholders
  3. `.dockstore.yml` contains entries for both `centrifuger_single.wdl` and `centrifuger_multi.wdl` with `subclass: WDL`
  4. Neither Dockstore entry contains a `testParameterFiles` key (placeholder paths must not be submitted for CI execution)
**Plans**: 1 plan
Plans:
- [x] 06-01-PLAN.md — Create test input JSONs and Dockstore registration entries

### Phase 7: join_read_classifications WDL Task, Workflow, and Registration
**Goal**: The join_read_classifications WDL task exists in tasks_metagenomics.wdl with the full 4-way DuckDB join logic embedded, a standalone workflow wrapper exists and passes miniwdl check, and the workflow is registered in .dockstore.yml
**Depends on**: Phase 6
**Requirements**: JRC-01, JRC-02, JRC-03, JRC-04
**Success Criteria** (what must be TRUE):
  1. `grep -c "task join_read_classifications" tasks_metagenomics.wdl` returns exactly 1
  2. The task accepts 4 optional File inputs (kallisto_summary, kraken2_reads, vnp_reads, genomad_virus_summary) and required String sample_id
  3. The task embeds the complete join_read_classifications.py logic as an inline python3 heredoc with DuckDB
  4. `miniwdl check pipes/WDL/tasks/tasks_metagenomics.wdl` exits 0
  5. `miniwdl check pipes/WDL/workflows/join_read_classifications.wdl` exits 0
  6. `.dockstore.yml` contains entry for `join_read_classifications.wdl` with `subclass: WDL` and no `testParameterFiles`
  7. `test/input/WDL/miniwdl-local/test_inputs-join_read_classifications-local.json` exists with workflow-level input keys
**Plans**: 2 plans
Plans:
- [ ] 07-01-PLAN.md — Append join_read_classifications task to tasks_metagenomics.wdl
- [ ] 07-02-PLAN.md — Create workflow wrapper, test input JSON, and Dockstore entry

## Progress

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1. Implement WDL Tasks and Workflows | v1.0 | 4/4 | Complete | 2026-04-01 |
| 2. parse_kraken2_reads Task, Workflow, and Registration | v1.1 | 2/2 | Complete | 2026-04-01 |
| 3. summarize_kb_extract_reads Task, Workflow, and Registration | v2.0 | 3/3 | Complete | 2026-04-01 |
| 4. Core centrifuger Task | v3.0 | 1/1 | Complete   | 2026-04-02 |
| 5. centrifuger_single and centrifuger_multi Workflow Wrappers | v3.0 | 1/1 | Complete   | 2026-04-02 |
| 6. Test Input JSONs and Dockstore Registration | v3.0 | 1/1 | Complete   | 2026-04-02 |
| 7. join_read_classifications Task, Workflow, and Registration | v3.1 | 0/2 | Planning | - |

---
*Last updated: 2026-04-02 — Phase 7 planned*
