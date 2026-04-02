---
phase: 06-test-input-jsons-and-dockstore-registration
verified: 2026-04-02T15:30:34Z
status: passed
score: 4/4 must-haves verified
re_verification: false
---

# Phase 06: Test Input JSONs and Dockstore Registration — Verification Report

**Phase Goal:** Create test input JSONs for centrifuger workflows and register them in Dockstore.
**Verified:** 2026-04-02T15:30:34Z
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| #  | Truth                                                                                                              | Status     | Evidence                                                                                                                      |
|----|-------------------------------------------------------------------------------------------------------------------|------------|-------------------------------------------------------------------------------------------------------------------------------|
| 1  | centrifuger_single test JSON exists with workflow-level input keys for reads_bam, centrifuger_db_tgz, and db_name | VERIFIED   | File at expected path; contains exactly 3 keys: `centrifuger_single.reads_bam`, `centrifuger_single.centrifuger_db_tgz`, `centrifuger_single.db_name` |
| 2  | centrifuger_multi test JSON exists with workflow-level input keys including array-typed reads_bams with 2 entries | VERIFIED   | File at expected path; `centrifuger_multi.reads_bams` is a 2-entry JSON array; exactly 3 keys total                          |
| 3  | Both centrifuger workflows are registered in .dockstore.yml with subclass WDL                                     | VERIFIED   | Lines 72-77 in `.dockstore.yml`; both entries have `subclass: WDL` and correct `primaryDescriptorPath`; YAML parses cleanly  |
| 4  | Neither Dockstore entry contains a testParameterFiles key                                                         | VERIFIED   | `testParameterFiles` does not appear in the centrifuger_multi or centrifuger_single sections; confirmed via grep and yaml parse |

**Score:** 4/4 truths verified

### Required Artifacts

| Artifact                                                                               | Expected                                             | Status   | Details                                                                                                  |
|----------------------------------------------------------------------------------------|------------------------------------------------------|----------|----------------------------------------------------------------------------------------------------------|
| `test/input/WDL/miniwdl-local/test_inputs-centrifuger_single-local.json`               | Placeholder test inputs for centrifuger_single       | VERIFIED | Exists, parses as valid JSON, 3 keys with `centrifuger_single.` prefix; no optional inputs               |
| `test/input/WDL/miniwdl-local/test_inputs-centrifuger_multi-local.json`                | Placeholder test inputs for centrifuger_multi        | VERIFIED | Exists, parses as valid JSON, 3 keys with `centrifuger_multi.` prefix; reads_bams is Array with 2 entries |
| `.dockstore.yml`                                                                        | Dockstore registration for both centrifuger workflows | VERIFIED | Contains both centrifuger entries; parses as valid YAML; entries confirmed programmatically              |

### Key Link Verification

| From                                               | To                                             | Via                                          | Status   | Details                                                                        |
|----------------------------------------------------|------------------------------------------------|----------------------------------------------|----------|--------------------------------------------------------------------------------|
| `test_inputs-centrifuger_single-local.json`        | `pipes/WDL/workflows/centrifuger_single.wdl`   | workflow-level input key prefix centrifuger_single.* | WIRED    | All 3 JSON keys use `centrifuger_single.` prefix; input names match WDL workflow `input {}` block (reads_bam, centrifuger_db_tgz, db_name) |
| `test_inputs-centrifuger_multi-local.json`         | `pipes/WDL/workflows/centrifuger_multi.wdl`    | workflow-level input key prefix centrifuger_multi.* | WIRED    | All 3 JSON keys use `centrifuger_multi.` prefix; input names match WDL workflow `input {}` block (reads_bams, centrifuger_db_tgz, db_name) |
| `.dockstore.yml`                                   | `pipes/WDL/workflows/centrifuger_single.wdl`   | primaryDescriptorPath                        | WIRED    | `primaryDescriptorPath: /pipes/WDL/workflows/centrifuger_single.wdl` present; target WDL file exists |
| `.dockstore.yml`                                   | `pipes/WDL/workflows/centrifuger_multi.wdl`    | primaryDescriptorPath                        | WIRED    | `primaryDescriptorPath: /pipes/WDL/workflows/centrifuger_multi.wdl` present; target WDL file exists   |

### Data-Flow Trace (Level 4)

Not applicable. These are static configuration/data files (JSON test inputs, YAML Dockstore config), not components that render or query dynamic data. No data-flow trace required.

### Behavioral Spot-Checks

| Behavior                                         | Command                                                                                                       | Result      | Status |
|--------------------------------------------------|---------------------------------------------------------------------------------------------------------------|-------------|--------|
| centrifuger_single JSON parses as valid JSON     | `python3 -c "import json; json.load(open(...))"` assertions on 3 keys                                        | ALL PASSED  | PASS   |
| centrifuger_multi JSON parses as valid JSON      | `python3 -c "import json; ..."` assertions on 3 keys; reads_bams is list of len 2                             | ALL PASSED  | PASS   |
| .dockstore.yml parses as valid YAML              | `python3 -c "import yaml; ..."` assertions on subclass, primaryDescriptorPath, no testParameterFiles          | ALL PASSED  | PASS   |
| Alphabetical ordering in .dockstore.yml          | Line position check: classify_krakenuniq(68) < centrifuger_multi(71) < centrifuger_single(74) < classify_multi(77) | CORRECT | PASS   |
| No optional inputs in either JSON                | Checked keys for machine_mem_gb, docker, cpu in both files                                                    | ABSENT      | PASS   |

### Requirements Coverage

| Requirement | Source Plan | Description                                                                                                    | Status    | Evidence                                                                          |
|-------------|-------------|----------------------------------------------------------------------------------------------------------------|-----------|-----------------------------------------------------------------------------------|
| CFGR-05     | 06-01-PLAN  | Test input JSON for centrifuger_single — placeholder paths, workflow-level input keys, v1.0-v2.0 conventions   | SATISFIED | File exists at correct path with correct keys; parses as valid JSON               |
| CFGR-06     | 06-01-PLAN  | Test input JSON for centrifuger_multi — placeholder Array inputs, workflow-level keys, v1.0-v2.0 conventions   | SATISFIED | File exists at correct path; reads_bams is Array[File]+ with 2 entries; valid JSON |
| CFGR-07     | 06-01-PLAN  | Dockstore registration entries in .dockstore.yml for both workflows; subclass WDL; no testParameterFiles       | SATISFIED | Both entries present in .dockstore.yml with subclass WDL; no testParameterFiles keys confirmed |

All 3 requirement IDs declared in the PLAN frontmatter are accounted for. No orphaned requirements for Phase 6 found in REQUIREMENTS.md.

### Anti-Patterns Found

None detected. All placeholder path values are intentional per project convention (D-04 from plan context: placeholder paths for local development test JSONs). The `placeholder` prefix in values is the established pattern, as verified by inspecting the reference file `test_inputs-parse_kraken2_reads-local.json` and confirmed in REQUIREMENTS.md ("Placeholder paths for index directory and FASTQ/BAM inputs").

### Human Verification Required

None. All items for this phase are static configuration artifacts (JSON and YAML files) that are fully verifiable programmatically.

### Gaps Summary

No gaps. All must-haves are fully satisfied:
- Both test input JSON files exist at their correct paths, parse as valid JSON, contain only required inputs with correct workflow-prefixed keys, and match the input contracts of their respective WDL workflows.
- Both centrifuger workflows are registered in `.dockstore.yml` with `subclass: WDL`, correct `primaryDescriptorPath` values, no `testParameterFiles`, and correct alphabetical placement (after `classify_krakenuniq`, before `classify_multi`).
- All 3 requirement IDs (CFGR-05, CFGR-06, CFGR-07) from the PLAN frontmatter are satisfied with evidence.
- Commits c1f858d8 and af018473 confirmed present in git history.

---

_Verified: 2026-04-02T15:30:34Z_
_Verifier: Claude (gsd-verifier)_
