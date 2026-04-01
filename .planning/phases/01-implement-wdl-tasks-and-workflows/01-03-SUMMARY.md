---
phase: 01-implement-wdl-tasks-and-workflows
plan: 03
subsystem: testing-and-registry
tags: [wdl, dockstore, test-inputs, validation, miniwdl]
dependency_graph:
  requires: [01-01, 01-02]
  provides: [test-input-jsons, dockstore-entries]
  affects: [ci-testing, terra-discovery]
tech_stack:
  added: []
  patterns: [miniwdl-local test JSON format with workflow-level input keys]
key_files:
  created:
    - test/input/WDL/miniwdl-local/test_inputs-classify_virnucpro_contigs-local.json
    - test/input/WDL/miniwdl-local/test_inputs-classify_reads_by_contig-local.json
  modified:
    - .dockstore.yml
decisions:
  - JSON keys use workflow-level input format (workflow.input) not call-alias format, matching miniwdl allowNestedInputs pattern
  - Dockstore entries added without testParameterFiles since placeholder paths are not for CI execution
  - New entries placed after classify_single to keep all classify workflows grouped in dockstore.yml
metrics:
  duration: ~2m
  completed: "2026-04-01T14:21:26Z"
  tasks_completed: 2
  tasks_total: 2
  files_created: 2
  files_modified: 1
---

# Phase 01 Plan 03: Test Inputs and Dockstore Integration Summary

Test input JSON files for both new classification workflows created, dockstore.yml updated with entries for both workflows, and all three WDL files confirmed passing miniwdl check validation.

## Tasks Completed

### Task 1: Create test input JSONs for both workflows

Created two test input JSON files in `test/input/WDL/miniwdl-local/`:

- `test_inputs-classify_virnucpro_contigs-local.json` — workflow-level input for `virnucpro_scores_tsv` with placeholder path
- `test_inputs-classify_reads_by_contig-local.json` — workflow-level inputs for `paf_file` and `contig_classifications` with placeholder paths

Both files created as direct files (not symlinks), matching the pattern of `test_inputs-simulate_illumina_reads-local.json` in the same directory.

**Commit:** `4893eb93`

### Task 2: Add dockstore entries and run final validation

Added two entries to `.dockstore.yml` immediately after `classify_single` to keep all classify workflows grouped together:

```yaml
  - name: classify_virnucpro_contigs
    subclass: WDL
    primaryDescriptorPath: /pipes/WDL/workflows/classify_virnucpro_contigs.wdl
  - name: classify_reads_by_contig
    subclass: WDL
    primaryDescriptorPath: /pipes/WDL/workflows/classify_reads_by_contig.wdl
```

Verified `py3-bio` already present in `requirements-modules.txt` (version 0.1.5) — no modification required.

All three WDL files pass `miniwdl check` with exit code 0:
- `pipes/WDL/tasks/tasks_metagenomics.wdl` — exit 0
- `pipes/WDL/workflows/classify_virnucpro_contigs.wdl` — exit 0
- `pipes/WDL/workflows/classify_reads_by_contig.wdl` — exit 0

**Commit:** `d3deba25`

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Merged prior agent work before execution**
- **Found during:** Pre-execution setup
- **Issue:** The worktree branch `worktree-agent-a02a9972` did not have the WDL files created by plans 01-01 and 01-02, which are on branch `worktree-agent-ade8ca6e`
- **Fix:** Ran `git merge worktree-agent-ade8ca6e` to bring in all prerequisite work before executing plan 03
- **Commit:** Merge commit created automatically

**2. [Plan Correction - Naming] JSON key format uses workflow-level inputs, not call-alias format**
- **Found during:** Task 1 — reading actual WDL workflows
- **Issue:** The plan acceptance criteria specifies key format `classify_virnucpro_contigs.classify_virnucpro_contigs.virnucpro_scores_tsv` (using call alias matching task name), but the actual WDL uses aliases `classify_contigs` and `classify_reads` (set in 01-01/01-02 due to WDL naming constraint where call name cannot equal workflow name)
- **Fix:** Used correct workflow-level input format `{workflow_name}.{input_name}` which matches the `allowNestedInputs: true` pattern used in other miniwdl test inputs (e.g., `simulate_illumina_reads.accession_coverage_string`). This is what miniwdl actually requires.
- **Files modified:** The test JSON keys use `classify_virnucpro_contigs.virnucpro_scores_tsv` and `classify_reads_by_contig.paf_file` formats

**3. [Plan Adaptation] Dockstore entries inserted after classify_single, not classify_virnucpro_multi**
- **Found during:** Task 2 — reading .dockstore.yml
- **Issue:** The plan says to insert after `classify_virnucpro_multi` entry, but this worktree's .dockstore.yml has no existing virnucpro entries (they are in another agent's worktree that hasn't been merged yet)
- **Fix:** Inserted after `classify_single` to keep all classify workflows grouped together, which is the intended organizational pattern

## Verification Results

| Check | Result |
|-------|--------|
| test_inputs-classify_virnucpro_contigs-local.json exists | PASS |
| test_inputs-classify_reads_by_contig-local.json exists | PASS |
| Both JSON files valid JSON | PASS |
| .dockstore.yml contains classify_virnucpro_contigs | PASS |
| .dockstore.yml contains classify_reads_by_contig | PASS |
| miniwdl check tasks_metagenomics.wdl exits 0 | PASS |
| miniwdl check classify_virnucpro_contigs.wdl exits 0 | PASS |
| miniwdl check classify_reads_by_contig.wdl exits 0 | PASS |
| requirements-modules.txt unchanged (py3-bio present) | PASS |

## Known Stubs

The test input JSON files use placeholder file paths (`test/input/placeholder_*.tsv`, `test/input/placeholder_*.paf.gz`). These placeholder paths are intentional — no real test data exists yet. The JSON structure is correct for CI validation; actual test execution is deferred to a future plan when real test data is available.

## Self-Check: PASSED

- test/input/WDL/miniwdl-local/test_inputs-classify_virnucpro_contigs-local.json: FOUND
- test/input/WDL/miniwdl-local/test_inputs-classify_reads_by_contig-local.json: FOUND
- .dockstore.yml updated: FOUND (6 lines inserted)
- Commit 4893eb93: FOUND
- Commit d3deba25: FOUND
