---
phase: 01-core-tasks-and-single-sample-workflow
plan: 03
subsystem: testing
tags: [genomad, testing, dockstore, ci-integration, platform-discovery]

# Dependency graph
requires:
  - phase: 01
    plan: 01
    provides: genomad_end_to_end and report_genomad_summary tasks
  - phase: 01
    plan: 02
    provides: genomad_single.wdl workflow
provides:
  - Test input JSON for genomad_single workflow
  - Dockstore registry entry for platform discovery
affects: [ci-testing, terra-integration, platform-users]

# Tech tracking
tech-stack:
  added: [test input JSON, Dockstore registration]
  patterns: [test input naming convention, miniwdl-local symlink pattern]

key-files:
  created:
    - test/input/WDL/test_inputs-genomad_single-local.json
    - test/input/WDL/miniwdl-local/test_inputs-genomad_single-local.json (symlink)
  modified:
    - .dockstore.yml

key-decisions:
  - "Use G5012.3.fasta as test assembly input (existing viral test data)"
  - "Reference genomad_db-tinytest.tar.zst for test database (to be created separately)"
  - "Follow existing pattern: JSON in parent WDL/ directory with symlink in miniwdl-local/"

patterns-established:
  - "Pattern: Test input JSON naming convention test_inputs-{workflow}-local.json"
  - "Pattern: Dockstore testParameterFiles reference parent directory path"
  - "Pattern: Document test database fixtures that need separate creation"

# Metrics
duration: 2min 14sec
completed: 2026-02-12
---

# Phase 01 Plan 03: Testing Infrastructure and Platform Registration Summary

**Created test input JSON and registered genomad_single workflow in Dockstore for CI discovery and platform integration**

## Performance

- **Duration:** 2 min 14 sec
- **Started:** 2026-02-12T17:12:12Z
- **Completed:** 2026-02-12T17:14:26Z
- **Tasks:** 2
- **Files created:** 2 (1 file + 1 symlink)
- **Files modified:** 1

## Accomplishments

- Created test_inputs-genomad_single-local.json following repository naming conventions
- Test JSON references G5012.3.fasta (existing viral genome test data)
- Test JSON points to genomad_db-tinytest.tar.zst (fixture to be created separately)
- Created symlink in miniwdl-local/ directory following repository pattern
- Registered genomad_single in .dockstore.yml in alphabetical order
- All WDL files pass miniwdl check validation
- Verified all meta and parameter_meta blocks exist in tasks and workflow
- Verified Docker image versions match requirements-modules.txt (3.0.4)

## Task Commits

Each task was committed atomically:

1. **Task 1: Create test input JSON for genomad_single workflow** - `ce965c53` (feat)
2. **Task 2: Register genomad_single in Dockstore and verify documentation** - `0be0878e` (feat)

## Files Created/Modified

**Created:**
- `test/input/WDL/test_inputs-genomad_single-local.json` - Test input JSON (4 lines)
- `test/input/WDL/miniwdl-local/test_inputs-genomad_single-local.json` - Symlink to parent directory

**Modified:**
- `.dockstore.yml` - Added genomad_single workflow entry (+5 lines)

## Decisions Made

**Use G5012.3.fasta as test assembly:** G5012.3.fasta is an existing ~19KB viral genome FASTA file in test/input/, making it suitable for genomad classification testing without requiring new test data creation.

**Reference genomad_db-tinytest.tar.zst:** Following the pattern of kraken2_db-tinytest.tar.zst, the test input JSON references a miniature genomad database fixture that needs to be created separately. The test will be discovered by CI but won't execute until this fixture is available.

**Follow miniwdl-local symlink pattern:** Observed that most test input files in miniwdl-local/ are symlinks to the parent WDL/ directory, maintaining a single source of truth while supporting both directory structures.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - test infrastructure created and validated successfully.

## User Setup Required

None for basic workflow usage. For running the test:

**Test database fixture needed:** The genomad_db-tinytest.tar.zst file referenced in the test input JSON needs to be created. Options:
1. Download full genomad database from https://zenodo.org/records/10594875 and subset it
2. Obtain pre-existing test database from genomad project
3. Create minimal test tarball with expected directory structure

Until this fixture is available, the CI test will be skipped (miniWDL won't find the database file). The test input JSON is correct so that when a database IS provided, tests run immediately.

## Next Phase Readiness

Test infrastructure complete. Phase 1 deliverables complete:
- ✓ genomad_end_to_end task (Plan 01-01)
- ✓ report_genomad_summary task (Plan 01-01)
- ✓ genomad_single workflow (Plan 01-02)
- ✓ Test input JSON (Plan 01-03)
- ✓ Dockstore registration (Plan 01-03)

Ready for:
- Phase 02: Batch processing workflows
- Creating genomad test database fixture for CI execution

**Blockers:** None (test database fixture is deferred work, not a blocker)

---
*Phase: 01-core-tasks-and-single-sample-workflow*
*Completed: 2026-02-12*

## Self-Check: PASSED

All claims verified:
- ✓ Created file exists: test/input/WDL/test_inputs-genomad_single-local.json
- ✓ Created symlink exists: test/input/WDL/miniwdl-local/test_inputs-genomad_single-local.json
- ✓ Modified file exists: .dockstore.yml
- ✓ Commit ce965c53 exists (Task 1: Create test input JSON)
- ✓ Commit 0be0878e exists (Task 2: Register in Dockstore)
- ✓ Test JSON is valid JSON with correct workflow key format
- ✓ Test JSON references existing G5012.3.fasta
- ✓ .dockstore.yml contains genomad_single entry with correct paths
- ✓ genomad_single.wdl passes miniwdl check validation
- ✓ All meta and parameter_meta blocks present
- ✓ Docker versions consistent (3.0.4-classify)
