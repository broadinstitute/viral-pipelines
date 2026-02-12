# Testing Patterns

**Analysis Date:** 2026-02-11

## Test Framework

**Runner:**
- miniWDL - Primary testing framework for local workflow execution
- Cromwell - Secondary testing framework for compatibility validation
- Both configured in CI: `.github/workflows/build.yml`

**Validation Tools:**
- miniwdl - WDL syntax validation: `miniwdl check pipes/WDL/workflows/*.wdl`
- womtool - Cromwell WDL validator: `womtool validate workflow.wdl`

**Run Commands:**
```bash
# Run all tests with miniWDL
github_actions_ci/tests-miniwdl.sh

# Run all tests with Cromwell
github_actions_ci/tests-cromwell.sh

# Validate WDL syntax (miniwdl)
miniwdl check pipes/WDL/workflows/*.wdl

# Validate WDL syntax (womtool)
github_actions_ci/validate-wdl-womtool.sh

# Run single workflow with miniWDL
miniwdl run pipes/WDL/workflows/assemble_refbased.wdl \
  -i test/input/WDL/miniwdl-local/test_inputs-assemble_refbased-local.json

# Run single task with miniWDL
miniwdl run --task task_name pipes/WDL/tasks/tasks_assembly.wdl \
  input1=value1 \
  input2=value2

# Check required inputs for any workflow
miniwdl run pipes/WDL/workflows/assemble_refbased.wdl

# Run miniWDL self-test
miniwdl run_self_test
```

## Test File Organization

**Location:**
- Test inputs organized by execution engine:
  - `test/input/WDL/miniwdl-local/` - miniWDL test inputs
  - `test/input/WDL/cromwell-local/` - Cromwell test inputs
  - `test/input/WDL/` - Shared test inputs (Cromwell uses these by default)
  - `test/input/` - Test data files (FASTA, BAM, etc.)

**Naming:**
- Input files: `test_inputs-{workflow_name}-local.json`
- Output files: `test_outputs-{workflow_name}-local.json`
- Test data: Descriptive names (e.g., `G5012.3.testreads.bam`, `ebov-makona.fasta`)

**Structure:**
```
test/
├── input/
│   ├── G5012.3.testreads.bam        # Test BAM data
│   ├── ebov-makona.fasta            # Test reference genome
│   └── WDL/
│       ├── miniwdl-local/
│       │   ├── test_inputs-assemble_refbased-local.json
│       │   ├── test_outputs-assemble_refbased-local.json
│       │   ├── test_inputs-assemble_denovo-local.json
│       │   └── test_outputs-assemble_denovo-local.json
│       └── cromwell-local/
│           └── test_inputs-assemble_refbased-local.json
```

**Coverage:**
- Not all workflows have tests (only 12 miniWDL test configs found)
- Key workflows tested: `assemble_refbased`, `assemble_denovo`, `classify_kraken2`, `augur_from_assemblies`, `load_illumina_fastqs`, `sarscov2_lineages`, `genbank_single`

## Test Structure

**Input JSON Format:**
```json
{
  "workflow_name.input_parameter": "value",
  "workflow_name.reads_unmapped_bams": ["test/input/G5012.3.testreads.bam"],
  "workflow_name.reference_fasta": "test/input/ebov-makona.fasta",
  "workflow_name.sample_name": "G5012.3"
}
```

**Expected Output JSON Format:**
```json
{
  "workflow_name.output_parameter": "expected_value"
}
```

**Workflow Test Execution (miniWDL):**
```bash
# From github_actions_ci/tests-miniwdl.sh
for workflow in ../pipes/WDL/workflows/*.wdl; do
  workflow_name=$(basename $workflow .wdl)
  input_json="test/input/WDL/miniwdl-local/test_inputs-$workflow_name-local.json"
  expected_output_json="test/input/WDL/miniwdl-local/test_outputs-$workflow_name-local.json"
  if [ -f $input_json ]; then
    # Run workflow
    miniwdl run -i $input_json -d $workflow_name/. --error-json --verbose $workflow

    # Validate outputs if expected outputs exist
    if [ -f $expected_output_json ]; then
      # Compare expected vs actual outputs
      diff expected actual
    fi
  fi
done
```

**Output Validation:**
- Extract subset of outputs matching expected output keys
- Use `jq` to parse JSON: `cat outputs.json | jq -r '.["key"]'`
- Compare with `diff` - non-zero exit fails test

## Mocking

**Framework:** No explicit mocking framework

**Patterns:**
- Use small real test data instead of mocks
- Test data designed to complete quickly (e.g., `G5012.3.testreads.bam`)
- Tasks support test-friendly parameters:
  - `spades_n_reads = 10000000` - Subsample for faster assembly
  - `always_succeed = false` - Allow tasks to succeed with empty output for testing

**What to Mock:**
- Not applicable - integration testing approach

**What NOT to Mock:**
- Real bioinformatics tools run in Docker containers
- Actual file I/O operations
- Real alignment/assembly operations

## Fixtures and Factories

**Test Data:**
- Located in `test/input/`
- Real biological data files (BAM, FASTA)
- Small enough for CI execution (completes in minutes)

**Example Test Data Files:**
```
test/input/G5012.3.testreads.bam      # Small BAM file for assembly testing
test/input/ebov-makona.fasta          # Reference genome for alignment testing
```

**Location:**
- `test/input/` - All test fixtures
- No separate factory patterns - JSON configs define workflow inputs directly

## Coverage

**Requirements:** None enforced programmatically

**Tested Workflows:**
- 12 workflows have miniWDL test configurations
- Primary workflows tested: assembly, classification, lineage calling, data loading
- Many utility workflows lack test coverage

**View Coverage:**
```bash
# List workflows with tests
ls test/input/WDL/miniwdl-local/test_inputs-*-local.json

# Count tested vs total workflows
echo "Tested: $(ls test/input/WDL/miniwdl-local/test_inputs-*-local.json | wc -l)"
echo "Total: $(ls pipes/WDL/workflows/*.wdl | wc -l)"
```

## Test Types

**Integration Tests:**
- Primary test type - full workflow execution with real tools
- Each test runs complete workflow from inputs to outputs
- Validates correct task composition and data flow
- Docker containers provide isolated, reproducible environments

**Unit Tests:**
- Individual task testing via `miniwdl run --task`
- Used for debugging specific tasks during development
- Not systematically run in CI

**Validation Tests:**
- WDL syntax validation with `miniwdl check` (all workflows)
- WDL syntax validation with `womtool validate` (all workflows)
- Docker version validation: `github_actions_ci/check-wdl-runtimes.sh`
- Documentation build test

**E2E Tests:**
- Not distinguished from integration tests
- Workflows are end-to-end by nature (reads → analysis → results)

## Common Patterns

**Test Isolation:**
```bash
# Create clean test directory
starting_dir="$(pwd)"
test_dir="miniwdl_testing"
mkdir $test_dir
cp -r test $test_dir
cd $test_dir

# Cleanup on exit
function cleanup(){
  cd "$starting_dir"
  if [ -d "$test_dir" ]; then
    rm -r "$test_dir"
  fi
}
trap cleanup EXIT SIGINT SIGQUIT SIGTERM
```

**Docker Image Management:**
```bash
# Prune images after each workflow to save disk space
docker image prune --all --force
```

**Output Validation:**
```bash
# Check for successful execution
if [ -f $workflow_name/outputs.json ]; then
  echo "$workflow_name SUCCESS"
else
  echo "$workflow_name FAILED"
  exit 1
fi
```

**miniWDL Specific:**
```bash
# Run workflow with fixed output directory (no timestamp)
miniwdl run -i $input_json -d $workflow_name/. --error-json --verbose $workflow

# outputs.json existence guarantees successful execution
[ -f $workflow_name/outputs.json ]
```

**Cromwell Specific:**
```bash
# Run with local configuration
java -Dconfig.file=../pipes/cromwell/cromwell.local-github_actions.conf \
  -jar cromwell.jar run workflow.wdl -i input.json

# Check pipe status for errors
if [ ${PIPESTATUS[0]} -gt 0 ]; then
  # Extract stderr logs for debugging
  error_logs=$(grep stderr cromwell.out | perl -lape 's/.*\s(\S+)$/$1/g')
  for log in $error_logs; do
    cat `dirname $log`/stderr*
    cat `dirname $log`/stdout*
  done
  exit 1
fi
```

## CI Configuration

**GitHub Actions Workflow:** `.github/workflows/build.yml`

**Jobs:**
1. `validate_wdl_miniwdl` - Syntax validation with miniwdl
2. `validate_wdl_womtool` - Syntax validation with womtool
3. `test_docs` - Documentation build test
4. `test_miniwdl` - Integration tests with miniWDL
5. `test_cromwell` - Integration tests with Cromwell
6. `deploy_dnanexus` - Deploy to DNAnexus (master/releases only)

**Setup:**
- Python 3.10 for miniWDL
- Java 17 for Cromwell/womtool
- Docker for task execution
- shellcheck for bash script validation

**Test Triggers:**
- All PRs
- All pushes to master
- Release creation

## Adding New Tests

**For New Workflows:**
1. Create test input JSON: `test/input/WDL/miniwdl-local/test_inputs-{workflow_name}-local.json`
2. Follow naming convention: `test_inputs-{workflow_name}-local.json`
3. Optionally create expected outputs: `test_outputs-{workflow_name}-local.json`
4. Test locally: `miniwdl run pipes/WDL/workflows/{workflow_name}.wdl -i test/input/WDL/miniwdl-local/test_inputs-{workflow_name}-local.json`
5. CI automatically picks up tests matching the naming pattern

**For Task Changes:**
1. Identify workflows importing the modified task
2. Run integration tests for those workflows
3. Use `miniwdl run --task` to test task in isolation during debugging
4. Verify outputs match expected results

**Test Data Requirements:**
- Small enough to complete in CI (minutes, not hours)
- Real biological data (not synthetic)
- Representative of actual use cases
- Committed to repository in `test/input/`

## Known Testing Gaps

**Missing Tests:**
- Many utility workflows lack test configurations
- No systematic unit test coverage for tasks
- No performance regression tests
- No security testing

**Manual Testing Required:**
- Platform-specific features (Terra, DNAnexus, AWS)
- Large-scale data processing
- Performance characteristics
- Resource usage optimization

---

*Testing analysis: 2026-02-11*
