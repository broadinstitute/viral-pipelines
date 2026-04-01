# Testing

## Framework Overview

This codebase uses **integration tests only** — no unit tests. Tests execute real WDL workflows end-to-end against small reference inputs using miniWDL or Cromwell.

## Test Engines

### miniWDL (Primary Local Engine)
- **Script:** `github_actions_ci/tests-miniwdl.sh`
- **Command:** `miniwdl run -i <input.json> -d <outdir>/. --error-json --verbose <workflow.wdl>`
- **Self-test:** `miniwdl run_self_test` (verifies Docker swarm mode works)
- **Output validation:** Compares `outputs.json` key-by-key against expected output JSON files

### Cromwell (Secondary Engine)
- **Script:** `github_actions_ci/tests-cromwell.sh`
- **Config:** `pipes/cromwell/cromwell.local-github_actions.conf`
- **Command:** `java -Dconfig.file=... -jar cromwell.jar run <workflow.wdl> -i <input.json>`
- **Error capture:** Parses `cromwell.out` for stderr logs on failure

### DNAnexus (Cloud Platform)
- **Script:** `github_actions_ci/tests-dx.sh`
- **Requires:** DNAnexus credentials and project setup

### WDL Validation (Non-Execution)
- **Script:** `github_actions_ci/validate-wdl-womtool.sh`
- **Tool:** Womtool (`womtool validate`)
- **Purpose:** Syntax validation without execution

## Test Structure

### Test Input Files
Location: `test/input/WDL/miniwdl-local/` (19 JSON files)

**Naming convention:**
- Inputs: `test_inputs-<workflow_name>-<platform>.json`
- Expected outputs: `test_outputs-<workflow_name>-<platform>.json`

**Workflows with miniWDL tests (12 of 101):**
- `align_and_plot`
- `assemble_denovo`
- `assemble_refbased`
- `augur_from_assemblies`
- `augur_from_beast_mcc`
- `beast_gpu` (non-auto/manual only)
- `classify_kraken2`
- `fastq_to_ubam`
- `genbank_single`
- `load_illumina_fastqs`
- `load_illumina_fastqs_deplete`
- `sarscov2_lineages`

**Workflows with DNAnexus tests:**
- `assemble_denovo` (dx)
- `classify_multi` (dx)
- `demux_plus` (dx)

### Test Data
Location: `test/input/`

| File | Description |
|------|-------------|
| `G5012.3.testreads.bam` (1.5M) | Primary test BAM |
| `G5012.3.mini.bam` (372K) | Minimal test BAM |
| `G5012.3.subset.bam` (16.7K) | Subset BAM |
| `G5012.3.fasta` (18.7K) | Reference FASTA |
| `ebov-makona.fasta` (18.8K) | Ebola reference |
| `kraken2_db-tinytest.tar.zst` (43.3K) | Tiny Kraken2 database |
| `krona.taxonomy-20200505.tab.zst` (15.9M) | Krona taxonomy |
| `beast-lasv.tree` (345.2K) | BEAST tree (LASV) |
| `beast-mers-mcc.tree` (266.1K) | BEAST MCC tree (MERS) |
| `clipDb.fasta` (18.5K) | Adapter clipping database |

## Output Validation

miniWDL tests optionally validate outputs against expected JSON:

```bash
for k in $(cat $expected_output_json | jq -r 'keys[]'); do
    echo -n "$k=" >> expected
    echo -n "$k=" >> actual
    cat $expected_output_json       | jq -r '.["'$k'"]' >> expected
    cat $workflow_name/outputs.json | jq -r '.["'$k'"]' >> actual
done
diff expected actual
```

Only keys present in the expected JSON are validated — extra output keys are ignored.

## CI/CD Pipeline

**GitHub Actions:** `.github/workflows/build.yml` (16.6K)

Typical CI flow:
1. `install-wdl.sh` — Install miniWDL, Cromwell, Womtool
2. `validate-wdl-womtool.sh` — Validate all WDL syntax
3. `check-wdl-runtimes.sh` — Verify Docker image versions
4. `tests-miniwdl.sh` — Run integration tests
5. `build-docs.sh` — Build documentation

## Coverage Gaps

- **12 of 101 workflows** have miniWDL integration tests (~12%)
- **0 unit tests** for individual WDL tasks
- No mocking — all tests use real Docker containers
- Large workflows (e.g., `sarscov2_illumina_full`) are not integration tested locally due to resource requirements
- DNAnexus and Terra tests require platform credentials

## Running Tests Locally

```bash
# Install dependencies
bash github_actions_ci/install-wdl.sh

# Validate all WDL files
bash github_actions_ci/validate-wdl-womtool.sh

# Run miniWDL integration tests
bash github_actions_ci/tests-miniwdl.sh

# Run Cromwell tests (requires cromwell.jar)
bash github_actions_ci/tests-cromwell.sh
```

## Python Environment

Managed via **pixi** (`pixi.toml`):
- Python >= 3.14.2, < 3.15 (conda-forge channel)
- `.coveragerc` present — suggests Python coverage tracking for any Python utilities
