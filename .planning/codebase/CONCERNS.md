# Concerns

## Technical Debt

### Ubiquitous DEBUG Logging
- **Severity:** Low (operational noise)
- **Where:** 60+ occurrences across task files (`pipes/WDL/tasks/*.wdl`)
- **Issue:** `loglevel` parameters or debug flags scattered throughout tasks, some causing issues (e.g., `stray loglevel param` commit `2ffcc5a6`)
- **Warning sign:** Tasks emit excessive log output in production runs
- **Mitigation:** Standardize log level as a runtime parameter with a default of `WARNING`

### Unresolved Kallisto Sequence ID Tagging
- **Severity:** Medium (workflow correctness)
- **Where:** `classify_kallisto_single.wdl` / `classify_kallisto_multi.wdl` workflow documentation
- **Issue:** TODO in workflow documentation for sequence ID tagging — may affect output traceability
- **Phase:** Should be addressed before using Kallisto outputs in downstream analyses

### Large Monolithic Task Files
- **Severity:** Low (maintainability)
- **Largest files:**
  - `tasks_nextstrain.wdl` (76.7K)
  - `tasks_ncbi.wdl` (72.7K)
  - `tasks_assembly.wdl` (55.4K)
  - `tasks_demux.wdl` (44.1K)
- **Issue:** Large files make it difficult to locate specific tasks, increase merge conflict risk
- **Mitigation:** Consider splitting by sub-domain (e.g., `tasks_ncbi_genbank.wdl`, `tasks_ncbi_biosample.wdl`)

### Hard-coded Docker Versions Requiring Manual Sync
- **Severity:** Medium (reproducibility risk)
- **Where:** `requirements-modules.txt` + individual task `runtime` blocks
- **Issue:** Docker versions must be manually kept in sync between `requirements-modules.txt` and task WDL files
- **Script:** `github_actions_ci/check-wdl-runtimes.sh` detects drift, but doesn't fix it
- **Warning sign:** CI failure on `check-wdl-runtimes.sh` with version mismatch

## Performance Concerns

### Demultiplexing OOM Errors
- **Severity:** High (workflow failures in production)
- **Where:** `tasks_demux.wdl` — demultiplexing tasks
- **Issue:** OOM errors with hardcoded thread limits when processing large flowcells
- **Mitigation:** Memory/CPU parameters should be made workflow inputs with sensible defaults; consider `memory_multiplier` pattern already used elsewhere

### Archive Extraction Bottlenecks
- **Severity:** Medium (throughput)
- **Where:** `unpack_archive_to_bucket.wdl` and archive-handling tasks in `tasks_utils.wdl`
- **Issue:** Sequential upload patterns after extraction slow down large archive processing
- **Mitigation:** Parallelize upload using scatter/gather where file list is known upfront

## Test Coverage Gaps
- **Severity:** Medium (regression risk)
- **Where:** 89 of 101 workflows have no automated integration tests
- **Tested:** 12 workflows via miniWDL, ~3 via DNAnexus
- **Risk areas:** Complex workflows like `sarscov2_illumina_full`, all Nextstrain augur workflows beyond `augur_from_assemblies`
- **Mitigation:** Prioritize adding test inputs for high-frequency production workflows first

## Security Considerations

### DNAnexus TAR Unpacking
- **Severity:** Medium (supply chain risk)
- **Where:** `tasks_terra.wdl`, `unpack_archive_to_bucket.wdl`
- **Issue:** Extracting user-provided TAR files to cloud buckets without content validation creates path traversal risk if archive originates from untrusted sources
- **Mitigation:** Only unpack archives from trusted internal sources; add checksum verification

### Docker Images from Multiple Registries
- **Where:** `requirements-modules.txt` references `broadinstitute/`, `quay.io/staphb/`, `nextstrain/`, `andersenlabapps/`
- **Issue:** Images from third-party registries (`quay.io/staphb/pangolin`, `andersenlabapps/ivar`) are not under Broad's control
- **Mitigation:** Mirror critical images to internal registry; pin by digest (`@sha256:...`) not tag

## Fragile Areas

### WDL Import Path Sensitivity
- **Where:** All `pipes/WDL/workflows/*.wdl` files use relative imports (`../tasks/tasks_*.wdl`)
- **Issue:** Moving any WDL file breaks imports; platform-specific path handling differs between Terra, DNAnexus, Cromwell
- **Script:** `github_actions_ci/relative-wdl-paths.sh` validates paths, but only catches absolute paths
- **Warning sign:** `womtool validate` failures after reorganization

### PAF Pipeline Naming Instability
- **Where:** `align_and_generate_PAF.wdl` and related tasks
- **Issue:** Recent commits (3 of last 5) involved renaming, fixing commas, and correcting output capture for PAF pipeline — suggests this workflow is under active development and may be unstable
- **Recommendation:** Add integration tests before relying on PAF outputs in production

### Pixi Environment Lock
- **Where:** `pixi.lock` (7.6K), `pixi.toml`
- **Issue:** Only `osx-arm64` platform declared — CI likely uses a different platform (Linux x86_64)
- **Risk:** Environment reproducibility broken for Linux CI runners
- **Mitigation:** Add `linux-64` to `[workspace] platforms` in `pixi.toml`

## Known Issues from Recent Git History

| Commit | Issue Fixed | Risk Area |
|--------|-------------|-----------|
| `eb78d172` | Incorrect PAF output capture for `BamToPAF` | PAF pipeline |
| `2ffcc5a6` | Stray `loglevel` param causing issues | Task parameter hygiene |
| `56a05282` | PAF pipeline naming | Workflow organization |
| `68b581f1` | Stray comma in PAF task | WDL syntax fragility |
