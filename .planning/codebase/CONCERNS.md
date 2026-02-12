# Codebase Concerns

**Analysis Date:** 2026-02-11

## Tech Debt

**WDL 1.0/1.1 Version Rollback Confusion:**
- Issue: Repository was upgraded to WDL 1.1 then rolled back to 1.0, but 6 task files still declare `version 1.1`
- Files: `pipes/WDL/tasks/tasks_ncbi_tools.wdl`, `pipes/WDL/tasks/tasks_intrahost.wdl`, `pipes/WDL/tasks/tasks_demux.wdl`, `pipes/WDL/tasks/tasks_reports.wdl`, `pipes/WDL/tasks/tasks_megablast.wdl`, `pipes/WDL/tasks/tasks_interhost.wdl`
- Impact: Version inconsistency could cause validation or compatibility issues with different WDL engines; confuses developers about which version is actually supported
- Fix approach: Change version headers to `version 1.0` in these 6 files to match the rest of the codebase

**Docker Registry Migration Incomplete:**
- Issue: Recent migration from quay.io to ghcr.io for viral-ngs images, but several tasks still use quay.io for other images
- Files: `pipes/WDL/tasks/tasks_16S_amplicon.wdl` (qiime2), `pipes/WDL/tasks/tasks_nextstrain.wdl` (snp-sites), `pipes/WDL/tasks/tasks_ncbi_tools.wdl` (ncbi-tools), `pipes/WDL/tasks/tasks_intrahost.wdl` (polyphonia), `pipes/WDL/tasks/tasks_demux.wdl` (py3-bio), `pipes/WDL/tasks/tasks_reports.wdl` (py3-bio, multiqc)
- Impact: Potential rate limiting from quay.io; inconsistent registry usage across codebase; some images may become unavailable if not also mirrored to ghcr.io
- Fix approach: Mirror all required images to ghcr.io and update WDL task references; document registry strategy

**PAPIv2 Deprecation Warning:**
- Issue: Comment in code notes "PAPIv2 is deprecated and will be removed in the future"
- Files: `pipes/WDL/tasks/tasks_terra.wdl`
- Impact: Future Terra backend changes may break workflows; no migration plan evident
- Fix approach: Monitor Terra announcements for PAPIv2 sunset timeline; prepare migration to newer backend when available

**Limited Test Coverage:**
- Issue: Only 19 test input JSON files exist for 93 workflow files (20% test coverage)
- Files: `test/input/WDL/miniwdl-local/` has 19 test files; `pipes/WDL/workflows/` has 93 workflows
- Impact: Many workflows lack automated testing; changes could introduce regressions that won't be caught in CI
- Fix approach: Add basic test cases for untested workflows following existing patterns; prioritize frequently-used workflows

**Hardcoded Resource Allocation:**
- Issue: Fixed disk sizes and memory allocations throughout tasks (e.g., 375GB disk, 32GB RAM defaults)
- Files: `pipes/WDL/tasks/tasks_assembly.wdl` (disk_size = 375), many tasks in `pipes/WDL/tasks/tasks_read_utils.wdl`, `pipes/WDL/tasks/tasks_nextstrain.wdl`
- Impact: Over-provisioning wastes resources and money; under-provisioning causes task failures; difficult to optimize across different sample sizes
- Fix approach: Implement dynamic resource scaling based on input file sizes; allow users to override defaults

## Known Bugs

**No explicit bugs found in code comments or issue markers**

## Security Considerations

**Unauthenticated GCS Copy Task:**
- Risk: `gcs_copy` task only works on Terra without additional authentication
- Files: `pipes/WDL/tasks/tasks_terra.wdl` (lines 3-32)
- Current mitigation: Task metadata warns "only works on Terra"; TODO comment acknowledges missing credential input
- Recommendations: Add optional GCP credentials input parameter to enable use outside Terra; document authentication requirements clearly

**Docker Image Version Pinning:**
- Risk: Some images use `:latest` tag instead of pinned versions
- Files: `pipes/WDL/tasks/tasks_16S_amplicon.wdl` (qiime2:latest), `pipes/WDL/tasks/tasks_intrahost.wdl` (polyphonia:latest)
- Current mitigation: Most images properly pinned via `requirements-modules.txt`
- Recommendations: Pin all Docker image versions for reproducibility and security; update `requirements-modules.txt` to include all images

## Performance Bottlenecks

**Tar Archive Unpacking Not Parallelized:**
- Problem: Archive extraction and upload to GCS happens sequentially
- Files: `pipes/WDL/tasks/tasks_utils.wdl` (lines 162-174)
- Cause: TODO comment explicitly notes parallelization opportunity not implemented
- Improvement path: Use GNU parallel or background tar processes to extract and upload files simultaneously; requires memory usage testing

**Single-Threaded CI:**
- Problem: Cromwell CI configured with `concurrent-job-limit = 1`
- Files: `pipes/cromwell/cromwell.local-github_actions.conf` (line 15)
- Cause: Likely resource constraints in CI environment
- Improvement path: If CI runners have sufficient resources, increase concurrency to speed up test suite

## Fragile Areas

**Cromwell Docker Container Timeout Workaround:**
- Files: `pipes/cromwell/cromwell.local-github_actions.conf` (lines 55-73)
- Why fragile: Custom `submit-docker` implementation works around detached container bug; uses `docker wait` and manual cleanup to get real exit codes
- Safe modification: Do not modify submit-docker section without understanding the timeout issue it addresses; test thoroughly with long-running tasks
- Test coverage: CI tests exercise this configuration

**WDL Version Rollback History:**
- Files: Multiple files upgraded to 1.1 then rolled back to 1.0 in recent commits (see git log: 74a291be, 3f79542b, d8b80659)
- Why fragile: Type coercion issues, optional type handling, string interpolation requirements were all fixed during 1.1 migration then reverted
- Safe modification: Exercise caution with type conversions, optional parameters, and string concatenation; these areas had breaking changes in 1.1
- Test coverage: Full WDL validation suite runs in CI with both miniwdl and womtool

**Memory Peak Detection Code:**
- Files: `pipes/WDL/tasks/tasks_assembly.wdl` (lines 90), `pipes/WDL/tasks/tasks_read_utils.wdl`, others
- Why fragile: Relies on specific cgroup paths that vary by kernel/container runtime; falls back to "0" if all paths fail
- Safe modification: Test memory reporting across different execution backends (Terra, miniWDL, Cromwell local)
- Test coverage: Outputs max_ram_gb but no validation of accuracy

## Scaling Limits

**Test Data Size Limitations:**
- Current capacity: Test inputs are small samples; no large-scale integration tests
- Limit: Cannot validate performance or correctness at production data sizes in CI
- Scaling path: Add performance benchmarking suite with larger test datasets; run periodically rather than on every commit

**Workflow Count vs. CI Time:**
- Current capacity: 93 workflows with 19 tested configurations
- Limit: Adding more tests increases CI time linearly; already requires concurrent-job-limit=1 to avoid resource exhaustion
- Scaling path: Implement selective testing (only test workflows affected by changes); use more powerful CI runners; parallelize across multiple runners

## Dependencies at Risk

**Quay.io Images Without ghcr.io Mirrors:**
- Risk: Several critical images only referenced from quay.io registry
- Impact: Rate limiting (100 pulls/hour unauthenticated); potential registry downtime affects builds
- Migration plan: Mirror qiime2, ncbi-tools, polyphonia, py3-bio, and biocontainers images to ghcr.io; update task files to use ghcr.io by default

**External Registry Dependencies:**
- Risk: Reliance on third-party registries (quay.io/biocontainers, quay.io/staphb, nextstrain/*, andersenlabapps/*)
- Impact: External organizations control image availability and updates
- Migration plan: Consider mirroring critical external images to ghcr.io/broadinstitute for reliability

## Missing Critical Features

**No Authentication Documentation for Offline Use:**
- Problem: GCS copy and other Terra-specific tasks lack documentation for local/offline execution
- Blocks: Running certain workflows outside Terra environment; reproducibility for external users

**No Resource Estimation Tool:**
- Problem: Users must guess appropriate memory/disk/CPU values for their data
- Blocks: Efficient resource utilization; frequently leads to over-provisioning or task failures

**Limited Retry Configuration:**
- Problem: Most tasks use `maxRetries: 2` without differentiation by failure type
- Blocks: Efficient handling of preemption (which should retry) vs. actual errors (which shouldn't)

## Test Coverage Gaps

**Untested Workflows:**
- What's not tested: 74 of 93 workflows lack integration tests (workflows without corresponding test JSON in `test/input/WDL/miniwdl-local/`)
- Files: Major gaps include many workflows in `pipes/WDL/workflows/` such as metagenomic pipelines, submission workflows, and specialty analysis workflows
- Risk: Changes to shared tasks could break these workflows without detection; no regression testing
- Priority: High - add tests for frequently-used workflows first

**No Unit Testing for Individual Tasks:**
- What's not tested: Task logic is only tested via full workflow integration tests
- Files: All tasks in `pipes/WDL/tasks/` lack isolated unit tests
- Risk: Task changes require running full workflows to validate; slow feedback loop; difficult to test edge cases
- Priority: Medium - implement `miniwdl run --task` tests for complex tasks

**No Performance Regression Tests:**
- What's not tested: Runtime performance, memory usage, or resource efficiency
- Files: Tasks output `max_ram_gb`, `runtime_sec`, `cpu_load_15min` but no validation
- Risk: Performance regressions go unnoticed; resource allocation changes could increase costs without detection
- Priority: Low - implement benchmarking suite for critical path tasks

---

*Concerns audit: 2026-02-11*
