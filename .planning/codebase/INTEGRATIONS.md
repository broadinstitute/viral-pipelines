# External Integrations

**Analysis Date:** 2026-02-11

## APIs & External Services

**NCBI APIs:**
- GenBank sequence download - Downloads genomic sequences via NCBI E-utilities
  - SDK/Client: viral-ngs ncbi module (in Docker images)
  - Auth: Email address required (passed as parameter)
  - Files: `pipes/WDL/tasks/tasks_ncbi.wdl`, `pipes/WDL/tasks/tasks_ncbi_tools.wdl`
  - Tasks: `download_fasta`, `download_fasta_from_accession_string`, `fetch_sra_to_bam`

- BioSample/SRA/GenBank submission - Submits sequences and metadata to NCBI databases
  - SDK/Client: viral-ngs ncbi module
  - Auth: NCBI submission credentials (passed via input files)
  - Files: `pipes/WDL/tasks/tasks_ncbi.wdl`, `pipes/WDL/workflows/submit_*.wdl`
  - Workflows: `submit_biosample`, `submit_genbank`, `submit_sra`

**GISAID API:**
- GISAID data ingest - Downloads and sanitizes SARS-CoV-2 sequences from GISAID
  - SDK/Client: Nextstrain ncov tools (in nextstrain/base Docker image)
  - Auth: GISAID credentials (workflow-specific)
  - Files: `pipes/WDL/workflows/sarscov2_gisaid_ingest.wdl`

**Nextstrain Services:**
- Phylogenetic tree visualization - Generates Nextstrain/Auspice JSON outputs
  - SDK/Client: Augur toolkit (in nextstrain/base:build-20240318T173028Z)
  - Auth: None required
  - Files: `pipes/WDL/tasks/tasks_nextstrain.wdl`
  - Tasks: Multiple augur-based tasks (filter, align, tree, refine, export, etc.)

- Nextclade typing - Viral clade assignment and QC
  - SDK/Client: nextclade CLI (in nextstrain/nextclade:3.18.1)
  - Auth: None required
  - Files: `pipes/WDL/tasks/tasks_sarscov2.wdl`

**Terra/Verily Platform:**
- Workspace data table operations - Read/write data to Terra workspace tables
  - SDK/Client: gcloud CLI, Terra API
  - Auth: Service account/pet account (automatic on Terra)
  - Files: `pipes/WDL/tasks/tasks_terra.wdl`, `pipes/WDL/workflows/terra_*.wdl`
  - Workflows: `terra_table_to_tsv`, `terra_tsv_to_table`, `terra_update_assemblies`
  - Environment detection via metadata.google.internal introspection

**DNAnexus Platform:**
- Workflow execution platform - Deploy and run WDL workflows
  - SDK/Client: dxWDL, dxCompiler, dx-toolkit
  - Auth: DX_API_TOKEN (stored in GitHub Actions secrets)
  - Project: project-F8PQ6380xf5bK0Qk0YPjB17P (specified in `.github/workflows/build.yml`)
  - Files: `pipes/dnax/` configuration files

**Dockstore Registry:**
- Workflow registry - Publish workflows for discoverability
  - Config: `.dockstore.yml`
  - Registered: 107 workflows
  - Auth: Not visible in codebase

## Data Storage

**Databases:**
- None - This is a workflow/pipeline repository, not a data service
  - Data storage handled by execution backends (GCS for Terra/GCP, local filesystem for miniWDL)

**File Storage:**
- **Google Cloud Storage (GCS):**
  - Used by Terra workflows for input/output data
  - Bucket paths: Workspace-specific (e.g., `gs://viral-temp-30d/USERNAME/cromwell-test`)
  - Client: gcloud storage CLI (in ghcr.io/broadinstitute/viral-ngs:3.0.4-baseimage)
  - Tasks: `gcs_copy` in `pipes/WDL/tasks/tasks_terra.wdl`
  - Auth: Service account tokens via `gcloud auth print-access-token`

- **Local Filesystem:**
  - Used by miniWDL and local Cromwell execution
  - No special configuration required

**Caching:**
- Cromwell call caching - Configured in `pipes/cromwell/*.conf` files
  - Duplication strategy: "copy" (specified in cromwell.gcid-viral-seq.conf)
  - Filesystem: GCS for cloud backends

## Authentication & Identity

**Auth Provider:**
- Custom per-platform
  - **Terra/GCP:** Service account ("pet" account) automatic authentication via metadata server
  - **DNAnexus:** API token-based (DX_API_TOKEN environment variable)
  - **NCBI:** Email-based for downloads, credential files for submissions
  - **Local execution:** None required

**Implementation:**
- Platform detection via `check_terra_env` task in `pipes/WDL/tasks/tasks_terra.wdl`
- Metadata introspection: Queries metadata.google.internal to detect GCP/Terra environment
- Conditional authentication: Auth tokens obtained dynamically when running on Terra

## Monitoring & Observability

**Error Tracking:**
- None - Workflow execution errors reported via WDL engine logs

**Logs:**
- Workflow logs captured by execution engine (miniWDL, Cromwell)
- Standard output/error streams from Docker containers
- Task-level stdout() captured as output files in many tasks

## CI/CD & Deployment

**Hosting:**
- **GitHub:** Source code repository (broadinstitute/viral-pipelines)
- **GitHub Container Registry (ghcr.io):** Docker image hosting
  - Primary registry for viral-ngs images (migrated from quay.io)
- **Quay.io:** Legacy Docker image hosting
  - quay.io/broadinstitute/viral-pipelines (production)
  - quay.io/broadinstitute/viral-pipelines-dev (development)
- **ReadTheDocs:** Documentation hosting
  - Configured in `.readthedocs.yml`
  - Build trigger: Automatic on git push
- **Dockstore:** Workflow registry hosting

**CI Pipeline:**
- **GitHub Actions** (`.github/workflows/build.yml`)
  - Triggers: push, pull_request, release, merge_group
  - Jobs:
    - `validate_wdl_miniwdl` - WDL validation with miniwdl
    - `validate_wdl_womtool` - WDL validation with womtool
    - `test_docs` - Documentation build testing
    - `test_cromwell` - Workflow testing with Cromwell
    - `test_miniwdl` - Workflow testing with miniWDL
    - `deploy_dnanexus` - Deploy workflows to DNAnexus platform
  - Platform: ubuntu-24.04 runners
  - Docker support: Docker Buildx configured

## Environment Configuration

**Required env vars:**
- `PYTHONIOENCODING=UTF8` - Set in GitHub Actions
- `DX_API_TOKEN` - DNAnexus authentication (GitHub Actions secret)
- `DX_PROJECT` - DNAnexus project ID for deployment

**Optional env vars:**
- Platform-specific variables detected at runtime:
  - `BATCH_JOB_UID` - Indicates GCP Batch execution
  - Terra workspace variables (workspace_id, workspace_name, etc.) - Auto-populated on Terra
  - `GOOGLE_PROJECT` - GCP project ID (detected via metadata server)

**Secrets location:**
- GitHub Actions Secrets: DX_API_TOKEN for DNAnexus
- NCBI credentials: Provided as workflow input files (not in repository)
- GCP/Terra: Service account keys managed by platform

## Webhooks & Callbacks

**Incoming:**
- None detected

**Outgoing:**
- None detected
- Data release workflows upload to external repositories (NCBI, GISAID) but do not use webhooks

## Reference Data & Databases

**Taxonomic Databases:**
- Kraken2/KrakenUniq databases - Custom built via `kraken2_build` workflow
  - Test database: `test/input/kraken2_db-tinytest.tar.zst`
- Krona taxonomy - Downloaded reference data
  - Test file: `test/input/krona.taxonomy-20200505.tab.zst`
- NCBI taxonomy dump - Reference taxonomic data
  - Test file: `test/input/taxdump-20231214.tar.gz`

**Virus Annotation Databases:**
- VADR (Viral Annotation DefineR) - For GenBank submissions
  - Configuration: `test/input/vadr-by-taxid.tsv`

---

*Integration audit: 2026-02-11*
