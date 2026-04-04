# External Integrations

**Analysis Date:** 2026-03-31

## APIs & External Services

**National Center for Biotechnology Information (NCBI):**
- GenBank sequence download and submission
- SRA (Sequence Read Archive) submission
- BioSample registration
- Taxonomy database access
- SDK/Client: `broadinstitute/ncbi-tools:2.11.1` Docker image
- Implementation: Tasks in `pipes/WDL/tasks/tasks_ncbi.wdl`, `tasks_ncbi_tools.wdl`
- Auth: `emailAddress` parameter (user credentials passed to NCBI API calls)
- Key tasks:
  - `download_fasta()` - Fetch sequences from GenBank accessions
  - `download_annotations()` - Download feature tables and annotations
  - `prepare_genbank_submission()` - Prepare data for GenBank submission (files: `pipes/WDL/tasks/tasks_ncbi.wdl` lines 1249+)
  - `validate_submission_table2asn()` - Validate GenBank submissions
  - `prepare_sra_submission()` - Prepare for SRA submission (file: `pipes/WDL/tasks/tasks_ncbi.wdl` line 621)

**Nextstrain/Augur Phylogenetics:**
- Phylogenetic tree generation and visualization
- Clade assignment via Nextclade
- SDK/Client: `nextstrain/nextclade:3.18.1`, `nextstrain/base:build-20240318T173028Z`
- Implementation: Tasks in `pipes/WDL/tasks/tasks_nextstrain.wdl`
- Key tasks:
  - `nextclade_one_sample()` - Single sample clade assignment (line 53)
  - `nextclade_many_samples()` - Batch clade assignment (line 144)
  - Dataset management via Nextclade public dataset API

**Pangolin (SARS-CoV-2 Lineage Assignment):**
- Lineage determination for SARS-CoV-2 samples
- SDK/Client: `quay.io/staphb/pangolin:4.3.3-pdata-1.36`
- Implementation: Tasks in `pipes/WDL/tasks/tasks_sarscov2.wdl` (lines 13, 96)
- Key tasks:
  - `pangolin_classify()` - Assign lineages to consensus sequences

**GISAID (Global Initiative on Sharing Avian Influenza Data):**
- Sequence data submission to GISAID database
- SDK/Client: `quay.io/broadinstitute/gisaid-cli:3.0`
- Implementation: Task in `pipes/WDL/tasks/tasks_sarscov2.wdl` (line 564)
- Auth: Credentials passed via GISAID CLI
- Key task: `gisaid_uploader()` - Uploads FASTA and metadata CSV to GISAID

**Dockstore (GA4GH Tool Registry Service):**
- Workflow registration and distribution
- Configuration: `.dockstore.yml` - Registers 100+ workflows
- Purpose: Enables easy import to Terra, DNAnexus, and other platforms supporting TRS API
- Deployment: Continuous deployment on master branch or releases

## Data Storage

**Databases:**
- NCBI Taxonomy Database
  - Connection: Downloaded as tarball in tasks
  - Access: Via taxonomy Python module from `broadinstitute/viral-ngs`
  - Used for: Taxid-to-dataset mapping, taxonomy hierarchy queries

- Nextclade Reference Datasets
  - Connection: Public API for dataset download
  - Access: `nextclade dataset get` command
  - Used for: Clade assignment, sequence annotations

- GenBank/SRA Public Databases
  - Connection: NCBI Entrez API
  - Auth: User email address (for courtesy logging)
  - Used for: Reference genome download, sequence submission

**File Storage:**
- GCS (Google Cloud Storage) - Primary on Terra/GCP
  - Implementation: `tasks_terra.wdl` contains `gcs_copy` task (line 3)
  - Method: `gcloud storage cp` (works natively only on Terra)
  - Used for: Output staging, workspace buckets

- Local filesystem - For local execution
  - Cromwell backend: `pipes/cromwell/cromwell.local-travis.conf`
  - Root directory: `cromwell-executions`

- DNAnexus Platform Storage
  - Via DNAnexus applet launcher (`pipes/dnax/dx-launcher/`)
  - Used for: Demultiplexing results, workflow outputs

**Caching:**
- Docker layer caching - Via Docker image registries (GHCR, Quay.io)
- No explicit caching service integration detected
- Workflow-level caching handled by Cromwell call caching

## Authentication & Identity

**Auth Provider:**
- Custom implementation (no unified auth system)

**Mechanisms by Service:**
- **Terra/GCP:** OAuth2 via gcloud CLI and FireCloud API
  - Detected via: `check_terra_env` task (line 34 in `tasks_terra.wdl`)
  - Credentials: Automatic via service account when running on Terra
  - Uses: `gcloud auth print-access-token` for API calls
  - FireCloud API endpoint: `https://api.firecloud.org/me?userDetailsOnly=true`

- **NCBI:** Email-based courtesy authentication
  - Parameter: `emailAddress` (passed to all NCBI tasks)
  - Implementation: Direct to NCBI Entrez API

- **GISAID:** Credentials via CLI
  - Method: GISAID username/password via gisaid-cli
  - Auth file location: Not hardcoded (passed at runtime)

- **DNAnexus:** API token
  - Implementation: Optional API token file for independent job launching
  - Location: Separate DNAnexus project for token isolation (see `pipes/dnax/dx-launcher/README.md`)

## Monitoring & Observability

**Error Tracking:**
- Not detected - No formal error tracking service integration

**Logs:**
- Stdout/Stderr capture: All WDL tasks capture logs to standard output
- Terra GCS logs: Automatically available in workspace buckets
- Local execution: Cromwell writes logs to `cromwell-executions/` directory
- Log analysis pattern documented in `AGENTS.md` (lines 204-327):
  - Batch job timing queries via `get_batch_job_status`
  - Log parsing from stderr for task timing
  - GCS file modification timestamps for accurate execution duration

## CI/CD & Deployment

**Hosting:**
- GitHub (source code repository)
- ReadTheDocs (documentation - `http://viral-pipelines.readthedocs.io`)
- Dockstore (workflow registry and distribution)
- DNAnexus CI project (continuous deployment)
- Docker registries: GHCR (GitHub Container Registry), Quay.io

**CI Pipeline:**
- GitHub Actions (`.github/workflows/build.yml`)
- Triggers: Push, pull request, release, merge_group
- Jobs:
  1. `validate_wdl_miniwdl` - WDL syntax check with miniwdl (line 45)
  2. `validate_wdl_womtool` - WDL validation with Cromwell's womtool
  3. `test_docs` - Sphinx documentation build
  4. `test_miniwdl` - Integration tests with miniWDL (line ~150)
  5. `test_cromwell` - Integration tests with Cromwell
  6. `deploy_dnanexus` - Deploy to DNAnexus on master/releases

**Deployment Targets:**
- DockerHub/Quay.io - Docker images pushed on releases
- Dockstore - Workflows registered automatically
- DNAnexus - dx-launcher applet deployed
- ReadTheDocs - Documentation auto-built on push/PR

## Environment Configuration

**Required env vars:**
- None hardcoded (platform-specific handling)
- For local testing: Can use `.env` equivalent but not tracked
- For Terra: Automatic service account credentials
- For NCBI tasks: `emailAddress` parameter (required input, not env var)
- For GISAID: Credentials passed via CLI interface

**Secrets location:**
- Not in repository (`.env*` and `secrets/` in `.gitignore`)
- Terra: Handled by service account (no user management)
- DNAnexus: API token stored in separate project (see `pipes/dnax/dx-launcher/README.md`)
- GISAID: Credentials passed at job submission time
- NCBI: Email passed as parameter

## Webhooks & Callbacks

**Incoming:**
- DNAnexus dx-streaming-upload - Triggers dx-launcher applet on sequence run upload completion
  - Endpoint: DNAnexus project-specific
  - Triggers: Demux_launcher automatically launches on run completion
  - See: `pipes/dnax/dx-launcher/README.md`

**Outgoing:**
- GISAID uploader - Submits sequences to GISAID
  - Destination: GISAID REST API
  - Triggered by: `gisaid_uploader` task in `tasks_sarscov2.wdl`

- NCBI Submission - GenBank and SRA submission
  - Destination: NCBI FTP/API
  - Triggered by: GenBank and SRA submission tasks in `tasks_ncbi.wdl`

## Integration Patterns

**Data Flow for GenBank Submission:**
1. Assemble workflow generates consensus FASTA files
2. `download_annotations()` fetches feature tables from reference genomes (GenBank)
3. `annot_transfer()` transfers annotations to new assemblies
4. `prepare_genbank_submission()` creates submission files (SBT template + source table)
5. `genbank_submission()` packages for FTP submission to NCBI

**Data Flow for SARS-CoV-2 Lineage Assignment:**
1. Nextclade task downloads latest dataset via public API
2. Pangolin task classifies samples
3. Results used in nextstrain workflows for phylogenetic visualization

**Data Flow for Terra Integration:**
1. Workflows executed on Terra platform via Dockstore import
2. `check_terra_env` task detects Terra/GCP environment
3. `gcs_copy` task stages outputs to workspace bucket
4. Metadata queries via FireCloud API (credentials automatic)

## Service Account & Credential Management

**Terra Service Account:**
- Automatic provisioning when workflows run on Terra
- Credentials passed via environment variables by Cromwell runtime
- Scoped to workspace and project permissions
- Used for: GCS access, FireCloud API calls

**Local Testing:**
- Credentials not required for local execution with miniWDL
- Docker images pulled from public registries

**Platform-Specific:**
- GCP/Terra: OAuth2 service account (automatic)
- DNAnexus: Manual API token (user-provided)
- NCBI: Email-based (no credential required, courtesy parameter)
- GISAID: Username/password via CLI interface

---

*Integration audit: 2026-03-31*
