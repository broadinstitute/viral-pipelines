# Stack Research

**Domain:** Genomad Integration for WDL Viral Genomics Pipelines
**Researched:** 2026-02-12
**Confidence:** HIGH

## Recommended Stack

### Core Technologies

| Technology | Version | Purpose | Why Recommended |
|------------|---------|---------|-----------------|
| WDL | 1.0 | Workflow definition language | Existing codebase standard, maximum compatibility across all execution platforms (Terra, DNAnexus, local), fully supported by all engines. WDL 1.1/1.2 add features but reduce platform compatibility. |
| miniWDL | 1.13.0+ | Local/HPC workflow execution and validation | Active development, WDL 1.0/1.1 support, lightweight Python-based runner ideal for local testing and HPC environments, supports Docker/Podman/Singularity for container flexibility. |
| Cromwell | 88+ | Cloud/enterprise workflow execution | Production-grade engine for Terra and cloud platforms, version 88+ includes critical database improvements and Google Batch backend (Cloud Life Sciences API removed July 2025). Version 91+ recommended for most recent fixes. |
| genomad | 1.8+ (database v1.8) | Viral discovery from assembled contigs | State-of-the-art viral/plasmid identification (Nature Biotech 2023), outperforms alternatives (MCC 95.3% for viruses), already installed in viral-ngs:3.0.4-classify Docker image, designed for metagenomics/viral discovery workflows. |

### Supporting Libraries

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| mmseqs2 | 15.6f452+ | Protein sequence search (genomad dependency) | Required by genomad for marker gene detection, automatically included in viral-classify Docker images. |
| pyrodigal-gv | 0.3.1+ | Gene prediction (genomad dependency) | Required by genomad for ORF finding in viral contigs, included in genomad installations. |
| tensorflow | 2.16+ | Deep learning framework (genomad dependency) | Required by genomad for neural network-based classification, must be compatible with keras >=3. |
| krona | (included in viral-classify) | Interactive visualization of taxonomic results | For generating HTML reports from genomad taxonomy assignments, follows existing viral-pipelines pattern. |

### Development Tools

| Tool | Purpose | Notes |
|------|---------|-------|
| miniwdl check | WDL syntax validation | Run on all .wdl files before commit, validates WDL 1.0 compliance |
| miniwdl run | Local workflow testing | Test with --task flag for individual task debugging, use test inputs from test/input/WDL/miniwdl-local/ |
| womtool | Cromwell-based WDL validation | Secondary validation tool, requires Java 11, installed via github_actions_ci/install-wdl.sh |
| gcloud storage | GCS file operations for Terra workflows | Preferred over gsutil (faster, more reliable), use wildcards for batch queries of workflow outputs |

## Installation

### For Workflow Development/Testing

```bash
# Install miniWDL for local development
pip install miniwdl>=1.13.0

# Install Cromwell for platform compatibility testing
github_actions_ci/install-wdl.sh

# Genomad is pre-installed in Docker images - no local installation needed
# Database will be user-provided as workflow input (~5GB)
```

### For Users Running Workflows

**No installation required** - workflows run via:
- **Terra**: Built-in Cromwell server, Docker images auto-pulled
- **DNAnexus**: dxCompiler + native workflow engine, containers managed by platform
- **Local/HPC**: miniWDL with Docker, Podman, or Singularity
- **AWS/GCP/Azure**: Cromwell with cloud backends, managed container execution

## Alternatives Considered

| Recommended | Alternative | When to Use Alternative |
|-------------|-------------|-------------------------|
| WDL 1.0 | WDL 1.1/1.2 | Only if targeting single platform (e.g., Terra-only) AND need new features (Directory type, requirements/hints). Not recommended due to reduced platform compatibility and DNAnexus partial support. |
| miniWDL 1.13.0+ | Cromwell for local testing | If testing cloud-specific features (GCP Batch, AWS Batch) locally. miniWDL preferred for faster iteration. |
| genomad | CheckV alone | If only validating already-identified viral sequences (CheckV is quality control, not discovery). Use both together: genomad for discovery + CheckV for quality assessment. |
| genomad | VirSorter2 | If database size is critical constraint (VirSorter2 smaller DB) OR legacy pipeline compatibility. Genomad has superior accuracy (MCC 95.3% vs ~85% for VirSorter2). |
| Docker (via viral-classify image) | Standalone conda/pip genomad | Only for non-WDL use cases. Docker ensures reproducibility and is required for WDL platform-agnostic design. |
| Cromwell 88+ | Cromwell <88 | Never - versions <88 lack Google Batch support (Cloud Life Sciences shutdown July 2025), older database schema. |

## What NOT to Use

| Avoid | Why | Use Instead |
|-------|-----|-------------|
| WDL draft-2 | Deprecated, unsupported by modern platforms (Terra, DNAnexus), lacks critical features (scatter/gather improvements, struct types) | WDL 1.0 (current codebase standard) |
| gsutil | Slower than gcloud storage, deprecated by Google, less reliable for large operations | gcloud storage (with wildcards for batch operations) |
| Cloud Life Sciences API (Papiv2) | Shut down by Google in July 2025, removed from Cromwell 91+ | Google Batch backend (standard in Cromwell 88+) |
| Docker image tags (e.g., :latest) | Non-reproducible, breaks workflows when images update | SHA256 digests or pinned version tags (e.g., viral-ngs:3.0.4-classify) |
| CheckV for viral discovery | CheckV is for quality control of known viruses, not discovery/classification | genomad for discovery, then CheckV for quality assessment |
| VirFinder/DeepVirFinder | Older ML approaches, lower accuracy, not designed for metagenomes | genomad (newer, more accurate, better metagenomic support) |

## Stack Patterns by Variant

### Single-Sample Genomad Workflow

**Use:**
- WDL 1.0 workflow calling genomad task
- viral-ngs:3.0.4-classify (or later) Docker image
- Input: FASTA file (assembled contigs)
- User-provided database path (File input)

**Pattern:**
```wdl
import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow genomad_single {
    call metagenomics.genomad {
        input:
            contigs_fasta = sample_contigs,
            genomad_db = user_database
    }
}
```

### Multi-Sample Genomad Workflow (Scatter-Gather)

**Use:**
- WDL 1.0 scatter block over array of FASTA inputs
- Same genomad task scattered across samples
- Database loaded once per shard (disk-efficient)
- Krona aggregation for merged visualization

**Pattern:**
```wdl
scatter(sample_fasta in contigs_fastas) {
    call metagenomics.genomad {
        input:
            contigs_fasta = sample_fasta,
            genomad_db = shared_database
    }
}
call metagenomics.krona_merge {
    input: reports = genomad.taxonomy_reports
}
```

**Rationale:**
- Matches existing classify_multi.wdl pattern (proven on Terra/DNAnexus)
- Independent scatter blocks enable partial failure tolerance
- Database reuse across scatter shards reduces I/O costs

### Platform-Specific Considerations

**If Terra:**
- Use Google Batch backend (Cromwell 88+)
- Use gcloud storage for timing analysis (see AGENTS.md Terra section)
- Database can be in GCS bucket (gs://path/to/genomad_db)

**If DNAnexus:**
- Use dxCompiler for WDL → native workflow conversion
- Test with WDL 1.0 only (WDL 2.0 support not production-ready)
- Database uploaded to DNAnexus project storage

**If Local/HPC:**
- Use miniWDL with Singularity on HPC (if Docker unavailable)
- Database on shared filesystem (NFS/Lustre)
- Configure miniWDL for HPC scheduler (SLURM/PBS) via plugins

## Version Compatibility

| Package A | Compatible With | Notes |
|-----------|-----------------|-------|
| genomad 1.8+ | mmseqs2 15.6f452+ | Minimum required version, included in viral-classify images |
| genomad 1.8+ | Python 3.9+ | Minimum Python version, tensorflow 2.16+ requires Python 3.9+ |
| genomad 1.8+ | tensorflow 2.16+, keras 3+ | TF 2.16+ required for keras 3 compatibility |
| viral-ngs:3.0.4-classify | genomad 1.8+ | Confirmed genomad included in classify flavor (user stated) |
| miniWDL 1.13.0 | WDL 1.0, 1.1, draft-2 | Full support for all versions |
| Cromwell 88+ | WDL 1.0, 1.1 (partial) | WDL 1.1 support in progress (development-1.1), use WDL 1.0 for production |
| Cromwell 91+ | Google Batch only (GCP) | Cloud Life Sciences API removed, no longer available |

## Docker Image Strategy

### Use Existing viral-ngs Images

**Recommended approach:** Use `viral-ngs:3.0.4-classify` (or later) which includes genomad.

**Why:**
- Genomad already installed (user confirmed)
- Consistent with existing viral-pipelines Docker strategy
- Version pinned in requirements-modules.txt (line: `broadinstitute/viral-ngs=3.0.4`)
- Mirrors available on ghcr.io and quay.io for registry reliability

**Image specification in WDL tasks:**
```wdl
String docker = "quay.io/broadinstitute/viral-classify:2.5.21.0"
# OR using requirements-modules.txt pattern:
String docker = "quay.io/broadinstitute/viral-ngs:3.0.4-classify"
```

### Docker Best Practices for WDL (2026)

1. **Pin to SHA256 digest for reproducibility:**
   ```wdl
   docker: 'quay.io/broadinstitute/viral-ngs@sha256:7a47ccc3bbe8a451...'
   ```
   - Ensures exact reproducibility across platforms
   - Prevents silent updates breaking workflows

2. **Specify docker in runtime block (WDL 1.0) or requirements block (WDL 1.2):**
   - WDL 1.0: `runtime { docker: docker_image }`
   - WDL 1.2: `requirements { container: docker_image }`
   - viral-pipelines uses WDL 1.0, stick with runtime block

3. **Make images publicly accessible:**
   - Required for Terra, DNAnexus, and most platforms
   - viral-ngs images already public on quay.io/ghcr.io

4. **Support Singularity for HPC compatibility:**
   - WDL engines auto-convert docker: specification to Singularity
   - No code changes needed, works transparently

### Database Management

**User-provided database as File input:**
```wdl
input {
    File genomad_db  # User provides path to ~5GB database
}
```

**Rationale:**
- Genomad database is ~5GB (too large for Docker image)
- Database updates independently of code
- Follows kraken2/krakenuniq pattern in existing viral-pipelines (see tasks_metagenomics.wdl)
- Platform-agnostic: works with GCS (gs://), S3 (s3://), local paths, DNAnexus file IDs

## Sources

### Core Tool Documentation (HIGH Confidence)
- [genomad GitHub repository](https://github.com/apcamargo/genomad) - Official source code and documentation
- [genomad official portal](https://portal.nersc.gov/genomad/) - Installation and usage guide
- [genomad Nature Biotechnology paper](https://www.nature.com/articles/s41587-023-01953-y) - Peer-reviewed publication (September 2023)
- [genomad Bioconda package](https://bioconda.github.io/recipes/genomad/README.html) - Dependency specifications

### WDL Specifications (HIGH Confidence)
- [OpenWDL specification repository](https://github.com/openwdl/wdl) - Official WDL language specs
- [WDL 1.0 specification](https://github.com/openwdl/wdl/blob/legacy/versions/1.0/SPEC.md) - Current production standard
- [WDL 1.2 release announcement](https://openwdl.org/wdl/bioinformatics/workflows/announcing-wdl-1-2-0/) - New features and changes
- [WDL versions documentation](https://docs.openwdl.org/language-guide/versions.html) - Version differences

### Execution Engines (HIGH Confidence)
- [miniWDL releases](https://github.com/chanzuckerberg/miniwdl/releases) - Version history and release notes
- [miniWDL PyPI](https://pypi.org/project/miniwdl/) - Official package repository
- [Cromwell releases](https://github.com/broadinstitute/cromwell/releases) - Version history including v88, v91
- [Cromwell documentation](https://cromwell.readthedocs.io/) - Official execution engine docs

### Platform Integration (MEDIUM Confidence)
- [Terra WDL support](https://support.terra.bio/hc/en-us/sections/360007274612-WDL-Resources) - Terra-specific guidance
- [DNAnexus dxCompiler](https://github.com/dnanexus/dxCompiler) - WDL compilation for DNAnexus
- [DNAnexus workflow importing](https://documentation.dnanexus.com/developer/workflows/importing-workflows) - WDL import process
- [Dockstore Terra integration](https://docs.dockstore.org/en/latest/launch-with/terra-launch-with.html) - Platform interoperability

### Docker Best Practices (MEDIUM Confidence)
- [WDL Docker best practices](https://sciwiki.fredhutch.org/datascience/wdl_workflows/) - Fred Hutch guidance (2026)
- [Cromwell container tutorial](https://cromwell.readthedocs.io/en/stable/tutorials/Containers/) - Container management
- [Docker best practices 2026](https://thinksys.com/devops/docker-best-practices/) - Current industry standards
- [Docker security 2026](https://thelinuxcode.com/docker-security-best-practices-2026-hardening-the-host-images-and-runtime-without-slowing-teams-down/) - Security considerations

### Scatter-Gather Patterns (HIGH Confidence)
- [WDL scatter-gather documentation](https://docs.openwdl.org/design-patterns/scatter-gather/index.html) - Official design pattern guide
- [Terra scatter-gather tutorial](https://support.terra.bio/hc/en-us/articles/360037128572-Scatter-gather-parallelism) - Platform-specific implementation
- [WDL parallelization guide](https://hutchdatascience.org/Developing_WDL_Workflows/optimization.html) - Optimization strategies

### Viral Classification Tools (MEDIUM Confidence)
- [Broad Institute viral-classify](https://github.com/broadinstitute/viral-classify) - Classification module repository
- [viral-ngs Docker Hub](https://hub.docker.com/r/broadinstitute/viral-ngs) - Public Docker registry
- [StaPH-B genomad container](https://hub.docker.com/r/staphb/genomad) - Alternative container source
- [MVP viromics pipeline paper](https://journals.asm.org/doi/10.1128/msystems.00888-24) - Modern viromics best practices using genomad

---
*Stack research for: genomad WDL pipeline integration*
*Researched: 2026-02-12*
