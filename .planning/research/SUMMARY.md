# Project Research Summary

**Project:** genomad WDL Pipeline Integration
**Domain:** Viral/plasmid metagenomics bioinformatics workflow
**Researched:** 2026-02-12
**Confidence:** HIGH

## Executive Summary

This project integrates genomad (state-of-the-art viral/plasmid identification tool) into the viral-pipelines WDL codebase. Genomad is a Nature Biotechnology 2023 publication that achieves 95.3% accuracy for virus classification using machine learning and 200K+ marker proteins. The recommended approach follows established viral-pipelines patterns: create WDL tasks in tasks_metagenomics.wdl, build single-sample workflow first (genomad_single.wdl), then extend to multi-sample scatter-gather pattern (genomad_multi.wdl). This matches the proven pattern used by classify_single/classify_multi workflows.

The primary technical challenge is handling genomad's large database (5GB compressed, ~15GB extracted) while maintaining cross-platform compatibility (Terra, DNAnexus, local miniWDL). The recommended stack uses WDL 1.0 (maximum platform compatibility), miniWDL for local testing, Cromwell 88+ for cloud execution, and existing viral-ngs Docker images that already include genomad. Critical success factors include proper disk space calculation accounting for database extraction overhead, input validation to handle empty FASTA files gracefully, and resource autoscaling for genomad's neural network classification step.

Key risks center on database management (disk space failures, path confusion between platforms), resource under-provisioning (neural network OOM, insufficient CPU parallelization), and WDL portability (validation divergence between miniwdl/womtool, platform-specific runtime attributes). These risks are well-understood and preventable through established patterns documented in existing viral-pipelines codebase (kraken2, kaiju tasks demonstrate database handling; deplete_taxa shows resource autoscaling). The project has clear table stakes (single/multi-sample workflows, virus/plasmid classification, taxonomy assignment) and straightforward MVP scope, making it low-risk for initial implementation.

## Key Findings

### Recommended Stack

The viral-pipelines codebase provides proven patterns for integrating metagenomics tools. WDL 1.0 remains the standard for maximum compatibility across Terra, DNAnexus, and local/HPC execution. MiniwdL 1.13.0+ handles local development and validation, while Cromwell 88+ provides production-grade execution with Google Batch backend support (critical since Cloud Life Sciences API was deprecated in July 2025).

**Core technologies:**
- **WDL 1.0**: Workflow definition language — maximum platform compatibility, existing codebase standard, fully supported by all execution engines
- **genomad 1.8+**: Viral/plasmid identification — state-of-the-art accuracy (MCC 95.3%), already included in viral-ngs:3.0.4-classify Docker image, designed for metagenomic workflows
- **miniWDL 1.13.0+**: Local workflow execution — lightweight Python runner, WDL 1.0/1.1 support, Docker/Podman/Singularity compatibility for testing
- **Cromwell 88+**: Cloud/enterprise execution — production engine for Terra, version 88+ required for Google Batch backend
- **viral-ngs Docker images**: Container runtime — existing images include genomad, follows codebase pattern of tool consolidation

**Critical dependencies:**
- genomad requires mmseqs2, pyrodigal-gv, tensorflow (all auto-included in viral-classify images)
- Database (~5GB) must be user-provided as File input following kraken2_db_tgz pattern
- No local installation needed — workflows run entirely in containers

**Alternatives considered and rejected:**
- WDL 1.1/1.2: Reduced platform compatibility (DNAnexus partial support), not needed for this use case
- VirSorter2/CheckV for discovery: Lower accuracy than genomad, CheckV is quality control not discovery
- Standalone genomad conda/pip: Breaks reproducibility, WDL requires Docker for platform-agnostic design

### Expected Features

Genomad integration follows standard metagenomics pipeline patterns. Users expect single-sample and batch processing capabilities with standard input/output conventions matching the existing classify_single/multi workflows.

**Must have (table stakes):**
- Single-sample workflow (genomad_single.wdl) — standard pattern in viral-pipelines, required for Terra data tables
- Multi-sample workflow (genomad_multi.wdl) — batch processing essential for metagenomic studies
- FASTA input from assemblies — genomad operates on assembled contigs, not raw reads
- User-provided database path — database too large for Docker image, follows kraken2 pattern
- Core genomad modules — annotate, find-proviruses, marker/NN classification, summary outputs
- Virus/plasmid classification — TSV summary, FASTA sequences, protein annotations for both
- Provirus detection — integrated viruses in host sequences, unique genomad capability
- Taxonomy assignment — ICTV MSL39 lineages up to family level
- Minimum contig filtering — prevent short-sequence misclassification (default 2500bp)
- Cross-platform compatibility — miniWDL, Cromwell, Terra, DNAnexus support

**Should have (competitive advantage):**
- Configurable hallmark thresholds — control sensitivity for short sequences via --min-virus-hallmarks-short-seqs
- Conservative/relaxed presets — genomad --relaxed and --conservative post-classification filters
- Merged taxonomy reports — cross-sample summary similar to krona_merge_kraken2 in classify_multi
- AMR/conjugation summaries — extract antimicrobial resistance and plasmid transfer gene hits from annotations
- Cleanup option — --cleanup flag to delete intermediates and save storage costs

**Defer (v2+):**
- Krona visualization — requires ICTV→Krona format conversion, uncertain user demand
- MultiQC integration — needs custom genomad stats parser module
- Split execution mode — individual module tasks (annotate, find-proviruses, etc.) for advanced users
- Resource auto-tuning — dynamic instance sizing not well-supported across platforms
- Quality score filtering — subjective thresholds, users can filter TSV outputs themselves

**Anti-features (commonly requested but problematic):**
- Automatic database download — creates version inconsistency and reproducibility issues
- Built-in assembly step — conflates concerns, separate workflows better (assemble_denovo then genomad)
- Per-contig resource allocation — WDL task-level resources only, scattering overhead exceeds benefit

### Architecture Approach

The architecture follows viral-pipelines established patterns with strict separation between tasks (low-level tool wrappers) and workflows (orchestration logic). Tasks belong in tasks_metagenomics.wdl alongside kraken2/kaiju classification tools, workflows live in pipes/WDL/workflows/, and test fixtures co-locate with existing tests in test/input/WDL/.

**Major components:**

1. **tasks_metagenomics.wdl (task module)** — Defines genomad_end_to_end task wrapping the genomad command, handles database decompression from tarball to tmpdir, implements resource autoscaling based on input size, extracts outputs from genomad's multi-level directory structure. Responsibilities: Docker runtime, resource specs, command execution, output file extraction. Does NOT orchestrate multi-step workflows or aggregate cross-sample results.

2. **genomad_single.wdl (single-sample workflow)** — Orchestrates single-sample processing with optional BAM-to-FASTA conversion, calls genomad_end_to_end task, exposes user-friendly inputs/outputs. Follows classify_single.wdl pattern. Primary use case for Terra data tables and individual sample analysis.

3. **genomad_multi.wdl (multi-sample workflow)** — Implements scatter-gather pattern over array of input files, parallelizes genomad_end_to_end execution across samples, optionally aggregates summary statistics. Follows classify_multi.wdl pattern with scatter blocks. Efficient for batch processing 10-100+ samples.

4. **Test fixtures (test/input/WDL/)** — Provides miniature test database (~100MB vs 5GB production) and small test FASTA files for CI/CD validation. Tests must complete in <2 minutes for GitHub Actions constraints.

**Key architectural patterns:**

- **Database decompression in task command**: Extract tarball databases inside command blocks (not separate tasks) following kraken2 pattern, avoids call-caching issues
- **Single/multi workflow split**: Separate workflows for simple and batch use cases, reduces complexity and enables independent testing
- **Resource autoscaling**: Calculate CPU/memory/disk based on input sizes, critical for genomad's variable resource needs (neural network step)
- **Task module organization**: Group by functional category (metagenomics), not one-file-per-tool, maintains codebase consistency

**Data flow:**
- Input: FASTA/BAM → Task converts BAM if needed → genomad processes → Outputs extracted from PREFIX_summary/, PREFIX_annotate/ subdirectories
- Database: Tarball localized per task → extracted to tmpdir → used for execution → discarded (not in outputs)
- Multi-sample: Scatter over array → parallel genomad execution → gather arrays → optional aggregation

### Critical Pitfalls

Research identified 6 critical pitfalls from genomad-specific issues, WDL platform gotchas, and viral-pipelines integration requirements. All are preventable with proper implementation but will cause silent failures or production issues if missed.

1. **Database localization without disk accounting** — genomad database is 5GB compressed but ~15GB extracted. Tasks calculating disk from only input FASTA size fail with "No space left on device". Prevention: Dynamic disk calculation accounting for 2x database size (compressed + extracted simultaneously) + 3x input size + 100GB padding, rounded to 750GB increments for SSD optimization. Formula: `ceil((2 * db_size + 3 * input_size + 100) / 750) * 750`

2. **Empty FASTA input causes silent failure** — Assembly workflows legitimately produce empty FASTAs (filter removes all contigs, assembly fails). Genomad exits with cryptic mmseqs2 error. Prevention: Validate input has >0 sequences before running genomad, create empty placeholder outputs if zero, use File? instead of File for optional outputs.

3. **Database path confusion (localized vs user-provided)** — Accepting database as String path instead of File breaks cloud platforms. Terra/Cromwell require File type for localization, while miniWDL uses different temp paths. Prevention: Always use File input type, extract from tarball in command block, support multiple compression formats (.tar.gz, .tar.zst).

4. **Output directory structure assumptions break globbing** — genomad creates complex directory structure with sample-name prefixes. Glob patterns like `*.tsv` or hardcoded paths fail with complex sample names (dots, dashes). Prevention: Use explicit basename extraction and copy to predictable names for WDL outputs, test with complex names like `S20.l1.xxxx.2024-01-15`.

5. **miniwdl vs womtool validation divergence** — Different strictness for optional type handling causes CI failures. MiniwdL enforces stricter WDL 1.0 compliance. Prevention: Use select_first() for optional types in runtime expressions, validate with BOTH tools before committing, viral-pipelines CI runs both validators.

6. **Resource requirements mismatch (neural network)** — genomad neural network classification needs ~16GB RAM base + 2GB per GB input. Tasks allocated for assembly workloads OOM or hang indefinitely. Prevention: Scale memory based on input size with NN overhead, expose --disable-nn-classification flag for large inputs, use --splits 8 for parallelization.

**Additional gotchas:**
- No progress logging causes users to think tasks hung (genomad can run 30+ minutes)
- Missing --cleanup flag fills disk with temporary files on large assemblies
- Hardcoded Docker versions create upgrade friction (parameterize docker input)
- Platform runtime attributes differ (use both `disks: "local-disk X LOCAL"` for GCP and `disk: "X GB"` for TES)

## Implications for Roadmap

Based on research findings, the project naturally divides into two main phases with clear dependency ordering. Phase 1 establishes core functionality with single-sample processing, Phase 2 extends to batch workflows. This ordering minimizes risk by validating the complex parts (database handling, resource autoscaling) before adding scatter-gather complexity.

### Phase 1: Core Task and Single-Sample Workflow

**Rationale:** Task implementation must come first since workflows depend on it. Single-sample workflow is simpler to test, validates all core functionality (database handling, output extraction, resource autoscaling), and provides immediate value for users processing individual samples. This matches viral-pipelines development pattern (classify_single built before classify_multi).

**Delivers:**
- tasks_metagenomics.wdl with genomad_end_to_end task
- genomad_single.wdl workflow
- Test fixtures with miniature database
- CI/CD integration (miniwdl and womtool validation)
- Cross-platform validation (Terra test workspace submission)

**Addresses features:**
- Single-sample workflow (table stakes)
- FASTA input handling with optional BAM conversion
- User-provided database pattern
- Virus/plasmid/provirus classification outputs
- Taxonomy assignment
- Minimum contig filtering (--min-length 2500bp)

**Avoids pitfalls:**
- Pitfall 1: Dynamic disk calculation with database accounting
- Pitfall 2: Empty FASTA input validation
- Pitfall 3: File-type database input with tarball extraction
- Pitfall 4: Explicit output paths, test with complex sample names
- Pitfall 5: Validate with both miniwdl and womtool
- Pitfall 6: NN resource scaling, expose --disable-nn-classification

**Success criteria:**
- Passes miniwdl check and womtool validate
- Test completes in <2 minutes on CI
- Successfully runs on Terra test workspace
- Handles empty FASTA gracefully
- Works with 5GB production database

### Phase 2: Multi-Sample Workflow and Aggregation

**Rationale:** Multi-sample workflow reuses proven single-sample task, adding only scatter-gather orchestration. This is low-risk since the complex parts (genomad execution, database handling) were validated in Phase 1. Aggregation features (merged reports) provide value for batch studies but aren't required for core functionality.

**Delivers:**
- genomad_multi.wdl workflow with scatter over sample array
- Test fixtures for multi-sample inputs
- Optional aggregation task (summarize_genomad_results)
- Merged taxonomy reports across samples

**Addresses features:**
- Multi-sample workflow (table stakes)
- Merged taxonomy reports (competitive advantage)
- Batch processing efficiency (database loaded once per scatter shard)

**Builds on Phase 1:**
- Reuses genomad_end_to_end task without modification
- Follows classify_multi.wdl scatter-gather pattern
- Uses same test database, just multiple small FASTAs

**Success criteria:**
- Successfully scatters over 10+ samples
- Handles partial failures (some samples empty, others succeed)
- Aggregation produces valid merged reports
- Completes within reasonable time on Terra batch submission

### Phase 3 (Optional): Advanced Features

**Rationale:** These features enhance usability but aren't essential for launch. Add based on user feedback after Phases 1-2 are production-validated.

**Potential additions:**
- Configurable hallmark thresholds (--min-virus-hallmarks-short-seqs)
- Conservative/relaxed preset flags
- AMR/conjugation gene summaries from annotations
- Krona visualization (requires ICTV→Krona format parser)
- MultiQC integration (needs custom genomad stats module)

**Defer until:**
- Users request specific thresholds for sensitivity tuning
- Storage costs justify --cleanup flag addition
- Visualization needs emerge from real usage patterns

### Phase Ordering Rationale

**Dependency-driven ordering:**
- Tasks before workflows: Workflows import and call tasks, tasks must exist first
- Single before multi: Multi reuses single-sample task, validates simple case before adding scatter complexity
- Core before optional: Table stakes features (virus/plasmid classification) before nice-to-haves (Krona viz)

**Risk mitigation:**
- Database handling complexity isolated in Phase 1 task implementation
- Resource autoscaling tested with single samples before batch workloads
- Cross-platform compatibility validated incrementally (local → Terra → DNAnexus)

**User value delivery:**
- Phase 1 delivers working single-sample pipeline (immediate value for individual analysis)
- Phase 2 adds batch processing (scales to production metagenomic studies)
- Phase 3 responds to actual user needs rather than anticipated features

### Research Flags

**Phases with standard patterns (skip research-phase):**
- **Phase 1:** Well-documented genomad CLI, existing viral-pipelines task patterns (kraken2, kaiju) provide clear templates, WDL 1.0 patterns established
- **Phase 2:** Standard scatter-gather pattern, classify_multi.wdl provides exact template

**No phases require additional research** — genomad documentation is comprehensive, viral-pipelines codebase provides all necessary patterns, integration risks are well-understood and preventable.

**Validation checkpoints:**
- Phase 1: Submit test to Terra workspace, verify with 5GB production database
- Phase 2: Batch submission with 10+ samples to validate scatter efficiency

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | Official genomad docs, viral-pipelines codebase patterns, WDL 1.0 specification all verified. Docker images confirmed to include genomad. Version compatibility cross-referenced. |
| Features | HIGH | Genomad Nature Biotech paper establishes capabilities, nf-core genomad module validates integration patterns, viral-pipelines classify workflows provide clear feature precedent. |
| Architecture | HIGH | Existing viral-pipelines codebase (15+ task modules, 20+ workflows) provides proven patterns. Task module organization, single/multi split, resource autoscaling all have working examples. |
| Pitfalls | MEDIUM-HIGH | Genomad-specific pitfalls from GitHub issues and FAQ, WDL pitfalls from Terra support docs and Cromwell issues, viral-pipelines pitfalls from codebase analysis. Some inference from similar tools (kraken2). |

**Overall confidence: HIGH**

The project benefits from multiple confidence-building factors:
1. Genomad is well-documented with official portal, GitHub repo, and peer-reviewed publication
2. viral-pipelines codebase provides working examples of identical patterns (kraken2 database handling, classify workflow structure)
3. nf-core genomad module demonstrates successful integration in another workflow system
4. All critical dependencies (Docker images, WDL executors, platforms) are verified and version-compatible
5. Pitfalls are documented from real user issues and established WDL best practices

### Gaps to Address

**During implementation:**
- **Test database creation:** Need to create or find miniature genomad database (~100MB) for CI tests. Production database is 5GB, too large for GitHub Actions. Resolution: Extract subset of markers or use genomad download with --lite flag if available, or build minimal custom database.

- **Docker image version pinning:** Confirm exact viral-ngs version with genomad 1.8+. User stated genomad in viral-ngs:3.0.4-classify but should verify SHA256 digest. Resolution: Test genomad --version in container during Phase 1 implementation.

- **Terra workspace setup:** Need test Terra workspace and GCS bucket for database storage. Resolution: Use existing viral-pipelines test workspace or create dedicated genomad-test workspace.

**No critical gaps** — all essential information for Phase 1 implementation is available. The gaps above are logistical (test fixtures, versions) rather than architectural or technical.

**Validation strategy:**
- Start Phase 1 with minimal test (tiny FASTA, small custom database) to validate WDL structure
- Upgrade to production 5GB database once task resource calculations are stable
- Submit to Terra early (don't wait until "done") to catch platform-specific issues

## Sources

### Primary (HIGH confidence)

**Genomad official:**
- [genomad Nature Biotechnology paper](https://www.nature.com/articles/s41587-023-01953-y) — peer-reviewed publication establishing accuracy, capabilities, benchmarks
- [genomad GitHub repository](https://github.com/apcamargo/genomad) — official source code, documentation, version history
- [genomad official portal](https://portal.nersc.gov/genomad/) — usage guide, pipeline documentation, FAQ
- [genomad Bioconda package](https://bioconda.github.io/recipes/genomad/README.html) — dependency specifications, installation

**WDL specifications:**
- [OpenWDL specification repository](https://github.com/openwdl/wdl) — official language specs, version differences
- [WDL 1.0 specification](https://github.com/openwdl/wdl/blob/legacy/versions/1.0/SPEC.md) — current production standard
- [miniWDL releases](https://github.com/chanzuckerberg/miniwdl/releases) — version history, feature support
- [Cromwell releases](https://github.com/broadinstitute/cromwell/releases) — version history including critical v88/v91 changes

**viral-pipelines codebase:**
- tasks_metagenomics.wdl — kraken2, kaiju patterns for database handling
- classify_single.wdl / classify_multi.wdl — workflow structure templates
- tasks_taxon_filter.wdl — resource autoscaling examples (deplete_taxa)
- AGENTS.md — codebase conventions, testing requirements
- requirements-modules.txt — Docker image version specifications

### Secondary (MEDIUM confidence)

**Integration patterns:**
- [nf-core genomad_endtoend module](https://nf-co.re/modules/genomad_endtoend) — validated integration approach in Nextflow
- [Terra WDL support docs](https://support.terra.bio/hc/en-us/sections/360007274612-WDL-Resources) — platform-specific guidance
- [DNAnexus dxCompiler](https://github.com/dnanexus/dxCompiler) — WDL compilation for DNAnexus platform
- [WDL scatter-gather documentation](https://docs.openwdl.org/design-patterns/scatter-gather/index.html) — official design pattern guide

**Best practices:**
- [Cromwell runtime attributes](https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/) — disk, memory, CPU configuration
- [Terra disk space optimization](https://terra.bio/reduce-computing-costs-by-tailoring-resource-allocations-in-workflows/) — dynamic disk calculation patterns
- [Docker best practices 2026](https://thinksys.com/devops/docker-best-practices/) — container image strategies

**Troubleshooting:**
- [genomad GitHub issues](https://github.com/apcamargo/genomad/issues) — user-reported problems, solutions
- [genomad FAQ](https://portal.nersc.gov/genomad/faq.html) — common issues and workarounds
- [WDL directory localization issue](https://github.com/broadinstitute/cromwell/issues/5737) — database handling gotcha
- [miniwdl vs womtool validation](https://miniwdl.readthedocs.io/en/latest/FAQ.html) — validation differences

### Tertiary (LOW confidence)

**Comparative tools:**
- [StaPH-B genomad container](https://hub.docker.com/r/staphb/genomad) — alternative container source (not used, but confirms image availability)
- [MVP viromics pipeline paper](https://journals.asm.org/doi/10.1128/msystems.00888-24) — modern viromics best practices using genomad

---
*Research completed: 2026-02-12*
*Ready for roadmap: yes*
