# Feature Research: geNomad WDL Pipeline

**Domain:** Viral/plasmid metagenomics classification pipeline
**Researched:** 2026-02-12
**Confidence:** HIGH

## Feature Landscape

### Table Stakes (Users Expect These)

Features users expect in a genomad pipeline. Missing these = pipeline feels incomplete.

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| Single-sample workflow | Standard pattern in viral-pipelines, required for Terra/DNAnexus | LOW | Follow existing *_single.wdl pattern (nextclade_single, classify_single) |
| Multi-sample workflow | Batch processing essential for metagenomic studies | LOW | Follow existing *_multi.wdl pattern with scatter blocks |
| FASTA input from assemblies | geNomad operates on assembled contigs, not raw reads | LOW | Standard WDL File input, matches genomad CLI expectations |
| User-provided database path | Database is large (~5GB), users maintain their own copy | LOW | Similar to kraken2_db_tgz pattern in classify_single |
| Core genomad modules | annotate, find-proviruses, marker-classification, aggregated-classification, summary | MEDIUM | Maps to genomad end-to-end command or individual modules |
| Virus classification outputs | TSV summary, FASTA sequences, protein annotations | LOW | Primary use case for genomad |
| Plasmid classification outputs | TSV summary, FASTA sequences, protein annotations | LOW | Secondary use case, same output structure as virus |
| Provirus detection | Integrated viruses in host sequences, unique to genomad | MEDIUM | Uses conditional random field model, outputs provirus coordinates |
| Taxonomy assignment | Virus taxonomy following ICTV MSL39 | LOW | Part of standard genomad output, up to family level |
| Gene-level annotations | Functional annotation of 200K+ marker proteins | LOW | Included in genomad annotate module output |
| Version tracking | Docker image version in outputs | LOW | Standard pattern: viralngs_version output in all tasks |
| Cross-platform compatibility | miniWDL, Cromwell, Terra, DNAnexus | LOW | Constraint from AGENTS.md, requires careful runtime/resource specs |

### Differentiators (Competitive Advantage)

Features that set this pipeline apart from basic genomad wrapper. Not required, but valuable.

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Minimum contig length filtering | Prevents short-sequence misclassification (genomad warns <2,500bp risky) | LOW | Add --min-length parameter, default 2500bp based on genomad docs |
| Configurable hallmark thresholds | Control sensitivity for short sequences via --min-virus-hallmarks-short-seqs | LOW | Genomad provides --relaxed and --conservative presets |
| MultiQC integration | Aggregate QC across samples in *_multi workflow | MEDIUM | Pattern exists in classify_multi, would need custom genomad stats parser |
| Split execution mode | Run individual modules vs end-to-end for debugging/advanced use | MEDIUM | Genomad supports both, end-to-end for simplicity, modules for control |
| Merged taxonomy reports | Cross-sample summary of virus/plasmid taxonomy hits | MEDIUM | Similar to krona_merge_kraken2 in classify_multi, aggregate TSV outputs |
| Optional cleanup flag | Delete intermediate files to save storage (--cleanup) | LOW | Pass-through to genomad, useful for cloud costs |
| Resource auto-tuning | Adjust memory/disk based on input size | HIGH | Pattern not common in viral-pipelines, most tasks use fixed resources |
| AMR gene detection | Antimicrobial resistance in marker annotations (NCBIfam-AMRFinder) | LOW | Already in genomad output (_genes.tsv annotation_amr column) |
| Conjugation gene detection | Plasmid transfer genes (CONJscan annotations) | LOW | Already in genomad output (_genes.tsv annotation_conjscan column) |
| Krona visualization | Interactive taxonomy pie charts like classify_single | MEDIUM | Would require custom parser to convert genomad taxonomy to Krona format |

### Anti-Features (Commonly Requested, Often Problematic)

Features that seem good but create problems in WDL context.

| Feature | Why Requested | Why Problematic | Alternative |
|---------|---------------|-----------------|-------------|
| Automatic database download | Convenience, users don't want to manually download 5GB | Database version inconsistency across runs, cloud egress costs, workflow bloat | Require database as input (like kraken2_db_tgz), document download in README |
| Real-time database updates | Keep taxonomy current with latest ICTV releases | WDL immutability principle, reproducibility loss, version drift | Pin database version to workflow version, explicit upgrades only |
| Built-in assembly step | "One-stop shop" from reads to classifications | Conflates two separate concerns (assembly quality vs classification), doubles workflow complexity | Separate workflows: assemble_denovo then genomad_single (Unix philosophy) |
| Quality filtering of outputs | Auto-remove low-score predictions | Users may want different thresholds, genomad provides --relaxed/--conservative for post-filters | Expose score thresholds as optional inputs, provide all classifications |
| Reference-based decontamination | Filter out host sequences before classification | genomad already handles integrated proviruses, additional filtering adds complexity | Document that users should run deplete/filter_bam_to_taxa first if needed |
| Per-contig resource allocation | Dynamic resources based on sequence count | WDL task-level resources only, scattering overhead exceeds benefit | Fixed resources based on typical metagenomic assembly sizes (~10-100K contigs) |

## Feature Dependencies

```
[Single-sample workflow]
    └──requires──> [Core genomad modules]
                       └──requires──> [Database input]
                       └──requires──> [FASTA input]

[Multi-sample workflow]
    └──requires──> [Single-sample workflow] (as scatter target)
    └──enables──> [Merged taxonomy reports]
    └──enables──> [MultiQC integration]

[Provirus detection] ──enhances──> [Virus classification]

[Taxonomy assignment] ──requires──> [Core genomad modules]

[AMR/Conjugation detection] ──included in──> [Gene annotations]

[Krona visualization] ──requires──> [Taxonomy assignment]
                      ──conflicts with──> [ICTV format] (needs conversion)

[Resource auto-tuning] ──conflicts with──> [Cross-platform compatibility] (Terra/DNAnexus expect fixed resources)
```

### Dependency Notes

- **Single-sample workflow requires core modules:** The genomad_single.wdl must wrap genomad end-to-end or compose annotate→find-proviruses→marker-classification→nn-classification→aggregated-classification→summary modules
- **Multi-sample requires single-sample:** Standard WDL pattern is scatter over single-sample workflow (see classify_multi.wdl lines 54-88)
- **Krona conflicts with ICTV format:** Genomad outputs ICTV taxonomy (semicolon-delimited ranks), Krona expects NCBI taxids. Would need taxonomy mapping file or parser.
- **Resource auto-tuning conflicts with cross-platform:** Terra/DNAnexus use fixed instance types (dx_instance_type in runtime), dynamic sizing not supported

## MVP Definition

### Launch With (v1)

Minimum viable pipeline to validate concept and match table stakes.

- [x] genomad_single.wdl workflow - Basic single-sample virus/plasmid classification
- [x] genomad_multi.wdl workflow - Batch processing with scatter
- [x] Core task: genomad_endtoend - Wrap genomad end-to-end command with database input
- [x] Standard outputs: virus/plasmid summary TSV, FASTA, protein FASTA, provirus TSV, genes TSV
- [x] Taxonomy assignment - ICTV-formatted lineages in output TSV
- [x] Minimum contig filtering - --min-length parameter (default 2500bp)
- [x] Docker image - Use existing quay.io/biocontainers/genomad or build custom
- [x] Test cases - test_inputs-genomad_single-local.json, test_inputs-genomad_multi-local.json
- [x] Cross-platform validation - Test on miniWDL (local) and Cromwell (CI)

**Why essential:** These are the minimum features to be considered a "genomad pipeline" vs a "bash script." Matches existing viral-pipelines patterns and satisfies core genomad use cases (virus/plasmid ID from assemblies).

### Add After Validation (v1.x)

Features to add once core is working and users provide feedback.

- [ ] Merged taxonomy report - Aggregate virus/plasmid taxonomy across samples in *_multi workflow (similar to metag_summary_report in classify_multi.wdl line 126-129)
- [ ] Configurable hallmark thresholds - Expose --min-virus-hallmarks-short-seqs and --min-plasmid-hallmarks-short-seqs as optional inputs
- [ ] Conservative/relaxed presets - Add workflow-level parameter for --relaxed or --conservative post-classification filtering
- [ ] AMR/Conjugation summary - Extract and aggregate annotation_amr and annotation_conjscan hits from genes.tsv across samples
- [ ] Cleanup option - Add cleanup parameter to delete intermediate files (--cleanup flag)

**Trigger for adding:** User requests for cross-sample analysis, storage cost concerns (cleanup), or sensitivity tuning (hallmark thresholds).

### Future Consideration (v2+)

Features to defer until product-market fit is established.

- [ ] Krona visualization - Interactive taxonomy charts (requires ICTV→Krona format conversion)
- [ ] MultiQC integration - Aggregate stats across samples (requires custom genomad stats module for MultiQC)
- [ ] Split execution mode - Individual module tasks (annotate, find-proviruses, etc.) for advanced users
- [ ] Quality score filtering - Auto-filter low-confidence predictions based on score thresholds
- [ ] Integration with classify_single - Optional genomad step after SPAdes assembly in existing workflows

**Why defer:**
- Krona/MultiQC require significant parser development with uncertain user demand
- Split execution adds 6+ tasks with minimal benefit (end-to-end covers 95% of use cases)
- Quality filtering is subjective (users can filter TSV outputs themselves)
- Integration with classify_single needs user validation that genomad adds value vs kraken2

## Feature Prioritization Matrix

| Feature | User Value | Implementation Cost | Priority |
|---------|------------|---------------------|----------|
| Single-sample workflow | HIGH | LOW | P1 |
| Multi-sample workflow | HIGH | LOW | P1 |
| Core genomad modules | HIGH | MEDIUM | P1 |
| Virus/plasmid classification outputs | HIGH | LOW | P1 |
| Provirus detection | HIGH | LOW | P1 |
| Taxonomy assignment | HIGH | LOW | P1 |
| Gene annotations | HIGH | LOW | P1 |
| Database input | HIGH | LOW | P1 |
| Version tracking | MEDIUM | LOW | P1 |
| Cross-platform compatibility | HIGH | MEDIUM | P1 |
| Minimum contig filtering | MEDIUM | LOW | P1 |
| Merged taxonomy report | MEDIUM | MEDIUM | P2 |
| Configurable hallmark thresholds | LOW | LOW | P2 |
| Conservative/relaxed presets | LOW | LOW | P2 |
| AMR/Conjugation summary | MEDIUM | LOW | P2 |
| Cleanup option | MEDIUM | LOW | P2 |
| Krona visualization | MEDIUM | HIGH | P3 |
| MultiQC integration | MEDIUM | HIGH | P3 |
| Split execution mode | LOW | HIGH | P3 |
| Quality score filtering | LOW | MEDIUM | P3 |
| Resource auto-tuning | LOW | HIGH | P3 |

**Priority key:**
- P1: Must have for launch (table stakes + essential differentiators)
- P2: Should have, add when possible (nice-to-have differentiators)
- P3: Nice to have, future consideration (high complexity/low value)

## Competitor Feature Analysis

| Feature | nf-core/phageannotator | Standalone genomad CLI | Our WDL Approach |
|---------|------------------------|------------------------|------------------|
| Virus identification | Yes (genomad/endtoend module) | Yes (end-to-end command) | Yes (genomad_endtoend task) |
| Plasmid identification | Yes | Yes | Yes |
| Provirus detection | Yes | Yes | Yes |
| Taxonomy assignment | Yes (ICTV MSL39) | Yes | Yes |
| Batch processing | Yes (Nextflow scatter) | Manual (bash loop) | Yes (WDL scatter in *_multi) |
| Database management | genomad/download module | genomad download-database | User-provided (like kraken2) |
| Module-level execution | No (end-to-end only) | Yes (7 modules) | No (end-to-end only, defer to P3) |
| Minimum length filtering | Yes (--min-length parameter) | Yes | Yes (default 2500bp) |
| Hallmark thresholds | Yes (--min-*-hallmarks-short-seqs) | Yes | Yes (P2 feature) |
| Preset filters | Yes (--relaxed, --conservative) | Yes | Yes (P2 feature) |
| Cleanup option | Yes (--cleanup) | Yes | Yes (P2 feature) |
| Output aggregation | Yes (per-sample outputs) | No (per-sample only) | Yes (merged reports in *_multi, P2) |
| Cross-platform | Nextflow executors (Terra, AWS, etc.) | Local/cluster only | miniWDL, Cromwell, Terra, DNAnexus |
| Container strategy | BioContainers (quay.io/biocontainers/genomad) | Conda or Docker | BioContainers or custom (requirements-modules.txt) |
| Testing framework | nf-test | Manual | miniWDL + Cromwell CI (test/input/WDL/) |

**Key observations:**
- nf-core uses genomad/endtoend module exclusively (no split execution), validates our P1 approach
- nf-core provides database download step, but we reject this as anti-feature (reproducibility)
- No competitor provides Krona/MultiQC integration, confirms P3 prioritization
- Hallmark thresholds and presets are supported but not heavily documented, suggests P2 is appropriate
- BioContainers genomad image is standard across Nextflow/Snakemake, we should use quay.io/biocontainers/genomad:VERSION

## Sources

**geNomad Overview and Capabilities:**
- [Identification of mobile genetic elements with geNomad | Nature Biotechnology](https://www.nature.com/articles/s41587-023-01953-y)
- [GitHub - apcamargo/genomad](https://github.com/apcamargo/genomad)
- [geNomad Portal](https://portal.nersc.gov/genomad/)

**geNomad Pipeline and Modules:**
- [The geNomad pipeline](https://portal.nersc.gov/genomad/pipeline.html)
- [geNomad NIH HPC Documentation](https://hpc.nih.gov/apps/genomad.html)
- [geNomad Quickstart](https://gensoft.pasteur.fr/docs/genomad/1.8.1/quickstart.html)

**geNomad Outputs and Features:**
- [Post-classification filtering](https://portal.nersc.gov/genomad/post_classification_filtering.html)
- [geNomad FAQ](https://portal.nersc.gov/genomad/faq.html)

**Workflow Integration Patterns:**
- [nf-core genomad/endtoend module](https://nf-co.re/modules/genomad_endtoend)
- [nf-core/phageannotator pipeline](https://nf-co.re/phageannotator/dev/)
- [EBI Metagenomics mobilome-annotation-pipeline](https://github.com/EBI-Metagenomics/mobilome-annotation-pipeline)

**Best Practices:**
- [metaFun metagenomic analysis pipeline](https://pmc.ncbi.nlm.nih.gov/articles/PMC12818822/)
- [MetaflowX scalable metagenomic workflow](https://academic.oup.com/nar/article/53/18/gkaf954/8271006)
- [MicrobiomeBestPracticeReview](https://github.com/grimmlab/MicrobiomeBestPracticeReview)

---
*Feature research for: geNomad WDL pipeline for viral-pipelines*
*Researched: 2026-02-12*
