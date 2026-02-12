# Genomad Pipeline Integration

## What This Is

A WDL pipeline for viral discovery from metagenomic assemblies using geNomad. Takes assembled FASTA contigs as input and produces virus/plasmid classification, functional annotations, and taxonomic assignments for identified viral sequences.

## Core Value

Reliable viral discovery from metagenomic assemblies must work. If genomad can identify viruses in the input, the pipeline must capture and annotate them correctly.

## Requirements

### Validated

<!-- Existing viral-pipelines capabilities -->

- ✓ Reference-based assembly from BAM files with multiple aligner support — existing
- ✓ De novo metagenomic assembly — existing
- ✓ Taxonomic classification via Kraken2/Kaiju — existing
- ✓ WDL 1.0 workflows with task-based modular architecture — existing
- ✓ Container-based execution via Docker images — existing
- ✓ Platform-agnostic design (miniWDL, Cromwell, Terra, DNAnexus) — existing
- ✓ Single-sample and multi-sample workflow patterns (*_single.wdl, *_multi.wdl) — existing
- ✓ CI/CD validation with miniWDL and Cromwell — existing
- ✓ Test input JSON files following naming conventions — existing

### Active

<!-- Genomad pipeline requirements -->

- [ ] Reusable genomad tasks in dedicated task module
- [ ] Single-sample workflow (genomad_single.wdl) accepting FASTA input
- [ ] Multi-sample workflow (genomad_multi.wdl) with scatter pattern
- [ ] Database path as user-provided input parameter
- [ ] Thread count parameter with sanitization
- [ ] Empty FASTA graceful handling (skip execution)
- [ ] All genomad outputs: virus/plasmid classification, annotations, taxonomy
- [ ] Integration test cases following test/input/WDL naming conventions
- [ ] WDL syntax validation with miniwdl check
- [ ] Docker image version matching requirements-modules.txt

### Out of Scope

- Auto-download of genomad database — users must provide database path to avoid runtime downloads and storage unpredictability
- Assembly step integration — pipeline operates on existing assemblies only, assembly workflows remain separate
- Pre-filtering or post-filtering — genomad runs end-to-end on all input sequences, downstream filtering is user's responsibility
- Nextclade or Pangolin integration — genomad provides taxonomy, specialized viral typing tools handled by separate workflows

## Context

**Viral-ngs integration:**
The viral-ngs project (v3.0.4+) includes genomad in the classify Docker image (`viral-ngs:3.0.4-classify`). The Python wrapper at `viral_ngs.classify.genomad` provides the `Genomad.end_to_end()` method that:
- Validates database path exists
- Creates output directory
- Handles empty FASTA files gracefully
- Executes: `genomad end-to-end INPUT OUTPUT DATABASE [--threads N]`
- Used via CLI: `metagenomics.py genomad`

**Use case:**
Primary use is viral discovery from metagenomic samples. Researchers assemble contigs from metagenomic sequencing data (SPAdes, MEGAHIT, etc.) and need to identify which contigs are viral/plasmid sequences and obtain functional annotations.

**Genomad database:**
Requires ~5GB reference database downloaded separately. Users typically download once and reuse across analyses. Database location varies by environment (local, Terra workspace, HPC shared storage).

**Existing patterns to follow:**
- Task organization: See `tasks_assembly.wdl`, `tasks_taxon_filter.wdl` for module structure
- Single/multi workflow: See `classify_single.wdl` / `classify_kaiju.wdl` for scatter pattern
- Test structure: See `test/input/WDL/miniwdl-local/test_inputs-classify_single-local.json`
- Docker parameterization: All tasks use `String docker` input with default from requirements-modules.txt

## Constraints

- **Tech stack**: WDL 1.0, viral-ngs:3.0.4-classify Docker image
- **Compatibility**: Must work on miniWDL (local), Cromwell (cloud), Terra, DNAnexus
- **Validation**: All WDL files must pass `miniwdl check` and `womtool validate`
- **Testing**: Must include integration tests for both single and multi workflows
- **Docker versions**: Must match versions in requirements-modules.txt exactly
- **Git practices**: No agent attribution in commit messages, no amending pushed commits

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| User-provided database path | Genomad database is ~5GB; auto-download would cause runtime delays and unpredictable storage usage across platforms | — Pending |
| Separate tasks file (tasks_genomad.wdl) | Genomad is a distinct classification tool; isolating tasks allows reuse and matches existing pattern (tasks_taxon_filter.wdl for other classifiers) | — Pending |
| FASTA input only (no assembly) | Genomad operates on assembled sequences; assembly workflows already exist; keeping concerns separated | — Pending |
| viral-ngs classify Docker image | Genomad already integrated in viral-ngs v3.0.4+ classify flavor; no need for separate image | — Pending |

---
*Last updated: 2026-02-12 after initialization*
