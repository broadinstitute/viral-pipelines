# Pitfalls Research

**Domain:** genomad WDL pipeline integration for viral-pipelines
**Researched:** 2026-02-12
**Confidence:** MEDIUM-HIGH

## Critical Pitfalls

### Pitfall 1: Database Localization Without Disk Accounting

**What goes wrong:**
Large genomad database (5GB compressed, ~15GB extracted) gets localized to VM disk without proper accounting. Task fails with "No space left on device" during database extraction, or worse, succeeds on small test inputs but fails unpredictably on production assemblies.

**Why it happens:**
Developers calculate disk space based only on input FASTA size, forgetting that:
1. Database must be extracted from tar.gz to local disk before use
2. Extraction creates temporary files (compressed + uncompressed simultaneously)
3. genomad creates multiple output files in working directory
4. cromwell/miniwdl localize inputs to same disk as working directory

**How to avoid:**
```wdl
# Calculate disk with explicit database accounting
Float genomad_db_size = size(genomad_db, "GB")
Float input_fasta_size = size(assembly_fasta, "GB")
# DB extraction: 2x compressed size (tar.gz + extracted),
# plus 3x input (FASTA + outputs), plus 100GB padding
Int disk_size = ceil((2 * genomad_db_size + 3 * input_fasta_size + 100) / 750.0) * 750
```

**Warning signs:**
- Disk size calculation doesn't reference database file size
- Using hardcoded disk values (e.g., `disk: "375 GB"`)
- Test inputs work but production inputs fail with ENOSPC
- `df -h` in failed task shows 100% disk usage

**Phase to address:**
Phase 1: Core WDL task implementation

---

### Pitfall 2: Empty FASTA Input Causes Silent Failure

**What goes wrong:**
genomad fails when given empty FASTA file (0 sequences) or FASTA with only short sequences (<200bp that get filtered). Task exits with cryptic mmseqs2 error or produces empty outputs that break downstream tasks expecting non-empty files.

**Why it happens:**
Assembly workflows can produce empty FASTAs legitimately (e.g., taxon filter removes all contigs, SPAdes fails to assemble). genomad doesn't handle empty inputs gracefully - it expects at least one valid sequence. WDL tasks blindly pass empty files to genomad without validation.

**How to avoid:**
```wdl
# In command block, validate input before running genomad
# Count sequences, skip genomad if zero, create empty placeholder outputs
num_sequences=$(grep -c "^>" "~{assembly_fasta}" || echo "0")
if [ "$num_sequences" -eq 0 ]; then
  echo "No sequences in input FASTA, skipping genomad" >&2
  # Create empty output files matching expected structure
  touch "~{sample_name}_summary/~{sample_name}_virus_summary.tsv"
  touch "~{sample_name}_summary/~{sample_name}_plasmid_summary.tsv"
  mkdir -p "~{sample_name}_annotate"
else
  # Run genomad normally
  genomad end-to-end ...
fi
```

**Warning signs:**
- No input validation before genomad invocation
- Output declared as `File` instead of `File?` (forces failure if missing)
- Error messages mention "mmseqs2: returned non-zero exit status 1"
- No test case with empty input FASTA

**Phase to address:**
Phase 1: Core WDL task implementation

---

### Pitfall 3: Database Path Confusion (Localized vs User-Provided)

**What goes wrong:**
Task expects `File genomad_db` (cloud path to tar.gz) but gets passed local directory path from workflow. Or vice versa - task expects pre-extracted directory but gets tar.gz. Manifests as "database not found" or "not a valid genomad database" errors.

**Why it happens:**
genomad can accept database as:
1. Cloud storage tar.gz that WDL localizes and task extracts
2. Pre-extracted local directory (for miniWDL local execution)
3. Path to download location (for `genomad download-database`)

Different execution platforms (Terra, miniWDL, DNAnexus) have different conventions. Developers write task for one platform, break others.

**How to avoid:**
```wdl
input {
  File genomad_db_archive  # Always accept File (cloud path or local)
  # NOT: String genomad_db_path
}

command <<<
  # Always extract from archive, don't assume pre-extracted
  DB_DIR=$(mktemp -d --suffix _genomad_db)

  # Support multiple compression formats like existing viral-pipelines tasks
  if [[ "~{genomad_db_archive}" == *.tar.gz ]]; then
    tar -xzf "~{genomad_db_archive}" -C $DB_DIR
  elif [[ "~{genomad_db_archive}" == *.tar.zst ]]; then
    tar -I zstd -xf "~{genomad_db_archive}" -C $DB_DIR
  else
    echo "Unsupported database format" >&2
    exit 1
  fi

  # genomad expects direct path to db directory
  genomad end-to-end "~{assembly_fasta}" output_dir $DB_DIR
>>>
```

**Warning signs:**
- Input declared as `String` instead of `File` for database
- No database extraction in command block
- Test works locally (miniWDL) but fails on Terra
- Error: "genomad_db is not a valid geNomad database"

**Phase to address:**
Phase 1: Core WDL task implementation

---

### Pitfall 4: Output Directory Structure Assumptions Break Globbing

**What goes wrong:**
genomad creates complex output directory structure: `{sample}_summary/`, `{sample}_annotate/`, `{sample}_find_proviruses/`, etc. WDL task uses hardcoded paths or glob patterns that break when genomad changes output structure, or when input basename has unexpected characters (dots, dashes).

**Why it happens:**
Developers test with simple sample names (`sample1`) but production has complex names (`S20.l1.xxxx.2024-01-15`). Glob patterns like `*.tsv` match too broadly, or hardcoded paths like `output_dir/summary/virus.tsv` fail because genomad prefixes with sample name.

**How to avoid:**
```wdl
# Use genomad's --prefix flag to control output naming
String sample_basename = basename(assembly_fasta, ".fasta")
command <<<
  genomad end-to-end \
    "~{assembly_fasta}" \
    output_dir \
    $DB_DIR \
    --cleanup \
    --splits 8

  # genomad creates: output_dir/{basename}_summary/, etc.
  # Use explicit paths from known structure
  OUT_BASE="output_dir/$(basename '~{assembly_fasta}' .fasta)"

  # Copy to predictable names for WDL output
  cp "${OUT_BASE}_summary/${OUT_BASE}_virus_summary.tsv" virus_summary.tsv
  cp "${OUT_BASE}_summary/${OUT_BASE}_plasmid_summary.tsv" plasmid_summary.tsv
>>>

output {
  File virus_summary = "virus_summary.tsv"
  File plasmid_summary = "plasmid_summary.tsv"
  # NOT: glob("*_summary/*_virus_summary.tsv")[0]  # Fragile!
}
```

**Warning signs:**
- Using `glob()[0]` to extract single file from array
- No test with complex sample names containing dots/dashes
- Output paths constructed with string concatenation instead of basename
- Task works in test but fails with real sample IDs

**Phase to address:**
Phase 1: Core WDL task implementation

---

### Pitfall 5: miniwdl vs womtool Validation Divergence

**What goes wrong:**
WDL passes `miniwdl check` but fails `womtool validate`, or vice versa. CI catches issues that local validation missed. Task uses features that work in one executor but not others.

**Why it happens:**
miniwdl and womtool have different strictness levels. miniwdl enforces stricter WDL 1.0 compliance (especially optional type handling), womtool is more permissive. viral-pipelines CI runs both validators, so tasks must satisfy both.

**How to avoid:**
```wdl
# BAD: womtool allows, miniwdl rejects
runtime {
  docker: "quay.io/broadinstitute/viral-classify:" + docker_tag  # String? concatenation
}

# GOOD: Use select_first for optional types
String docker_tag_actual = select_first([docker_tag, "2.1.33.0"])
runtime {
  docker: "quay.io/broadinstitute/viral-classify:~{docker_tag_actual}"
}
```

**Warning signs:**
- Only running one validator locally
- Optional types (`String?`) used in runtime attribute expressions
- Test passes locally but CI validation fails
- Different validation results between `miniwdl check` and `womtool validate`

**Phase to address:**
Phase 1: Core WDL task implementation (validation step)

---

### Pitfall 6: Resource Requirements Mismatch (CPU, Memory for Neural Network)

**What goes wrong:**
genomad neural network classification step hangs indefinitely or fails with OOM. Task allocated for assembly workload (low memory) but genomad NN module needs substantial RAM. CPU under-allocation causes 10x slowdown.

**Why it happens:**
genomad has two classification modes:
1. Marker-based: Low memory, fast, parallelizable with `--splits`
2. Neural network: High memory (~16GB), GPU-optional, single-threaded

Developers test with small inputs where NN is fast, miss that large metagenomic assemblies (10K+ contigs) cause NN to consume excessive RAM and CPU.

**How to avoid:**
```wdl
input {
  Boolean disable_nn_classification = false  # Allow users to disable NN
  Int? machine_mem_gb
  Int? cpu
}

# Base memory on input size + NN overhead
Float input_size_gb = size(assembly_fasta, "GB")
# NN needs ~16GB base + 2GB per GB of input
Int mem_for_nn = if disable_nn_classification then 0 else ceil(16 + 2 * input_size_gb)
Int machine_mem_gb_actual = select_first([
  machine_mem_gb,
  max(32, mem_for_nn)  # At least 32GB if NN enabled
])

command <<<
  genomad end-to-end \
    "~{assembly_fasta}" output_dir $DB_DIR \
    ~{if disable_nn_classification then "--disable-nn-classification" else ""} \
    --splits 8
>>>
```

**Warning signs:**
- Hardcoded memory value (e.g., `memory: "8 GB"`) not based on input size
- No `--disable-nn-classification` option exposed
- Task "hangs" at "Classifying sequences" step
- Memory not scaled for NN classification overhead

**Phase to address:**
Phase 1: Core WDL task implementation

---

## Technical Debt Patterns

| Shortcut | Immediate Benefit | Long-term Cost | When Acceptable |
|----------|-------------------|----------------|-----------------|
| Hardcoded disk size (e.g., 375GB) | Quick to implement | Fails on large inputs, wastes money on small inputs | Never - always calculate dynamically |
| Skip empty FASTA validation | Fewer lines of code | Silent failures, broken pipelines | Never - validation is 3 lines |
| Accept `String` database path instead of `File` | Works for local testing | Breaks on cloud platforms, no localization | Never - use `File` for portability |
| Use `glob()[0]` for single outputs | Simpler than explicit paths | Breaks on naming changes, no failure if empty | Never - use explicit paths or proper optional handling |
| Disable womtool validation to pass CI | Faster local dev cycle | Breaks on Cromwell/Terra | Never - fix to satisfy both validators |
| Expose only `--disable-nn-classification` flag | Simpler interface | Users can't tune splits, other params | Phase 1 MVP only, expand in Phase 2 |
| Single Docker image version (no parameterization) | Less complexity | Hard to test version upgrades | Phase 1 MVP only, parameterize in Phase 2 |

## Integration Gotchas

| Integration | Common Mistake | Correct Approach |
|-------------|----------------|------------------|
| Terra/Cromwell | Assume database stays localized across tasks | Each task must re-localize database (File input), or pass extracted directory as String between tasks in same VM |
| miniWDL local | Use relative paths for database | Always use File inputs with absolute paths; miniWDL localizes to temp directories |
| DNAnexus | Assume local-disk runtime attribute | Use both `disks: "local-disk X LOCAL"` (GCP) and `disk: "X GB"` (TES/DNAnexus) |
| Terra data tables | Pass database as workspace data | Use workflow input, not workspace data - database is large and shouldn't be per-sample |
| Dockstore import | Hardcode database path in test JSON | Use Terra workspace bucket for database, document in README |

## Performance Traps

| Trap | Symptoms | Prevention | When It Breaks |
|------|----------|------------|----------------|
| No --splits parameter | Single-threaded execution, 10x slower | Always use `--splits 8` for marker annotation | Any assembly >1000 contigs |
| Database not decompressed in parallel | 2-5 min wait before genomad starts | Extract database in background (`&`) while preparing inputs | Database >2GB compressed |
| Neural network on large inputs | Task hangs for hours | Expose `--disable-nn-classification` flag, document tradeoff | Assemblies >10K contigs |
| Excessive disk allocation | High GCP costs | Use dynamic calculation with 750GB rounding for SSD optimization | Any production usage |
| Not using --cleanup flag | Disk fills with temp files | Always use `--cleanup` to remove intermediate files | Assemblies >1GB with small disk |

## Security Mistakes

| Mistake | Risk | Prevention |
|---------|------|------------|
| Accepting database as String path | Path traversal if user-controlled | Always use File type, let WDL handle localization |
| No database integrity check | Corrupted/malicious database could execute code | Check for expected files (genomad_db/genomad_db.dbtype), validate structure |
| Running genomad as root in container | Container escape vulnerabilities | viral-pipelines uses non-root user in Docker images (inherited) |
| Exposing arbitrary genomad flags | Command injection via shell metacharacters | Only expose vetted flags, validate enum values |

## UX Pitfalls

| Pitfall | User Impact | Better Approach |
|---------|-------------|-----------------|
| Task fails with "mmseqs2 error" | Cryptic message, no actionable info | Catch empty inputs, provide clear error: "Input FASTA is empty" |
| No progress indication for long tasks | User thinks task hung | Log progress: "Classifying X sequences", "Annotating proteins" |
| Require pre-extracted database | User confusion about format | Accept tar.gz, extract automatically |
| Output only TSV summaries | Users want FASTA sequences too | Also output `*_virus.fna` and `*_plasmid.fna` |
| No summary statistics | Users must parse TSV manually | Add simple metrics: `virus_count`, `plasmid_count` in output |

## "Looks Done But Isn't" Checklist

- [ ] **Empty input handling:** Task validates input FASTA has >0 sequences before running genomad
- [ ] **Both validators pass:** `miniwdl check` AND `womtool validate` both succeed
- [ ] **Dynamic disk calculation:** Disk size accounts for database extraction (2x compressed size) + inputs + outputs
- [ ] **Database extraction tested:** Test case with tar.gz database, not pre-extracted directory
- [ ] **Complex sample names:** Test with basename containing dots, dashes, underscores
- [ ] **Optional output handling:** Outputs are File? or Array[File], not File (to allow empty results)
- [ ] **Platform compatibility:** Test on miniWDL local AND submit to Terra test workspace
- [ ] **Version pinning:** Docker image version matches requirements-modules.txt
- [ ] **Integration test exists:** Test input JSON in test/input/WDL/miniwdl-local/
- [ ] **Resource scaling:** Memory/CPU scale with input size, not hardcoded
- [ ] **NN classification option:** `--disable-nn-classification` exposed as input parameter
- [ ] **Progress logging:** genomad output visible in stderr for monitoring
- [ ] **Cleanup enabled:** `--cleanup` flag used to remove temp files

## Recovery Strategies

| Pitfall | Recovery Cost | Recovery Steps |
|---------|---------------|----------------|
| Disk space failure | LOW | Re-run with correct disk calculation; no code changes needed |
| Empty FASTA failure | MEDIUM | Add validation logic, create empty outputs, update tests |
| Database path confusion | MEDIUM | Change input from String to File, update all test JSONs |
| Output glob breaking | LOW | Fix output paths, add test case with complex names |
| Validation divergence | MEDIUM | Fix optional type handling, re-validate with both tools |
| NN OOM failure | LOW | Add `--disable-nn-classification` input, re-run with flag |
| Wrong Docker version | LOW | Update runtime.docker, run version-wdl-runtimes.sh script |
| Missing integration test | MEDIUM | Create test input JSON, add to CI test list |

## Pitfall-to-Phase Mapping

| Pitfall | Prevention Phase | Verification |
|---------|------------------|--------------|
| Database localization disk | Phase 1: Core task | Test with 5GB database + 1GB assembly, check disk usage |
| Empty FASTA handling | Phase 1: Core task | Test case with 0-sequence FASTA, verify empty outputs created |
| Database path confusion | Phase 1: Core task | Test on miniWDL AND Terra, both should work |
| Output directory structure | Phase 1: Core task | Test with sample name `S20.l1.xxxx.2024-01-15` |
| Validation divergence | Phase 1: Core task | CI runs both miniwdl check and womtool validate |
| NN resource requirements | Phase 1: Core task | Test with/without `--disable-nn-classification` |
| No progress logging | Phase 1: Core task | Check stderr shows genomad output during execution |
| Missing --cleanup flag | Phase 1: Core task | Monitor disk usage, ensure temp files removed |
| Not using --splits | Phase 1: Core task | Benchmark with/without splits on 1000-contig assembly |
| Hardcoded Docker version | Phase 1: Core task | Run check-wdl-runtimes.sh, verify matches requirements-modules.txt |
| No integration test | Phase 1: Core task | Test JSON exists in test/input/WDL/miniwdl-local/ |
| Platform incompatibility | Phase 1: Core task | Test on miniWDL, Cromwell, submit to Terra workspace |
| Optional output misuse | Phase 2: Workflow | Workflow handles empty arrays from genomad task |
| Terra workspace database | Phase 2: Workflow | Document database localization in Terra README |

## Sources

### genomad-specific
- [geNomad GitHub repository](https://github.com/apcamargo/genomad) - Official tool repository
- [geNomad quickstart](https://portal.nersc.gov/genomad/quickstart.html) - Usage documentation
- [geNomad FAQ](https://portal.nersc.gov/genomad/faq.html) - Common issues and troubleshooting
- [geNomad issues](https://github.com/apcamargo/genomad/issues) - User-reported problems
- [geNomad ViWrap integration issue](https://github.com/AnantharamanLab/ViWrap/issues/43) - Pipeline integration gotchas

### WDL best practices
- [Cromwell runtime attributes](https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/) - Disk, memory, CPU configuration
- [Terra disk space optimization](https://terra.bio/reduce-computing-costs-by-tailoring-resource-allocations-in-workflows/) - Dynamic disk calculation
- [WDL optional outputs handling](https://support.terra.bio/hc/en-us/community/posts/360056207052-Handling-optional-outputs) - File? patterns
- [Broad WDL style guide](https://broadinstitute.github.io/long-read-pipelines/development_guide/wdl_style_guide/) - Code conventions
- [Terra workflow system overview](https://support.terra.bio/hc/en-us/articles/360055105051-Overview-How-the-workflow-system-works) - File localization
- [WDL directory localization issue](https://github.com/broadinstitute/cromwell/issues/5737) - Database handling gotcha
- [miniwdl vs womtool validation](https://miniwdl.readthedocs.io/en/latest/FAQ.html) - Validation differences
- [womtool validation strictness issue](https://github.com/broadinstitute/cromwell/issues/5354) - Version-specific validation

### Codebase patterns (viral-pipelines)
- `/home/unix/carze/projects/viral-pipelines/pipes/WDL/tasks/tasks_metagenomics.wdl` - krakenuniq task shows database extraction pattern
- `/home/unix/carze/projects/viral-pipelines/pipes/WDL/tasks/tasks_taxon_filter.wdl` - deplete_taxa task shows disk calculation and autoscaling
- `/home/unix/carze/projects/viral-pipelines/pipes/WDL/tasks/tasks_assembly.wdl` - align_reads task shows empty input handling
- `/home/unix/carze/projects/viral-pipelines/AGENTS.md` - Codebase conventions and testing requirements

---
*Pitfalls research for: genomad WDL pipeline integration*
*Researched: 2026-02-12*
