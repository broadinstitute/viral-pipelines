# Technology Stack: Docker Image and Python Package Strategy

**Project:** VirNucPro Contig & Read Classification WDL Tasks
**Researched:** 2026-03-31
**Overall confidence:** HIGH — all findings are from direct codebase inspection

---

## Executive Summary

This project adds two WDL tasks to `tasks_metagenomics.wdl`. Both tasks run Python scripts that require third-party packages (`pandas` for contig classification, `duckdb` for read classification). The repo has a well-established pattern for this: use `quay.io/broadinstitute/py3-bio:0.1.3` as the base image and embed the Python script inline using the `python3<<CODE ... CODE` heredoc idiom. When a package is missing from py3-bio, the convention (explicitly logged in PROJECT.md) is to `pip install` it at task runtime rather than build a new Docker image.

No task in the current codebase uses runtime `pip install` yet — this project introduces that pattern for the first time. There is no `cat > script.py` pattern anywhere; the codebase exclusively uses the inline `python3<<CODE` heredoc approach.

---

## Findings by Question

### 1. Do any existing tasks run `pip install` at runtime?

**No.** A full search across all files under `pipes/WDL/tasks/` found zero occurrences of `pip install`.

- Search performed: `grep -rn "pip install"` — no matches.
- The decision to use `pip install duckdb` at runtime is recorded in `PROJECT.md` (Key Decisions table, row: "py3-bio + pip install duckdb") but has not yet been implemented in any task file.

This means the new `classify_reads_by_contig` task will introduce the first runtime pip-install pattern in the repo. The pattern to follow should be:

```bash
pip install duckdb --quiet
python3<<CODE
import duckdb
...
CODE
```

Placed inside the WDL `command <<<` block, before the `python3<<CODE` heredoc.

---

### 2. What does a py3-bio task command block look like?

The canonical pattern appears in at least 7 locations. The full structure is:

**Minimal form** — `tasks_utils.wdl:551-562` (`sanitize_fasta_headers`):

```wdl
task sanitize_fasta_headers {
  input {
    File   in_fasta
    String out_filename = "~{basename(in_fasta, '.fasta')}-sanitized.fasta"
  }
  String docker = "quay.io/broadinstitute/py3-bio:0.1.3"
  Int    disk_size = 375
  command <<<
    python3<<CODE
    import re
    import Bio.SeqIO
    with open('~{in_fasta}', 'rt') as inf:
      with open('~{out_filename}', 'wt') as outf:
        for seq in Bio.SeqIO.parse(inf, 'fasta'):
          seq.id = re.sub(r'[^0-9A-Za-z!_-]', '-', seq.id)
          seq.description = seq.id
          seq.name = seq.id
          Bio.SeqIO.write(seq, outf, 'fasta')
    CODE
  >>>
  output {
    File sanitized_fasta = out_filename
  }
  runtime {
    docker: docker
    memory: "2 GB"
    cpu:    2
    disks: "local-disk ~{disk_size} LOCAL"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 1
  }
}
```

**Extended form with `set -e` guard** — `tasks_utils.wdl:918-942` (`tsv_drop_cols`):

```wdl
command <<<
    set -e
    python3<<CODE
    import pandas as pd
    in_tsv = "~{in_tsv}"
    df = pd.read_csv(in_tsv, sep='\t', dtype=str).dropna(how='all')
    drop_cols = list(x for x in '~{sep="*" drop_cols}'.split('*') if x)
    if drop_cols:
        df.drop(columns=drop_cols, inplace=True)
    df.to_csv("~{out_filename}", sep='\t', index=False)
    CODE
>>>
```

**Observations:**

- The delimiter between shell and Python is `python3<<CODE` (no space before `<<`, no quotes around `CODE`). Two variant spellings appear in the codebase: `python3<<CODE` and `python3 << CODE` — both work identically in bash.
- WDL interpolation (`~{variable}`) works inside the heredoc body because WDL processes it before the shell sees the text. This is the mechanism used to inject input file paths and parameter values into the Python script.
- `set -e` is used in multi-step tasks where the Python block is preceded by shell setup steps. Single-block tasks sometimes omit it.
- The `docker` input is declared either as a task-level variable (`String docker = "..."` outside `input {}`) or as an input parameter with a default, making the image overridable per-task-call. Both patterns are in active use.

Key files and lines:
- `tasks_utils.wdl:549-575` — `sanitize_fasta_headers` (Bio.SeqIO, minimal pattern)
- `tasks_utils.wdl:915-942` — `tsv_drop_cols` (pandas, set -e pattern)
- `tasks_assembly.wdl:644-694` — `ivar_trim_stats` (pandas + plotly, medium-length script)
- `tasks_reports.wdl:232-284` — `merge_coverage_per_position` (pandas + Bio.SeqIO, ~30 lines)
- `tasks_sarscov2.wdl:257-386` — `sars_cov_2_meta_etl` (pandas + numpy, ~100 lines inline)
- `tasks_sarscov2.wdl:407-~550` — `crsp_meta_etl` (pandas + numpy, ~100+ lines inline)
- `tasks_ncbi.wdl:1097-~1250` — `generate_author_sbt_file` (pyyaml + jinja2, multi-function inline script)
- `tasks_demux.wdl:102-141` — `revcomp_i5` (Bio.Seq, short script)

The `tasks_sarscov2.wdl` and `tasks_ncbi.wdl` examples confirm that **multi-hundred-line Python scripts embedded inline are normal in this codebase**.

---

### 3. Is there a `cat > *.py << 'EOF'` heredoc pattern anywhere?

**No.** A full search for `cat >`, `python3 - <<`, and `<< 'EOF'` / `<<EOF` patterns found no matches in the tasks directory.

The codebase uses only the `python3<<CODE ... CODE` inline execution pattern. There is no convention of writing a `.py` file to disk first and then calling it.

---

### 4. Recommended approach for embedding a multi-hundred-line Python script

Based on direct evidence from `tasks_sarscov2.wdl` (both ETL tasks, ~100-150 lines each) and `tasks_ncbi.wdl:generate_author_sbt_file` (multi-function script with helpers), the repo-established approach is:

**Embed the entire script inline using `python3<<CODE ... CODE`.**

Template for the new tasks:

```wdl
command <<<
    set -e
    pip install duckdb --quiet   # only needed where duckdb is required
    python3<<CODE
    # --- full script body here ---
    # WDL interpolations like ~{in_file} work normally
    CODE
>>>
```

Key rules inferred from existing tasks:

1. Use `set -e` at the top of any `command <<<` block that has pre-Python shell steps (like pip install).
2. The `CODE` terminator must appear at the start of the line with no leading whitespace (WDL heredoc rule).
3. Indentation inside the heredoc is cosmetic — Python sees it, so use consistent indentation but do not add leading spaces to the terminator.
4. WDL `~{...}` interpolations are resolved before the shell executes, so they work inside the Python block. Use `"~{in_file}"` for file paths, `~{true="True" false="False" bool_var}` for booleans, and `~{default='None' '"' + optional_var + '"'}` for optional strings.
5. No size limit: the `tasks_sarscov2.wdl` and `tasks_ncbi.wdl` examples show scripts in the 100-200+ line range are acceptable.

---

## Recommended Stack

### Docker Images

| Image | Registry | Purpose | Used By |
|-------|----------|---------|---------|
| `quay.io/broadinstitute/py3-bio:0.1.3` | quay.io | Python 3 + biopython + pandas + numpy + plotly + pyyaml + jinja2 + epiweeks | All Python-heavy tasks |
| `ghcr.io/broadinstitute/viral-ngs:3.0.6-core` | ghcr.io | Full viral-ngs toolchain (bash-heavy tasks) | Assembly, alignment tasks |
| `python:slim` | Docker Hub | Minimal Python (no scientific libraries) | Simple CSV manipulation |
| `ubuntu` | Docker Hub | Shell-only tasks | Text processing |
| `quay.io/broadinstitute/virnucpro-cuda:1.0.9` | quay.io | VirNucPro GPU scoring | `classify_virnucpro` task |

Source: `requirements-modules.txt` and task `runtime.docker` declarations.

### For the New Tasks (v1.0–v2.0)

| Task | Image | Additional Runtime Step |
|------|-------|------------------------|
| `classify_virnucpro_contigs` | `quay.io/broadinstitute/py3-bio:0.1.3` | None — pandas is already present |
| `classify_reads_by_contig` | `quay.io/broadinstitute/py3-bio:0.1.3` | `pip install duckdb --quiet` before `python3<<CODE` |

**Rationale:** py3-bio already ships pandas (confirmed by its use in `tsv_drop_cols`, `ivar_trim_stats`, `sars_cov_2_meta_etl`, `crsp_meta_etl`). duckdb is not present in py3-bio but is lightweight and fast to install. Avoiding a custom image build eliminates Docker registry overhead and build maintenance. This tradeoff is explicitly captured in PROJECT.md Key Decisions.

---

## v3.0 Milestone: Centrifuger Docker Image

**Researched:** 2026-04-01
**Confidence:** HIGH

### Recommended Image

```
quay.io/biocontainers/centrifuger:1.1.0--hf426362_0
```

- **Size:** 133.9 MB compressed
- **Source:** BioContainers — auto-built from Bioconda package `centrifuger==1.1.0`
- **Registry:** `quay.io/biocontainers/` (same registry as many other BioContainers images used in viral-pipelines-adjacent tools)
- **No custom build required.** The BioContainers pipeline automatically publishes a container for every merged Bioconda PR. Centrifuger 1.1.0 is available today.

### Centrifuger Version

**v1.1.0** — released February 18, 2025 by `mourisl` (Li Song). This is the current stable release.

Centrifuger is NOT the same as Centrifuge. They share a conceptual ancestor but differ in:
- Centrifuger uses FM-index + run-block compressed BWT (more memory-efficient)
- Centrifuger binary names all use the `centrifuger-*` prefix
- Centrifuger index format is incompatible with Centrifuge indexes

### Binaries in the Image

All executables use the `centrifuger-*` prefix. Do not confuse with `centrifuge-*` (the older tool's binaries).

| Binary | Role in WDL task |
|--------|-----------------|
| `centrifuger` | Primary classifier — call for classification in the WDL command block |
| `centrifuger-kreport` | Converts per-read output to Kraken-style report — call after `centrifuger` |
| `centrifuger-build` | Index construction — out of scope for v3.0 classification task |
| `centrifuger-download` | NCBI reference downloader — out of scope for v3.0 |
| `centrifuger-quant` | Abundance estimation / profiling — out of scope for v3.0 |
| `centrifuger-inspect` | Index inspection — out of scope for v3.0 |
| `centrifuger-promote` | Genome promotion — out of scope for v3.0 |

The v3.0 WDL task needs only `centrifuger` and `centrifuger-kreport`.

### Why Not Existing Broad Images

| Image | Verdict | Reason |
|-------|---------|--------|
| `quay.io/broadinstitute/viral-ngs:3.0.10-classify` | Do NOT use | Does not include centrifuger binary; wraps Kraken2 via Python metagenomics CLI |
| `quay.io/broadinstitute/viral-classify:2.1.33.0` | Do NOT use | KrakenUniq-focused; archived repo; no centrifuger |
| `quay.io/broadinstitute/py3-bio:0.1.3` | Do NOT use | Python-only image; no compiled C++ binaries |

None of the existing Broad images include centrifuger. The BioContainers image is the correct choice and requires zero custom work.

### Why Not Build a Custom Image

The PROJECT.md Key Decisions table already records: "Custom Docker image build — will use existing tool image or runtime install" as out of scope. The BioContainers image satisfies the requirement without any custom build. The precedent in this repo is to use specialist images for specialist tools (e.g. `quay.io/broadinstitute/ncbi-tools:2.12.0` for NCBI tasks, `quay.io/broadinstitute/qiime2:latest` for 16S tasks). Centrifuger follows the same pattern.

### WDL Task Docker Declaration Pattern

```wdl
String docker = "quay.io/biocontainers/centrifuger:1.1.0--hf426362_0"
```

Use a pinned tag with the full `version--buildhash_buildnum` format. BioContainers does not provide a `latest` tag; always pin explicitly.

### Tag Format Note

BioContainers tags follow the Bioconda build naming convention:
```
<version>--<compiler_hash>_<build_number>
```
For centrifuger 1.1.0: `1.1.0--hf426362_0`

When a newer version is released, the tag will change. Do not drop the `--hf426362_0` suffix and expect the image to resolve — it will not.

### Compatibility with Existing Repo Patterns

| Concern | Assessment |
|---------|------------|
| WDL runtime block | Standard — `docker: docker`, `memory`, `cpu`, `disks` |
| Output from `centrifuger-kreport` | Kraken-style report format; compatible with the existing `parse_kraken2_reads` Python task if downstream parsing is needed |
| Multi-sample batch pattern | Same node-level DB load pattern as `krakenuniq` task in `tasks_metagenomics.wdl` — scatter classification across samples after single DB decompression |
| `miniwdl check` | No constraints from image choice; validation is purely WDL syntax |

---

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Inline heredoc pattern | HIGH | Direct codebase inspection, 7+ task examples |
| py3-bio package contents | HIGH | Confirmed via tasks using pandas, Bio, plotly, yaml, jinja2 |
| pip install at runtime | HIGH (decision), MEDIUM (mechanics) | Decision confirmed in PROJECT.md; exact shell invocation inferred from standard pip usage — not yet present in codebase |
| No `cat > script.py` pattern | HIGH | Full grep across tasks dir found zero matches |
| duckdb not in py3-bio | MEDIUM | Absence of evidence; py3-bio 0.1.3 not independently inspected |
| Centrifuger BioContainers image | HIGH | BioContainers TRS API confirms `1.1.0--hf426362_0` tag; Bioconda page confirms availability |
| Centrifuger binary names | HIGH | GitHub repo README and release notes confirm `centrifuger-*` prefix for all binaries |
| Broad images lack centrifuger | MEDIUM | viral-classify repo archived; confirmed via search + absence of centrifuger in any task file; no direct Dockerfile inspection |

---

## Sources

All v1.0–v2.0 findings are from direct file inspection. No web search was required.

- `/Users/carze/Documents/work/Broad_EMI/projects/viral-pipelines/.planning/PROJECT.md` — project decisions and validated patterns
- `/Users/carze/Documents/work/Broad_EMI/projects/viral-pipelines/requirements-modules.txt` — registered Docker images
- `pipes/WDL/tasks/tasks_utils.wdl:549-575` — `sanitize_fasta_headers` (minimal py3-bio pattern)
- `pipes/WDL/tasks/tasks_utils.wdl:915-942` — `tsv_drop_cols` (set -e + pandas pattern)
- `pipes/WDL/tasks/tasks_assembly.wdl:637-694` — `ivar_trim_stats` (pandas + plotly)
- `pipes/WDL/tasks/tasks_reports.wdl:225-285` — `merge_coverage_per_position` (Bio + pandas)
- `pipes/WDL/tasks/tasks_sarscov2.wdl:257-387` — `sars_cov_2_meta_etl` (~130 lines inline)
- `pipes/WDL/tasks/tasks_sarscov2.wdl:390-~550` — `crsp_meta_etl` (~140 lines inline)
- `pipes/WDL/tasks/tasks_ncbi.wdl:1086-~1250` — `generate_author_sbt_file` (multi-function inline)
- `pipes/WDL/tasks/tasks_demux.wdl:98-141` — `revcomp_i5` (Bio.Seq, short form)

v3.0 Centrifuger sources:
- [Bioconda centrifuger recipe](https://bioconda.github.io/recipes/centrifuger/README.html) — versions list, docker pull confirmed
- [BioContainers TRS API](https://api.biocontainers.pro/ga4gh/trs/v2/tools/centrifuger/versions) — exact tag `1.1.0--hf426362_0`, image size 133.9 MB
- [mourisl/centrifuger GitHub releases](https://github.com/mourisl/centrifuger/releases) — v1.1.0 released Feb 18 2025; binary names confirmed
- [mourisl/centrifuger GitHub README](https://github.com/mourisl/centrifuger) — full binary list, installation instructions
- `pipes/WDL/tasks/tasks_metagenomics.wdl` (direct read) — confirmed no centrifuger reference in any existing task
