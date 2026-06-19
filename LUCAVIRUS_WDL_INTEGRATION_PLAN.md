# LucaVirus WDL Integration Plan

Date: 2026-06-17

Branch: `feature/lucavirus`

Status: **Final implementation plan, pending external image tags.**

## Summary

Add LucaVirus support to viral-pipelines using the same split-runner design used for GPU classifiers such as VirNucPro:

- viral-ngs classify image owns preflight and post-processing helpers.
- `lucavirus-cuda` owns GPU model execution and bundled LucaVirus assets.
- viral-pipelines owns WDL task orchestration and backend runtime declarations.
- One workflow invocation runs one LucaVirus profile for one input FASTA.
- Multi-sample or multi-profile fanout happens upstream in a larger workflow.

The first WDL pass adds reusable tasks to `pipes/WDL/tasks/tasks_metagenomics.wdl` and a single user-facing workflow at `pipes/WDL/workflows/classify_lucavirus.wdl`.

## Source Contracts Reviewed

- LucaVirus CUDA container: `https://github.com/broadinstitute/lucavirus-cuda` at commit `98ed6e7f7ee51169705b220313285c562b33655b`.
- viral-ngs feature checkout: `/home/unix/carze/projects/viral-ngs`.
- In-tree GPU runtime precedent: `pipes/WDL/tasks/tasks_interhost.wdl::beast`.
- Historical split-runner classifier precedent: `feature-virnucpro` branch in this repository.

The user-provided upstream tool URL used singular `LucaVirusTask`, but the container repo currently pins and documents `LucaOne/LucaVirusTasks` with a trailing `s`. Treat the plural repository as the effective upstream contract unless the container repo changes.

## Design Decisions

- Do not install LucaVirus, PyTorch, CUDA libraries, or model assets into the viral-ngs classify image.
- Do not implement a viral-ngs `core.Tool` model runner for LucaVirus in this pass.
- Do not add a `_multi` workflow. Scatter over samples or profiles upstream.
- Do not translate nucleotide contigs or call ORFs in viral-pipelines.
- Treat the input FASTA as already containing the sequences LucaVirus should score. Current supported LucaVirus profiles are protein profiles, so WDL parameter metadata must call this an amino-acid/protein FASTA.
- Do not validate protein-vs-nucleotide content in WDL. This is a documented input contract because the current viral-ngs helper only converts FASTA records to LucaVirus CSV rows.
- Keep output LucaVirus-specific. Do not convert to Kraken2, Krona, or NCBI taxonomy reports in this pass.
- Do not support external LucaVirus model/profile localization in the first WDL pass. Use the bundled production `lucavirus-cuda` image.
- Do not expose a non-GPU LucaVirus execution mode. Non-empty LucaVirus
  classification requires CUDA/GPU acceleration; empty inputs skip the CUDA task
  before execution.

## External Prerequisites

These must be resolved before the WDL can run end to end:

1. Add and publish viral-ngs support for `metagenomics lucavirus_empty_predictions`.
2. Publish or identify a viral-ngs `*-classify` feature image containing `lucavirus_prepare`, `lucavirus_empty_predictions`, and `lucavirus_normalize`.
3. LucaVirus CUDA image: `ghcr.io/broadinstitute/lucavirus-cuda:v1.0`.

The viral-ngs helper image line in WDL must include `#skip-global-version-pin` while it points at a feature tag, because `requirements-modules.txt` currently pins `broadinstitute/viral-ngs=3.0.16` and the released `3.0.16-classify` image does not contain LucaVirus helpers.

## Remaining External Task Tracker

### viral-ngs Feature Branch

These are the remaining changes to complete on the viral-ngs side before the viral-pipelines WDL can be tested end to end:

- Add `metagenomics lucavirus_empty_predictions output.tsv`.
- Implement the command as a thin wrapper over `viral_ngs.classify.lucavirus.write_empty_predictions()`.
- Register the parser in `src/viral_ngs/metagenomics.py` using the same `cmd.common_args(...)`, `cmd.attach_main(...)`, and `__commands__.append(...)` pattern as `lucavirus_prepare` and `lucavirus_normalize`.
- Add parser coverage in `tests/unit/classify/test_metagenomics.py` that patches `viral_ngs.metagenomics.lucavirus.write_empty_predictions` and verifies the output argument is forwarded.
- Keep the existing `tests/unit/classify/test_lucavirus.py::test_write_empty_predictions` coverage for the helper function itself.
- Update `src/viral_ngs/classify/DESIGN.md` to list `metagenomics.py lucavirus_empty_predictions` as a WDL-facing entry point.
- Confirm `docker/Dockerfile.classify` import checks include `lucavirus` and that no CUDA, PyTorch, or LucaVirus model assets are added to the viral-ngs classify image.
- Build and publish a feature classify image containing all three WDL-facing commands:
  `lucavirus_prepare`, `lucavirus_empty_predictions`, and `lucavirus_normalize`.
- Record the published feature image tag so viral-pipelines can replace `quay.io/broadinstitute/viral-ngs:<feature-lucavirus-tag>-classify`.

Expected viral-ngs verification:

```bash
pytest tests/unit/classify/test_lucavirus.py tests/unit/classify/test_metagenomics.py
docker build -f docker/Dockerfile.classify -t <feature-viral-ngs-classify-image> .
docker run --rm <feature-viral-ngs-classify-image> metagenomics lucavirus_empty_predictions /tmp/empty.lucavirus.tsv
```

### lucavirus-cuda Container

The LucaVirus CUDA image is now available as:

```text
ghcr.io/broadinstitute/lucavirus-cuda:v1.0
```

Remaining container/backend checks:

- Pull the image on the intended execution hosts.
- Run a direct `docker run --gpus all` smoke with a prepared LucaVirus CSV and at least one supported profile.
- Confirm the target Cromwell backend can resolve GHCR image hashes. If it cannot, mirror the image to Quay or Artifact Registry and override `classify_lucavirus.lucavirus_docker`.
- Record the final production image reference, ideally including digest, once the backend path is confirmed.

### viral-pipelines Follow-Up

The WDL implementation has landed in commit `593f0176`, but it still uses the placeholder viral-ngs helper image. After the viral-ngs feature image exists:

- Replace `quay.io/broadinstitute/viral-ngs:feature-lucavirus-classify` in `tasks_metagenomics.wdl` and `classify_lucavirus.wdl` with the published feature image tag.
- Keep `#skip-global-version-pin` while using a feature tag not represented in `requirements-modules.txt`.
- Add an empty protein FASTA fixture and `test/input/WDL/miniwdl-local/test_inputs-classify_lucavirus-local.json`.
- Run the empty-input workflow test to validate prepare -> skip GPU -> normalize wiring without pulling the LucaVirus CUDA image.
- Run a non-empty local Cromwell GPU smoke with a real protein FASTA.
- Decide whether to add `.dockstore.yml` registration immediately after the image path is stable, or wait until the helper image is part of a released viral-ngs tag.

## Runtime Contract

The LucaVirus container provides:

```bash
/opt/lucavirus_cli.py input.csv output.tsv \
  --task-profile rdrp \
  --use-gpu \
  --gpu-id 0 \
  --verbose
```

The CLI also supports `--no-gpu`, `--model-path` / `--asset-path`,
`--profile-path`, and `--version`, but this WDL intentionally always passes
`--use-gpu` and uses the bundled model/profile paths in the image.

Supported profile values:

```text
rdrp
viral_capsid
virus_ec4
```

Normalized output header:

```text
seq_id	seq	prob	label_index	label
```

## viral-ngs Contract

The viral-ngs feature branch already contains:

- `metagenomics lucavirus_prepare`
- `metagenomics lucavirus_normalize`
- `viral_ngs.classify.lucavirus.prepare_contigs()`
- `viral_ngs.classify.lucavirus.normalize_output()`
- `viral_ngs.classify.lucavirus.write_empty_predictions()`

Before implementing WDL, add this stable CLI wrapper in viral-ngs:

```bash
metagenomics lucavirus_empty_predictions output.tsv
```

That command should call `write_empty_predictions()` so the LucaVirus TSV header remains single-sourced in viral-ngs. Do not call Python modules inline from WDL and do not hardcode the output header in shell.

## WDL Tasks To Add

Add the following tasks to `pipes/WDL/tasks/tasks_metagenomics.wdl`.

### `lucavirus_prepare`

Purpose: Convert input FASTA to LucaVirus CSV, write an empty raw-prediction fallback TSV, and produce stats for WDL branching.

Runtime image: viral-ngs classify image.

Inputs:

```wdl
File    input_fasta
String? sample_id
String  seq_type = "prot"

Int     machine_mem_gb = 4
String  docker = "quay.io/broadinstitute/viral-ngs:<feature-lucavirus-tag>-classify" #skip-global-version-pin
```

Private declaration:

```wdl
String out_basename = select_first([sample_id, basename(basename(basename(input_fasta, ".fasta"), ".fa"), ".fna")])
```

Key command shape:

```bash
set -e -o pipefail

metagenomics --version | tee VERSION

metagenomics lucavirus_prepare \
  "~{input_fasta}" \
  "~{out_basename}.lucavirus_input.csv" \
  "~{out_basename}.lucavirus_prepare_stats.tsv" \
  --seq-type "~{seq_type}" \
  --loglevel=DEBUG

tail -n +2 "~{out_basename}.lucavirus_prepare_stats.tsv" | cut -f 1 > N_SEQUENCES
tail -n +2 "~{out_basename}.lucavirus_prepare_stats.tsv" | cut -f 2 > HAS_INPUT

metagenomics lucavirus_empty_predictions \
  "~{out_basename}.raw.lucavirus.empty.tsv" \
  --loglevel=DEBUG
```

Outputs:

```wdl
String  output_basename           = out_basename
File    lucavirus_input_csv       = "~{out_basename}.lucavirus_input.csv"
File    prepare_stats_tsv         = "~{out_basename}.lucavirus_prepare_stats.tsv"
File    empty_raw_predictions_tsv = "~{out_basename}.raw.lucavirus.empty.tsv"
Int     n_sequences               = read_int("N_SEQUENCES")
Boolean has_lucavirus_input       = read_boolean("HAS_INPUT")
String  viralngs_version          = read_string("VERSION")
```

Implementation notes:

- Derive a private `out_basename` declaration from `sample_id` when present; otherwise strip `.fasta`, `.fa`, and `.fna` from `input_fasta`.
- Emit `output_basename` and consume it downstream. Do not recompute it in the workflow.
- Parameter metadata must state that `input_fasta` is amino-acid/protein FASTA and that ORF calling/translation must happen upstream.
- Parameter metadata should enumerate `seq_type = "prot"` as the only supported value for this pass.

### `lucavirus`

Purpose: Run the LucaVirus CUDA model for one profile.

Runtime image: standalone `lucavirus-cuda` image.

Inputs:

```wdl
File    lucavirus_input_csv
String  output_basename
String  task_profile = "rdrp"
Int     gpu_id = 0

Int     machine_mem_gb = 64
Int     cpu = 8
String? accelerator_type
Int?    accelerator_count
String? gpu_type
Int?    gpu_count
String? vm_size
Int     boot_disk_size_gb = 100
Int     disk_size_gb = 50
Int     preemptible_attempts = 0
String  docker = "ghcr.io/broadinstitute/lucavirus-cuda:v1.0"
```

Do not add `model_path_archive`, `profile_json`, `asset_path`, or `profile_path` inputs in the first pass. If externalized model assets are needed later, add a separate design because `--model-path` expects an extracted directory and has different disk/localization behavior.

Key command shape:

```bash
set -euo pipefail

/opt/lucavirus_cli.py --version | tee VERSION

export PATH="/usr/local/nvidia/bin:${PATH}"
if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi
else
  echo "WARNING: nvidia-smi is not available inside the container; letting LucaVirus validate CUDA availability." >&2
fi

/opt/lucavirus_cli.py \
  "~{lucavirus_input_csv}" \
  "~{output_basename}.raw.lucavirus.tsv" \
  --task-profile "~{task_profile}" \
  --gpu-id "~{gpu_id}" \
  --use-gpu \
  --verbose

test -s "~{output_basename}.raw.lucavirus.tsv"

{ if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi; } > MEM_BYTES
```

Outputs:

```wdl
File   raw_predictions_tsv = "~{output_basename}.raw.lucavirus.tsv"
Int    max_ram_gb          = ceil(read_float("MEM_BYTES")/1000000000)
String lucavirus_version   = read_string("VERSION")
```

Runtime:

```wdl
runtime {
  docker: docker
  memory: "~{machine_mem_gb} GB"
  cpu: cpu
  disks: "local-disk ~{disk_size_gb} LOCAL"
  disk: "~{disk_size_gb + boot_disk_size_gb} GB"
  bootDiskSizeGb: boot_disk_size_gb
  gpu: true
  acceleratorType: select_first([accelerator_type, "nvidia-tesla-t4"])
  acceleratorCount: select_first([accelerator_count, gpu_count, 1])
  gpuType: select_first([gpu_type, "nvidia-tesla-t4"])
  gpuCount: select_first([gpu_count, accelerator_count, 1])
  nvidiaDriverVersion: "410.79"
  vm_size: select_first([vm_size, "Standard_NC6s_v3"])
  dx_instance_type: "mem2_ssd1_gpu1_x8"
  preemptible: preemptible_attempts
  maxRetries: 2
}
```

Implementation notes:

- Use `preemptible_attempts = 0` by default for the first cut. LucaVirus has no checkpoint/resume behavior, and preemption repeats model load and inference from scratch.
- Keep `maxRetries: 2` for transient infrastructure failures.
- Keep `nvidiaDriverVersion` because PAPIv2 may require it; Google Batch can ignore or supersede it depending on backend behavior.
- Confirm `dx_instance_type: "mem2_ssd1_gpu1_x8"` before relying on DNAnexus deployment. If unverified, use the known `beast` value or mark DNAnexus unsupported for the first pass.
- Validate `boot_disk_size_gb = 100` against the published image size. Increase it if the bundled image plus worker overhead approaches the boot disk limit.
- Default image is `ghcr.io/broadinstitute/lucavirus-cuda:v1.0`. If a Cromwell backend cannot resolve GHCR image hashes, override `lucavirus_docker` with a mirrored Quay/Artifact Registry tag.

### `lucavirus_normalize`

Purpose: Validate raw LucaVirus output and emit the durable user-facing TSV.

Runtime image: viral-ngs classify image.

Inputs:

```wdl
File   raw_predictions_tsv
String output_basename
String task_profile = "rdrp"

Int    machine_mem_gb = 4
String docker = "quay.io/broadinstitute/viral-ngs:<feature-lucavirus-tag>-classify" #skip-global-version-pin
```

Key command shape:

```bash
set -e -o pipefail

metagenomics --version | tee VERSION

metagenomics lucavirus_normalize \
  "~{raw_predictions_tsv}" \
  "~{output_basename}.lucavirus.tsv" \
  --task-profile "~{task_profile}" \
  --loglevel=DEBUG
```

Outputs:

```wdl
File   predictions_tsv  = "~{output_basename}.lucavirus.tsv"
String viralngs_version = read_string("VERSION")
```

## Workflow To Add

Add `pipes/WDL/workflows/classify_lucavirus.wdl`.

Workflow inputs:

```wdl
File    input_fasta
String? sample_id
String  task_profile = "rdrp"
Int     gpu_id = 0

Int     helper_machine_mem_gb = 4
String  viralngs_docker = "quay.io/broadinstitute/viral-ngs:<feature-lucavirus-tag>-classify" #skip-global-version-pin

Int     machine_mem_gb = 64
Int     cpu = 8
String  lucavirus_docker = "ghcr.io/broadinstitute/lucavirus-cuda:v1.0"
String? accelerator_type
Int?    accelerator_count
String? gpu_type
Int?    gpu_count
String? vm_size
Int     boot_disk_size_gb = 100
Int     disk_size_gb = 50
Int     preemptible_attempts = 0
```

Workflow body:

```wdl
call metagenomics.lucavirus_prepare {
  input:
    input_fasta    = input_fasta,
    sample_id      = sample_id,
    machine_mem_gb = helper_machine_mem_gb,
    docker         = viralngs_docker
}

if (lucavirus_prepare.has_lucavirus_input) {
  call metagenomics.lucavirus {
    input:
      lucavirus_input_csv = lucavirus_prepare.lucavirus_input_csv,
      output_basename     = lucavirus_prepare.output_basename,
      task_profile        = task_profile,
      gpu_id              = gpu_id,
      machine_mem_gb      = machine_mem_gb,
      cpu                 = cpu,
      docker              = lucavirus_docker,
      accelerator_type    = accelerator_type,
      accelerator_count   = accelerator_count,
      gpu_type            = gpu_type,
      gpu_count           = gpu_count,
      vm_size             = vm_size,
      boot_disk_size_gb   = boot_disk_size_gb,
      disk_size_gb        = disk_size_gb,
      preemptible_attempts = preemptible_attempts
  }
}

call metagenomics.lucavirus_normalize {
  input:
    raw_predictions_tsv = select_first([
      lucavirus.raw_predictions_tsv,
      lucavirus_prepare.empty_raw_predictions_tsv
    ]),
    output_basename     = lucavirus_prepare.output_basename,
    task_profile        = task_profile,
    machine_mem_gb      = helper_machine_mem_gb,
    docker              = viralngs_docker
}
```

Workflow outputs:

```wdl
File    lucavirus_predictions_tsv      = lucavirus_normalize.predictions_tsv
File    lucavirus_prepare_stats_tsv    = lucavirus_prepare.prepare_stats_tsv
Int     lucavirus_input_sequences      = lucavirus_prepare.n_sequences
String  viralngs_version               = lucavirus_normalize.viralngs_version
String? lucavirus_version              = lucavirus.lucavirus_version
Int?    lucavirus_max_ram_gb           = lucavirus.max_ram_gb
```

Do not expose the prepared CSV or raw TSV by default.

## Parameter Metadata Requirements

Add explicit metadata for:

- `input_fasta`: amino-acid/protein FASTA only; nucleotide contigs must be translated or ORF-called upstream.
- `task_profile`: choices `rdrp`, `viral_capsid`, `virus_ec4`.
- `gpu_id`: GPU index passed to LucaVirus. Non-empty execution always uses GPU acceleration.
- `viralngs_docker`: feature classify image until LucaVirus helpers are in a released viral-ngs image.
- `lucavirus_docker`: standalone CUDA image with bundled assets.
- GPU runtime fields: backend-specific fields for GCP/PAPIv2, Terra, Azure/TES, and DNAnexus.

## CI And Test Inputs

The full GPU workflow should not run in GitHub CI. The basic CI path should exercise the empty-input branch, which skips the GPU call while still validating workflow wiring.

Add a tiny empty FASTA fixture and a miniwdl-local workflow input JSON:

```text
test/input/lucavirus.empty.protein.fasta
test/input/WDL/miniwdl-local/test_inputs-classify_lucavirus-local.json
```

The miniwdl-local JSON should use the feature viral-ngs classify image and any placeholder `lucavirus_docker` value. The GPU call is skipped because `has_lucavirus_input=false`, so Docker should not pull the LucaVirus image.

Example empty-input JSON shape:

```json
{
  "classify_lucavirus.input_fasta": "test/input/lucavirus.empty.protein.fasta",
  "classify_lucavirus.sample_id": "lucavirus_empty",
  "classify_lucavirus.task_profile": "rdrp",
  "classify_lucavirus.viralngs_docker": "quay.io/broadinstitute/viral-ngs:<feature-lucavirus-tag>-classify",
  "classify_lucavirus.lucavirus_docker": "ghcr.io/broadinstitute/lucavirus-cuda:v1.0"
}
```

Add expected outputs only for stable scalar outputs. Avoid matching localized output file paths.

Optional local-only task checks:

```bash
miniwdl run --task lucavirus_prepare pipes/WDL/tasks/tasks_metagenomics.wdl \
  input_fasta=test/input/lucavirus.empty.protein.fasta \
  sample_id=lucavirus_empty \
  docker=<feature-viral-ngs-classify-image>

miniwdl run --task lucavirus_normalize pipes/WDL/tasks/tasks_metagenomics.wdl \
  raw_predictions_tsv=<header-only-or-two-row-raw-tsv> \
  output_basename=lucavirus_normalize_test \
  task_profile=rdrp \
  docker=<feature-viral-ngs-classify-image>
```

Add `classify_lucavirus` to `.dockstore.yml` only after a public/default LucaVirus Docker image is available.

## Validation Plan

Minimum repository validation:

```bash
miniwdl check pipes/WDL/workflows/classify_lucavirus.wdl
womtool validate pipes/WDL/workflows/classify_lucavirus.wdl
MODULE_VERSIONS=./requirements-modules.txt github_actions_ci/check-wdl-runtimes.sh
github_actions_ci/tests-miniwdl.sh
git diff --check
```

Manual container validation before full WDL:

```bash
docker run --rm --gpus all \
  -v "$PWD":/data \
  ghcr.io/broadinstitute/lucavirus-cuda:v1.0 \
  /opt/lucavirus_cli.py /data/input.csv /data/output.tsv \
  --task-profile rdrp \
  --use-gpu \
  --gpu-id 0 \
  --verbose
```

Local/GPU Cromwell validation outside GitHub CI:

```bash
cromwell -Dconfig.file=pipes/cromwell/cromwell.local-gpu.conf \
  run \
  pipes/WDL/workflows/classify_lucavirus.wdl \
  -i test/input/WDL/cromwell-gpu/test_inputs-classify_lucavirus-local_gpu.json \
  -m lucavirus_gpu_smoke.metadata.json
```

This is the required path for non-empty LucaVirus validation. MiniWDL can
exercise the empty branch because it skips the CUDA task, but a non-empty
LucaVirus run requires GPU acceleration and should be validated through
Cromwell with Docker GPU exposure.

Local GPU setup files:

- `pipes/cromwell/cromwell.local-gpu.conf`: local Cromwell backend config that
  passes `--gpus all` to Docker.
- `test/input/WDL/cromwell-gpu/test_inputs-classify_lucavirus-local_gpu.json`:
  non-empty LucaVirus Cromwell input JSON.
- `test/input/lucavirus.gpu-smoke.protein.fasta`: small synthetic protein FASTA
  for smoke testing pipeline wiring. Replace with a known positive/real protein
  fixture once available.
- `test/input/lucavirus_capsid_panel/`: UniProt-derived `viral_capsid` smoke
  panel split into capsid-positive, viral non-capsid, and non-viral protein
  FASTA controls, plus metadata.
- `test/input/WDL/cromwell-gpu/test_inputs-classify_lucavirus-capsid_panel-local_gpu.json`:
  local Cromwell input JSON for running the combined capsid panel with
  `task_profile = "viral_capsid"`.

Current local readiness:

- `cromwell 89` and `womtool 89` are available on `PATH`.
- Docker reports an `nvidia` runtime.
- Codex sandboxed command execution cannot see `/dev/nvidia*`, so sandboxed
  `nvidia-smi` fails. Unsandboxed host execution sees two NVIDIA GeForce RTX
  4090 GPUs with driver `590.48.01` and CUDA `13.1`.
- `ghcr.io/broadinstitute/lucavirus-cuda:v1.0` is not currently present locally
  and must be pulled or otherwise made available before the full smoke run. An
  anonymous `docker pull ghcr.io/broadinstitute/lucavirus-cuda:v1.0` currently
  returns `denied`, so authenticate to GHCR or mirror the image to a registry
  this host can access.

## Remaining Implementation Sequence

Current viral-pipelines status:

- Done: `lucavirus_prepare`, `lucavirus`, and `lucavirus_normalize` are implemented in `tasks_metagenomics.wdl`.
- Done: `classify_lucavirus.wdl` is implemented.
- Done: WDL validation passes with `miniwdl check`, `womtool validate`, `check-wdl-runtimes.sh`, and `git diff --check`.
- Done locally: viral-ngs now provides `metagenomics lucavirus_prepare`, `metagenomics lucavirus_empty_predictions`, and `metagenomics lucavirus_normalize`.
- Done locally: a local viral-ngs classify image was built as `quay.io/broadinstitute/viral-ngs:feature-lucavirus-classify`.
- Done locally: an empty-input miniwdl fixture exercises prepare -> skip GPU -> normalize using the local viral-ngs image; the expected scalar output JSON matches the run output.

Remaining sequence:

1. Push or otherwise publish a viral-ngs feature classify image after the branch has a stable commit/tag.
2. Replace any temporary local-only image references with the published viral-ngs feature image tag.
3. Confirm the target Cromwell backend can resolve `ghcr.io/broadinstitute/lucavirus-cuda:v1.0`; mirror and override `lucavirus_docker` if needed.
4. Validate the miniwdl-local test through the repository test harness, not only direct `miniwdl run`.
5. Run a manual `docker run --gpus all` LucaVirus smoke with a real non-empty protein FASTA.
6. Run a local Cromwell GPU smoke with a real non-empty protein FASTA.
7. Add Dockstore registration after the production image path is stable.

## Open Questions

- Should a later pass support external model/profile localization through `--asset-path` and `--profile-path`?
- Should future reporting add a LucaVirus-specific summary table over `seq_id`, `prob`, `label_index`, and `label`?
- Should a later end-to-end workflow run all three profiles by scattering over `["rdrp", "viral_capsid", "virus_ec4"]`, or should profile selection remain user-driven?
- What final image registry and tag should become the default for Cromwell/Terra?
- Is `dx_instance_type: "mem2_ssd1_gpu1_x8"` valid for DNAnexus, or should the first pass use the known BEAST GPU instance type?
