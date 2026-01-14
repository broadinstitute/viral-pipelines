# VirNucPro GPU Task Architecture

## Architecture

```
User invokes workflow (Terra/GCP/DNAnexus/Azure)
         |
         v
+------------------------------------------+
|  classify_virnucpro_multi.wdl            |
|                                          |
|  scatter(bam in bams) {                 |
|    call virnucpro.classify_virnucpro    |
|  }                                       |
+------------------------------------------+
         | N parallel task invocations
         | (1 GPU instance per BAM)
         v
+------------------------------------------+
|  tasks_metagenomics.wdl                  |
|  task classify_virnucpro {               |
|                                          |
|    runtime {                             |
|      gpu: true                           |
|      acceleratorType: nvidia-tesla-p4    |
|      gpuType: nvidia-tesla-p4            |
|      dx_instance_type: mem1_ssd1_gpu2_x8 |
|      vm_size: Standard_NC6               |
|    }                                     |
|                                          |
|    command {                             |
|      virnucpro_cli.py $bam $tsv \        |
|        --expected-length 300             |
|    }                                     |
|  }                                       |
+------------------------------------------+
         | Each instance runs:
         |   BAM -> FASTA (pysam)
         |   FASTA dedup
         |   VirNucPro prediction.py
         |   (DNABERT_S + ESM2-3B models)
         |   GPU acceleration
         v
+------------------------------------------+
|  Output: TSV per BAM                     |
|  Sequence_ID  Prediction  score1 score2 |
|  read1        virus       0.95   0.05   |
|  read2        non-virus   0.12   0.88   |
+------------------------------------------+
```

## Data Flow

```
Input: Array[BAM] (unaligned reads)
    |
    | WDL scatter block creates N parallel jobs
    v
N x GPU Instances (independent, isolated)
    |
    | Each instance:
    |   1. BAM -> FASTA conversion (pysam)
    |   2. FASTA ID deduplication (_ensure_unique_fasta_ids)
    |   3. VirNucPro subprocess: prediction.py
    |      - DNABERT_S embedding (nucleotide)
    |      - ESM2-3B embedding (amino acid)
    |      - MLP classification
    |   4. Parse results to TSV
    v
Output: Array[TSV] (1:1 correspondence with input BAMs)
```

## Why This Structure

**Task file organization**: VirNucPro is a viral sequence classifier, grouped with other classification tools (Kraken2, KrakenUniq, Kaiju, kb, BLASTx) in `tasks_metagenomics.wdl`. While VirNucPro uses GPU acceleration (unlike the CPU-based classifiers), the conceptual similarity as a classifier outweighs the runtime difference. This co-location improves discoverability: users looking for classification tools find all options in one file.

**Two workflows**: Single-sample workflow (`classify_virnucpro_single.wdl`) simplifies testing and debugging. Multi-sample workflow (`classify_virnucpro_multi.wdl`) provides production capability. Mirrors existing pattern: `classify_kraken2.wdl` (single) vs `classify_multi.wdl` (multi).

**Scatter over batch**: Each BAM gets independent GPU instance. Rationale: BEAST pattern prioritizes reliability over cost. Failure isolation critical for production (one corrupt BAM doesn't fail entire batch). GPU startup overhead (~1-3 min) acceptable compared to unknown VirNucPro runtime (requires benchmarking). Batch processing can be added later if cost analysis justifies development effort.

**Multi-platform runtime attributes**: Different platforms use different attribute names for GPU specification. WDL runtime block includes ALL platform attributes; each platform ignores attributes it doesn't recognize. `select_first([user_override, default_value])` pattern allows user control while providing sensible defaults. Proven pattern from BEAST task (pipes/WDL/tasks/tasks_interhost.wdl).

## Invariants

**Platform attribute completeness**: Runtime block MUST include GPU attributes for ALL supported platforms (GCP: acceleratorType/acceleratorCount, Terra: gpuType/gpuCount, DNAnexus: gpu/dx_instance_type, Azure: vm_size). Missing platform attributes cause silent fallback to CPU or workflow failure.

**Expected length must match model**: VirNucPro has two models (300_model.pth, 500_model.pth). `expected_length` parameter (300 or 500) must match model file. Mismatch produces invalid predictions. VirNucPro prediction.py silently accepts wrong length but results are meaningless.

**GPU driver version pinned**: `nvidiaDriverVersion: "410.79"` matches BEAST task. Driver version must be compatible with CUDA 11.8 (VirNucPro Docker requirement). Newer drivers are backward-compatible, but older drivers may lack CUDA 11.8 support.

**Task outputs correspond to inputs**: Scatter workflow produces `Array[File]` outputs in same order as `Array[File]` inputs. Downstream workflows depend on 1:1 BAM-to-TSV correspondence. basename preservation maintains traceability.

## Tradeoffs

**Speed vs Cost**: Scatter approach (1 BAM = 1 GPU instance) maximizes wall-clock speed through parallelism, but incurs N x GPU startup overhead. Batch approach (N BAMs = 1 GPU instance) would reduce cost but serialize processing and worsen failure isolation. Decision: prioritize speed and reliability; cost optimization deferred pending runtime benchmarks.

**Multi-platform portability vs Platform optimization**: Using platform-agnostic runtime attributes (same GPU type across all platforms) simplifies maintenance but may not leverage platform-specific optimizations (e.g., Terra-optimized instances, DNAnexus mem1_ssd1_gpu2_x8 vs mem2_ssd1_gpu2_x8). Decision: portability prioritized; users can override via optional inputs for platform-specific tuning.

**Fail-fast vs Graceful degradation**: No CPU fallback means workflow fails if GPU unavailable. Alternative would be automatic retry with --no-gpu flag. Decision: explicit failure provides clear feedback and prevents silent performance degradation (CPU can be 10-100x slower). Users expecting GPU resources should know immediately if unavailable.

**Test coverage vs CI speed**: GitHub Actions CI lacks GPU hardware, so VirNucPro workflows are excluded from automated execution tests. WDL syntax validation runs in CI via `miniwdl check`, but functional testing requires manual execution on GPU-enabled infrastructure (Terra, GCP, DNAnexus, or Azure). Test with small BAM files first to validate workflow structure before production runs.

## Testing

**CI Validation**: WDL syntax is validated automatically via `miniwdl check` in GitHub Actions.

**Manual Testing Required**: Functional testing must be performed on GPU-enabled platforms:

```bash
# Terra/Cromwell: Upload workflow and use Terra UI
# GCP with Cromwell:
java -jar cromwell.jar run pipes/WDL/workflows/classify_virnucpro_single.wdl \
  -i test_inputs.json

# DNAnexus:
dx run classify_virnucpro_single -i reads_bam=project-xxx:file-yyy

# Local with GPU (requires nvidia-docker):
miniwdl run --runtime.gpu=true pipes/WDL/workflows/classify_virnucpro_single.wdl \
  reads_bam=test.bam expected_length=300
```

**Test Input Example**:
```json
{
  "classify_virnucpro_single.reads_bam": "gs://bucket/sample.bam",
  "classify_virnucpro_single.expected_length": 300
}
```
