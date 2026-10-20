# Project · Quantizing a Small Language Model: from RTN to GPTQ / AWQ

!!! abstract "What this project solves"
    You already do FP8 / int8 quantization on MCUs; this project moves the same skill to large models. On a 0.5B model, **compare several quantization algorithms on accuracy, memory and speed by hand, and explain where the accuracy loss comes from**. It fills the gap between "tuning llama.cpp flags" and "understanding quantization algorithms". One NVIDIA GPU with 8 GB+ memory is enough; about 5–7 days.

## How this differs from F5

| | F5 · Inference comparison | This project |
|---|---|---|
| Focus | Runtime: how fast each quantized build runs | Algorithms: why each method works or fails |
| Where quantization comes from | Pre-converted GGUF files | Your own calibration and quantization |
| Deliverable | Speed and memory table | Accuracy–memory–speed comparison + per-layer error analysis |

## Target level

| Skill | After the project you can |
|---|---|
| Quantization basics | Explain symmetric / asymmetric, per-tensor / per-channel / group-wise, W8A8 vs W4A16 |
| Post-training quantization (PTQ) | Quantize once each with RTN, GPTQ, AWQ and SmoothQuant and say what each fixes |
| Error analysis | Find the most sensitive layers; explain why activation outliers make activation quantization hard |
| Inference | Measure prefill and per-token decode speed; know what KV-cache quantization saves |
| Quantization-aware training (QAT, optional) | Compare PTQ and QAT on a small task |

## Steps

1. **Baseline**: Qwen2.5-0.5B-Instruct in FP16. Fixed evaluation: WikiText-2 perplexity + one small task (e.g. ARC-Easy or a HellaSwag subset in lm-evaluation-harness) + a manual check on 20 fixed prompts. Record memory, prefill and decode speed.
2. **Plain rounding (RTN)**: write a few dozen lines yourself to quantize linear-layer weights per channel and per group (group size 128) to INT8 and INT4, dequantize, and evaluate. No library — the point is to see where the error comes from.
3. **GPTQ and AWQ**: W4A16 with llm-compressor, 128 calibration samples. Compare perplexity and task score, then swap the calibration set for text from another domain and watch what changes.
4. **Activation quantization and SmoothQuant**: W8A8 first without smoothing, then with SmoothQuant. Plot per-channel activation maxima for a few layers to see the outlier channels.
5. **Per-layer sensitivity**: quantize one layer at a time (others FP16), record the perplexity change and plot sensitivity; keep the most sensitive layers in FP16 and try mixed precision.
6. **Speed**: serve two quantized results with vLLM (and once on a Jetson if available); measure prefill, decode and memory; then enable FP8 / INT8 KV cache and see the memory saved at long context.
7. **Optional · QAT**: fine-tune with torchao's QAT flow on a small downstream task; compare INT4 PTQ with QAT.

## Deliverables

- One summary table: method × bit width → perplexity, task score, memory, prefill / decode speed.
- Two plots: activation outlier channels, per-layer sensitivity.
- A short "accuracy-drop debugging" note: which layer is most sensitive, why, and how to fix it.
- Code repository and README with model version, calibration set, hardware and repetitions; mark unfinished parts as "planned".

## Curated resources

<div class="grid cards" markdown>

-   :material-tools:{ .lg .middle } __Main tool · llm-compressor__

    ---

    vLLM's quantization library with ready examples for GPTQ, AWQ, SmoothQuant and FP8; results load directly in vLLM.

    [:octicons-arrow-right-24: Open repository](https://github.com/vllm-project/llm-compressor)

-   :material-fire:{ .lg .middle } __QAT and low-bit · torchao__

    ---

    PyTorch's quantization and sparsity library, with INT4 / INT8 weight quantization and a QAT flow.

    [:octicons-arrow-right-24: Open repository](https://github.com/pytorch/ao)

-   :material-clipboard-check:{ .lg .middle } __Evaluation · lm-evaluation-harness__

    ---

    One harness for perplexity and common benchmarks, so before and after use the same ruler.

    [:octicons-arrow-right-24: Open repository](https://github.com/EleutherAI/lm-evaluation-harness)

</div>

**Three must-read papers** (method section and main figures only):

<div class="grid cards" markdown>

-   __GPTQ__

    ---

    Column-by-column quantization with second-order error compensation.

    [:octicons-arrow-right-24: arXiv 2210.17323](https://arxiv.org/abs/2210.17323)

-   __AWQ__

    ---

    Protects salient weight channels found from activation magnitude.

    [:octicons-arrow-right-24: arXiv 2306.00978](https://arxiv.org/abs/2306.00978)

-   __SmoothQuant__

    ---

    Migrates activation outliers into the weights so W8A8 becomes feasible.

    [:octicons-arrow-right-24: arXiv 2211.10438](https://arxiv.org/abs/2211.10438)

</div>

## Self-check

??? question "Why do LLMs usually use weight-only quantization (W4A16) instead of quantizing weights and activations?"
    LLM activations contain a few channels with very large values (outliers); per-tensor quantization starves the other channels of precision. Decode is mostly memory-bandwidth bound, so compressing weights alone captures most of the speed gain.

??? question "What fundamentally separates GPTQ from RTN?"
    RTN rounds each weight independently; after quantizing a column, GPTQ uses a Hessian approximation to spread that column's error onto columns not yet quantized, reducing overall output error.

??? question "Why do results change when the calibration set comes from another domain?"
    GPTQ, AWQ and SmoothQuant all rely on activation statistics from calibration data; if that distribution differs from real use, the chosen salient channels and scales are off.

??? question "What does KV-cache quantization mainly save?"
    At long context and large batch the KV cache can exceed the weights in memory; quantizing it mainly saves memory and enlarges usable context and batch, with little effect on single short-prompt speed.

## Related

<div class="grid cards" markdown>

-   __Inference & deployment__

    ---

    Principles of quantization, KV cache and serving systems.

    [:octicons-arrow-right-24: Open](../INFERENCE/inference.md)

-   __F5 · Inference comparison__

    ---

    Runtime-level comparison; an extension of step 6 here.

    [:octicons-arrow-right-24: Open](project-based-learning.md#f5)

-   __AI · Efficient training & engineering__

    ---

    Timing methods, compute-bound vs memory-bound.

    [:octicons-arrow-right-24: Open](../../AI/SYSTEMS/systems.md)

</div>
