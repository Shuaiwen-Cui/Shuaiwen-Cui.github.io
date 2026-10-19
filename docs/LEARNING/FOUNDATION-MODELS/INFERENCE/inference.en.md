# Inference & Deployment

!!! abstract "What this page solves"
    Once a model is trained, how to make it fast, small and still good on given hardware. This is your home ground: quantisation and memory work on MCUs carries over to LLMs, only at a different scale. **Basics required, depth by direction**; prerequisite: you can load and run a pretrained model.

## How far to go

| Level | Learn | Done when |
|---|---|---|
| Orientation | Prefill vs decode, time-to-first-token, throughput, memory | You can explain why TTFT and generation speed differ |
| Foundation | Measuring under fixed conditions, quantisation comparison | F5 done with quality, latency and memory reported together |
| Advanced | KV-cache management, continuous batching, speculative decoding, kernels | A bottleneck located by measurement, gains and limits explained |

## Core knowledge

**Two phases**: prefill processes the whole prompt in parallel, mostly compute-bound, and sets time-to-first-token (TTFT); decode generates one token at a time, reads all weights every step, is mostly memory-bandwidth-bound, and sets tokens/s.

**Back-of-envelope decode speed**: at batch 1, tokens/s ≲ memory bandwidth ÷ weight bytes. A 7B model in FP16 is ~14 GB, so on 100 GB/s the ceiling is ~7 tokens/s; at 4-bit, ~3.5 GB and ~28 tokens/s. That is why quantisation speeds things up.

**KV cache**: caches past tokens' K and V to avoid recomputation; size = 2 × layers × KV heads × head dim × length × batch × bytes. With long contexts or many users it outgrows the weights. GQA / MQA share KV heads to shrink it; PagedAttention removes fragmentation.

**Decoding strategies**: greedy, temperature sampling, top-k, top-p (nucleus), beam search. Higher temperature means more diversity and more errors; low for factual tasks, high for creative ones.

**Quantisation**: weight-only quantisation (GPTQ, AWQ, GGUF k-quants) is most common and 4-bit usually costs acceptable quality; weight-and-activation quantisation (e.g. W8A8) is harder because activations have outliers — SmoothQuant moves the difficulty from activations to weights.

**Distillation and pruning**: distillation trains a small model on a large model's outputs; pruning removes unimportant weights or structures (only structured pruning truly speeds up).

**Speculative decoding**: a small draft model guesses several tokens, the large model verifies them in one parallel pass and accepts the correct prefix; the output distribution is exactly the large model's, and speed depends on acceptance rate.

**Continuous batching and serving**: servers batch requests at token granularity to keep the GPU busy; vLLM, TensorRT-LLM and SGLang are typical.

**Edge runtimes**: llama.cpp (CPU / GPU / Apple Silicon, GGUF), ONNX Runtime, TensorRT, MLC-LLM. At the edge, peak memory, power and TTFT matter most.

## Debate · general acceleration vs. hardware-specific optimisation

General: rely on mature frameworks and the GPU ecosystem, build once, deploy anywhere. Specific: operator fusion and custom quantisation for an NPU or chip win orders of magnitude in power and cost, at a high porting cost. Edge and automotive increasingly go specific — where an embedded background pays off.

## Curated resources

<div class="grid cards" markdown>

-   :material-github:{ .lg .middle } __Project entry · llama.cpp__

    ---

    Run small models on a PC and a Jetson; compare quantisations of one model (F5).

    [:octicons-arrow-right-24: Open repository](https://github.com/ggml-org/llama.cpp)

-   :material-school:{ .lg .middle } __Direction course · MIT 6.5940 TinyML__

    ---

    Quantisation, pruning, distillation and efficient inference, systematically; once efficient inference is your direction.

    [:octicons-arrow-right-24: Open course](https://hanlab.mit.edu/courses/2024-fall-65940)

</div>

## Follow-along experiment

<div class="grid cards" markdown>

-   :material-hammer-wrench:{ .lg .middle } __F5 · Inference comparison on a fixed task__

    ---

    Several quantisations of one model; record TTFT, tokens/s, peak memory and quality regressions.

    [:octicons-arrow-right-24: See steps](../PROJECT-BASED-LEARNING/project-based-learning.md)

</div>

## Self-check

??? question "Why does total decode throughput rise with batch size?"
    Batch-1 decode is bandwidth-bound: one read of the weights yields one token. Larger batches serve more tokens per read and raise compute utilisation until compute becomes the limit.

??? question "A 7B model at 4-bit on a 50 GB/s device: roughly what batch-1 decode speed?"
    About 3.5 GB of weights, so a ceiling of ~50 / 3.5 ≈ 14 tokens/s; real numbers are lower.

??? question "Does speculative decoding change output quality?"
    No. Rejection sampling keeps the output distribution identical to the large model's; only speed changes.

??? question "Why is activation quantisation harder than weight quantisation?"
    Activations change with the input and some channels have huge outliers, so one scale leaves most values with too little precision.

## Related topics

<div class="grid cards" markdown>

-   __Efficient Training & Engineering (AI)__

    ---

    Quantisation, distillation and pruning for vision and general models.

    [:octicons-arrow-right-24: Open](../../AI/SYSTEMS/systems.md)

-   __Sim-to-Real (Embodied)__

    ---

    Onboard compute and real-time constraints.

    [:octicons-arrow-right-24: Open](../../EMBODIED-INTELLIGENCE/SIM-TO-REAL/sim-to-real.md)

</div>
