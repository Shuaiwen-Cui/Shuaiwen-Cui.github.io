# Efficient Training, Inference & AI Engineering

!!! abstract "What this page solves"
    Turn algorithms into systems that are reproducible, measurable and deployable. This is your home ground: embedded and edge-deployment experience carries over directly to resource budgets, profiling and quantization; fill GPU training and inference runtimes as real gaps appear. **P0: engineering and reproducibility required; P1: compression and performance.**

## Knowledge map

| Level | Topic | What to understand |
|---|---|---|
| P0 | Engineering and reproducibility | Data versions, random seeds, configs, dependencies, checkpoints, experiment logs |
| P1 | Performance and training efficiency | Profiling, mixed precision, batch size, gradient accumulation, data loading and memory |
| P1 | Model compression | Quantization, distillation, pruning and sparsity; verify quality and real speed separately |
| P2 | Large-model systems | Efficient attention, KV cache, batching, parallelism and communication |
| P2 | Deployment and lifecycle | ONNX and runtimes, edge devices, monitoring, drift, regression evaluation |

## Core mechanics

**Reproducibility**: fix random seeds; record code version, dependencies, config and data version. Without these, "one hyperparameter change made it better" cannot be confirmed.

**Profiling**: measure before optimising. Separate data loading, forward, backward and optimiser step; GPU timing needs synchronisation, and the first iterations are warm-up and excluded. Many training runs are actually bottlenecked by data loading.

**Mixed precision**: forward and backward in FP16 / BF16 with FP32 master weights. BF16 has FP32's exponent bits, rarely overflows, and is preferred on modern GPUs; FP16 needs loss scaling.

**Quantization**: move weights (and activations) from float to INT8 or lower, shrinking the model and bandwidth needs. Post-training quantization (PTQ) sets scales with calibration data; quantization-aware training (QAT) simulates quantization error during training. Per-channel is more accurate than per-tensor.

**Distillation**: a small model imitates a large model's soft outputs (probability distributions), which carry more information than hard labels.

**Pruning**: remove unimportant weights. Unstructured pruning compresses a lot but is hard for hardware to accelerate; structured pruning (removing channels, attention heads) is what actually speeds things up.

**Compute-bound vs memory-bound**: whether an operator is fast depends on arithmetic intensity (operations per byte read). Many layers are bandwidth-bound, so fewer FLOPs need not mean lower latency — the same lesson as optimising kernels on an MCU.

**Deployment**: export to ONNX, then let runtimes such as TensorRT or ONNX Runtime fuse operators and quantize; on device also check peak memory, operator support and power.

## Curated resources

<div class="grid cards" markdown>

-   :material-chip:{ .lg .middle } __Direction course · MIT 6.5940 TinyML__

    ---

    Efficient deep learning: quantization, pruning, distillation, neural architecture search and hardware deployment.

    [:octicons-arrow-right-24: Open course page](https://hanlab.mit.edu/courses/2024-fall-65940)

-   :material-timer-outline:{ .lg .middle } __Measurement · PyTorch Profiler__

    ---

    Locate compute bottlenecks, memory cost and operator latency on a real model.

    [:octicons-arrow-right-24: Open guide](https://docs.pytorch.org/tutorials/recipes/recipes/profiler_recipe.html)

-   :material-speedometer:{ .lg .middle } __Quantization practice · ONNX Runtime__

    ---

    Dynamic and static quantization, calibration and post-quantization accuracy debugging.

    [:octicons-arrow-right-24: Open docs](https://onnxruntime.ai/docs/performance/model-optimizations/quantization.html)

</div>

## Follow-along · AI-6 · An explainable performance experiment

1. Reuse a model already running from F1, AI-2 or F4; no new task.
2. Fix device, input and precision; separate loading, warm-up, forward and data processing.
3. Time with proper device synchronisation, repeat, and record the latency distribution and peak memory.
4. Change only batch, precision or one optimisation, then recheck quality and real resources.

**Done when**: the report states environment, measurement scope, repetitions and quality limits; theoretical FLOPs or parameter count are not deployment latency.

## Self-check

??? question "What is wrong with timing a GPU forward pass with time.time()?"
    GPUs execute asynchronously; without synchronisation you time only the submission. Call torch.cuda.synchronize() first and exclude warm-up iterations.

??? question "BF16 or FP16?"
    BF16 has FP32's dynamic range, trains stably and needs no loss scaling — prefer it when hardware supports it; FP16 has more mantissa bits but overflows easily and needs loss scaling.

??? question "Pruning 50% of weights unstructured — why might inference not get faster?"
    Sparse weights are still stored and computed as dense matrices unless the hardware or runtime supports a specific sparse format; structured pruning reduces compute directly.

## Related topics

<div class="grid cards" markdown>

-   __Inference & deployment (Foundation Models)__

    ---

    LLM quantization, KV cache and serving systems.

    [:octicons-arrow-right-24: Open](../../FOUNDATION-MODELS/INFERENCE/inference.md)

-   __Software systems & simulation (Embodied)__

    ---

    Timing and execution in the robot loop.

    [:octicons-arrow-right-24: Open](../../EMBODIED-INTELLIGENCE/SYSTEMS/systems.md)

</div>
