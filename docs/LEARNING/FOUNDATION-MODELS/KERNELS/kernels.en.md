# GPU Programming & Kernel Development

!!! abstract "What this page solves"
    Inference speed ends up in kernels: the same attention or matrix multiply can differ several- to tens-fold depending on how it is written. This page covers how GPUs execute, how to diagnose bottlenecks, how classic kernels are optimised, and where Triton, compilers and profilers fit. Your experience hand-writing DSP / NN kernels and optimising memory on MCUs maps directly onto the central lesson here: memory traffic costs more than arithmetic.

## Target level

| Level | Learn | Done when |
|---|---|---|
| Literacy | GPU execution model, memory hierarchy, compute- vs memory-bound | You can use a roofline to explain why a kernel is slow |
| Foundation | Reduction, softmax, LayerNorm, GEMM tiling; kernel fusion | You write a fused kernel in Triton and benchmark it against PyTorch |
| Advanced | FlashAttention, Tensor Cores, quantized GEMM, compilers | You can explain FlashAttention's tiling and online softmax and why it is fast |

## Core knowledge

**Execution model**: a kernel launches as a grid of thread blocks; blocks are scheduled onto SMs; an SM executes threads in groups of 32 (warps) in lockstep. Divergent branches within a warp execute serially and waste compute.

**Memory hierarchy**: registers → shared memory (per block, tens to a few hundred KB, near-register speed) → L2 → HBM (large, high bandwidth, high latency). The core of optimisation is **read data from HBM once and reuse it several times in shared memory and registers**.

**Memory access**: neighbouring threads must access neighbouring addresses for coalesced loads; shared memory is organised in banks and simultaneous hits on one bank conflict; occupancy, set by register and shared-memory use, hides memory latency.

**Roofline and arithmetic intensity**: intensity = FLOPs ÷ bytes moved. Below the hardware's compute-to-bandwidth ridge a kernel is memory-bound and speeds up by moving fewer bytes (fusion, quantisation, reuse); above it, it is compute-bound and speeds up with Tensor Cores and better tiling. Matrix multiplies in LLM decoding are essentially memory-bound.

**Classic kernels**:

- **Reduction**: tree reduction plus warp shuffles, avoiding shared-memory bank conflicts.
- **Softmax**: max then sum needs two passes; **online softmax** updates the running max and normaliser in one pass.
- **LayerNorm / RMSNorm**: one block per row, mean and variance in a single pass (Welford), often fused with the residual add.
- **GEMM**: tile matrices into shared memory, each thread computes a small output tile, then hand the inner product to Tensor Cores (WMMA / MMA) in mixed precision.

**Kernel fusion**: fold element-wise operations (activation, bias, residual, quantise / dequantise) into neighbouring kernels to avoid writing intermediates to HBM. It is the first lever of every inference engine.

**FlashAttention**: standard attention writes the N×N score matrix to memory; FlashAttention streams Q, K, V tiles into shared memory and accumulates the output with online softmax, never materialising the intermediate. Same FLOPs, far less memory traffic, so faster and leaner.

**Triton vs CUDA**: CUDA controls individual threads, offering the highest ceiling at the highest cost; Triton is written per block in Python and handles coalescing and shared memory automatically. Many FlashAttention and vLLM kernels have Triton versions.

**Compilers**: `torch.compile` (TorchInductor generates Triton), TensorRT (graph optimisation, kernel selection, quantisation), TVM / MLIR (multi-target compiler stacks). Edge NPU vendor toolchains sit at this layer too.

**Profiling**: Nsight Systems for the whole timeline (CPU launch, copies, gaps between kernels); Nsight Compute for one kernel's bandwidth, occupancy and bottleneck; PyTorch Profiler to locate the operator. Measure before you change.

## Debate · hand-written kernels vs compiler-generated

Hand-tuners (CUDA / CUTLASS) chase peak performance, and hot kernels are still expert-tuned; compilers (Triton, TorchInductor, TVM) cover the long tail with auto-tuning and code generation and iterate faster. In practice both coexist: hot paths by hand, the rest by compiler.

## Curated resources

<div class="grid cards" markdown>

-   :material-language-python:{ .lg .middle } __Main practice · Triton tutorials__

    ---

    Vector add → fused softmax → matmul → LayerNorm → FlashAttention, each with benchmarks.

    [:octicons-arrow-right-24: Open tutorials](https://triton-lang.org/main/getting-started/tutorials/index.html)

-   :material-video:{ .lg .middle } __Course · GPU MODE lectures__

    ---

    Community lectures and code: CUDA basics, profiling, Triton, FlashAttention, quantized kernels.

    [:octicons-arrow-right-24: Open repository](https://github.com/gpu-mode/lectures)

-   :material-book-open-variant:{ .lg .middle } __Reference · CUDA C++ Programming Guide__

    ---

    The official account of the execution model, memory hierarchy and best practices.

    [:octicons-arrow-right-24: Open docs](https://docs.nvidia.com/cuda/cuda-c-programming-guide/)

</div>

## Follow-along · write a fused kernel

1. Run Triton's fused-softmax tutorial and confirm why it beats PyTorch's two-pass version.
2. Write a fused "residual add + RMSNorm" kernel covering several row lengths.
3. Compare bandwidth (GB/s) against native PyTorch with `triton.testing`, plotted against row length.
4. Inspect one kernel in Nsight Compute; record bandwidth utilisation and occupancy and explain the bottleneck.

**Done when**: correctness test (error vs PyTorch), performance curve and bottleneck explanation; a cloud notebook works without a local GPU.

## Self-check

??? question "Why does quantizing weights directly speed up LLM decoding?"
    Each generated token reads all weights from memory once; arithmetic intensity is low, so decoding is memory-bound. Going from 16-bit to 4-bit weights cuts bytes read about fourfold and raises the speed ceiling accordingly.

??? question "Does FlashAttention reduce FLOPs?"
    No, compute stays O(N²). It reduces memory traffic: the score matrix never goes back to HBM and tiles are processed in shared memory, so it is faster and uses less memory.

??? question "How does online softmax read the data only once?"
    Keep a running max m and normaliser l; when a larger value appears, rescale the existing sum by exp(m_old − m_new) and continue. Normalise by l at the end.

??? question "What is warp divergence and how do you avoid it?"
    Threads in one warp take different branches, so the hardware runs the branches one after another while others idle. Keep a warp on one path, or replace branches with data reordering and masking.

??? question "An element-wise kernel reaches 30% of peak bandwidth. What do you check first?"
    Whether loads are coalesced, whether each thread does too little work (launch overhead dominates), and whether there are redundant intermediate reads and writes; then whether it can be fused with neighbours.

## Related topics

<div class="grid cards" markdown>

-   __Inference & Deployment__

    ---

    Prefill and decode, KV cache, quantisation and serving.

    [:octicons-arrow-right-24: Open](../INFERENCE/inference.md)

-   __Project · LLM Quantization__

    ---

    Quantized GEMM and the accuracy–speed trade-off.

    [:octicons-arrow-right-24: Open](../PROJECT-BASED-LEARNING/llm-quant.md)

-   __AI · Efficient Training & Engineering__

    ---

    Measurement methods; compute- vs memory-bound.

    [:octicons-arrow-right-24: Open](../../AI/SYSTEMS/systems.md)

</div>
