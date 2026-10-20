# GPU 编程与算子开发

!!! abstract "这一页解决什么"
    推理加速最终落在算子上：同一个注意力或矩阵乘，写法不同，速度可以差几倍到几十倍。这一页讲 GPU 怎样执行、性能瓶颈怎样判断、几个经典算子怎样优化，以及 Triton、编译器和性能分析工具的位置。你在 MCU 上手写 DSP / NN 算子、做内存优化的经验，正好对应这里的「访存比计算更贵」。

## 学到什么程度

| 层次 | 学什么 | 完成标准 |
|---|---|---|
| 扫盲 | GPU 执行模型、内存层次、计算受限与访存受限 | 能用 roofline 解释一个算子为什么慢 |
| 基础 | 归约、softmax、LayerNorm、GEMM 分块；算子融合 | 用 Triton 写一个融合算子，并和 PyTorch 对比速度 |
| 进阶 | FlashAttention、Tensor Core、量化 GEMM、编译器 | 读懂 FlashAttention 的分块与在线 softmax，能讲清它为什么快 |

## 核心知识

**执行模型**：一个 kernel 启动成网格（grid），网格由线程块（block）组成，块被调度到 SM 上；SM 以 32 个线程为一组（warp）同步执行同一条指令。warp 内分支不一致（divergence）会串行执行，浪费算力。

**内存层次**：寄存器 → 共享内存（每个 block 可见，几十到上百 KB，接近寄存器速度）→ L2 → 显存 HBM（容量大、带宽高但延迟大）。优化的核心是<strong>让数据从 HBM 只读一次，在共享内存和寄存器里多用几次</strong>。

**访存要点**：相邻线程访问相邻地址才能合并访存（coalescing）；共享内存按 bank 组织，多个线程撞同一 bank 会冲突；占用率（occupancy）由寄存器和共享内存用量决定，用来掩盖访存延迟。

**Roofline 与算术强度**：算术强度 = 浮点运算数 ÷ 访存字节数。低于硬件的「算力 ÷ 带宽」拐点就是访存受限，提速靠减少访存（融合、量化、复用）；高于拐点是计算受限，提速靠 Tensor Core 和更好的分块。LLM 解码阶段的矩阵乘基本是访存受限。

**经典算子**：

- **归约**：树形归约 + warp shuffle，避免共享内存 bank 冲突。
- **Softmax**：先求最大值再求和需要读两遍；**在线 softmax** 边读边更新最大值和归一化因子，只读一遍。
- **LayerNorm / RMSNorm**：一行一个 block，均值方差一次遍历（Welford），常与残差加法融合。
- **GEMM**：把矩阵切成小块载入共享内存，每个线程算一个小块结果；再交给 Tensor Core（WMMA / MMA 指令）做混合精度乘加。

**算子融合**：把逐元素操作（激活、偏置、残差、量化 / 反量化）合进前后算子，省掉中间结果写回 HBM。这是推理引擎提速的第一手段。

**FlashAttention**：标准注意力要把 N×N 的分数矩阵写回显存；FlashAttention 把 Q、K、V 分块读进共享内存，用在线 softmax 逐块累加输出，不落地中间矩阵。计算量不变，访存从 O(N²) 降下来，所以更快、更省显存。

**Triton 与 CUDA**：CUDA 控制到线程级，上限最高也最费工；Triton 以 block 为单位写 Python，自动处理合并访存和共享内存，开发效率高，FlashAttention、vLLM 的很多算子都有 Triton 版本。

**编译器**：`torch.compile`（TorchInductor 生成 Triton 代码）、TensorRT（图优化 + 算子选择 + 量化）、TVM / MLIR（面向多种硬件的编译栈）。端侧 NPU 厂商的工具链也属于这一层。

**性能分析**：Nsight Systems 看整条时间线（CPU 发射、数据拷贝、kernel 间空隙）；Nsight Compute 看单个 kernel 的带宽、占用率、瓶颈；PyTorch Profiler 定位到算子。先测再改。

## 路线之争 · 手写算子 vs. 编译器生成

手写派（CUDA / CUTLASS）追求极致性能，关键算子仍靠专家手调；编译器派（Triton、TorchInductor、TVM）用自动调优和代码生成覆盖长尾算子，迭代快。现实是两者并存：热点算子手写，其余交给编译器。

## 精选资源

<div class="grid cards" markdown>

-   :material-language-python:{ .lg .middle } __主线实践 · Triton 官方教程__

    ---

    向量加法 → 融合 softmax → 矩阵乘 → LayerNorm → FlashAttention，循序渐进，每个都有性能对比。

    [:octicons-arrow-right-24: 打开教程](https://triton-lang.org/main/getting-started/tutorials/index.html)

-   :material-video:{ .lg .middle } __系统课程 · GPU MODE 讲座__

    ---

    社区讲座与配套代码：CUDA 基础、性能分析、Triton、FlashAttention、量化 kernel。

    [:octicons-arrow-right-24: 打开仓库](https://github.com/gpu-mode/lectures)

-   :material-book-open-variant:{ .lg .middle } __查阅 · CUDA C++ 编程指南__

    ---

    执行模型、内存层次与最佳实践的官方说明，遇到具体问题时查。

    [:octicons-arrow-right-24: 打开文档](https://docs.nvidia.com/cuda/cuda-c-programming-guide/)

</div>

## 跟随实验 · 写一个融合算子

1. 跑通 Triton 教程的融合 softmax，确认它比 PyTorch 两次读写版本快的原因。
2. 自己写「残差加法 + RMSNorm」融合算子，覆盖不同行长度。
3. 用 `triton.testing` 对比 PyTorch 原生实现的带宽（GB/s），画出随行长变化的曲线。
4. 用 Nsight Compute 看一次 kernel，记录带宽利用率与占用率，解释瓶颈。

**验收**：正确性测试（与 PyTorch 结果误差）、性能曲线、瓶颈解释；没有 GPU 时可在云端 Notebook 完成。

## 自测

??? question "为什么 LLM 解码阶段量化权重能直接提速？"
    解码时每生成一个 token 都要把全部权重从显存读一遍，算术强度很低，是访存受限；权重从 16 位降到 4 位，读的字节少了约四倍，速度上限随之提高。

??? question "FlashAttention 降低了计算量吗？"
    没有，计算量还是 O(N²)。它减少的是显存读写：中间的分数矩阵不写回 HBM，分块在共享内存里完成，所以更快、更省显存。

??? question "在线 softmax 是怎样做到只读一遍的？"
    遍历时维护当前最大值 m 和归一化和 l；遇到更大的值时，用 exp(m_old − m_new) 修正已有的和，再继续累加。最后用 l 归一化。

??? question "warp divergence 是什么？怎么避免？"
    同一 warp 的线程走了不同分支，硬件只能先后执行各分支，其他线程空等。尽量让同一 warp 走同一路径，或把分支改成按数据重排、掩码计算。

??? question "一个逐元素算子的带宽只有理论值的 30%，先查什么？"
    先查访存是否合并、每个线程处理的数据量是否太少（启动开销占比高）、是否有多余的中间读写；再看能否和前后算子融合。

## 相关专题

<div class="grid cards" markdown>

-   __推理与部署__

    ---

    预填充与解码、KV Cache、量化与服务系统。

    [:octicons-arrow-right-24: 进入](../INFERENCE/inference.md)

-   __实战 · 大模型量化__

    ---

    量化 GEMM 与精度、速度的取舍。

    [:octicons-arrow-right-24: 进入](../PROJECT-BASED-LEARNING/llm-quant.md)

-   __AI · 高效训练与工程__

    ---

    性能测量方法、计算受限与访存受限。

    [:octicons-arrow-right-24: 进入](../../AI/SYSTEMS/systems.md)

</div>
