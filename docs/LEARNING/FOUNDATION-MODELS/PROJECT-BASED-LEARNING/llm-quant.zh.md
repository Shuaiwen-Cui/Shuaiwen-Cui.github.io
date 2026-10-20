# 实战 · 小语言模型量化：从 RTN 到 GPTQ / AWQ

!!! abstract "这个项目解决什么"
    你已经会在 MCU 上做 FP8 / int8 量化；这个项目把同一套能力搬到大模型上：在一个 0.5B 小模型上，**亲手比较几种量化算法的精度、显存和速度，并解释掉点的来源**。它补的是"只会调 llama.cpp 参数"和"懂量化算法"之间的那一段。一张 8 GB 以上显存的 NVIDIA 显卡即可，约 5–7 天。

## 和 F5 的区别

| | F5 · 推理对比 | 本项目 |
|---|---|---|
| 关注点 | 运行时：不同量化版本跑起来多快 | 算法：每种量化方法为什么好或差 |
| 量化从哪来 | 用别人转好的 GGUF | 自己做校准、自己量化 |
| 交付 | 速度与内存表 | 精度–显存–速度对比 + 逐层误差分析 |

## 学到什么程度

| 能力 | 做完应该能 |
|---|---|
| 量化基础 | 讲清对称 / 非对称、逐张量 / 逐通道 / 分组、W8A8 与 W4A16 的区别 |
| 训练后量化（PTQ） | 用 RTN、GPTQ、AWQ、SmoothQuant 各量化一次，说出各自解决什么问题 |
| 误差分析 | 找出量化最敏感的层，解释激活离群值为什么让激活量化更难 |
| 推理侧 | 测预填充（prefill）与逐 token 解码（decode）的速度，知道 KV Cache 量化省的是什么 |
| 量化感知训练（QAT，选做） | 在一个小任务上比较 PTQ 与 QAT 的差距 |

## 步骤

1. **基线**：Qwen2.5-0.5B-Instruct，FP16。固定评测：WikiText-2 困惑度 + 一个小任务（如 lm-evaluation-harness 里的 ARC-Easy 或 HellaSwag 子集）+ 20 条固定提示的人工抽查。记录显存、prefill 和 decode 速度。
2. **最朴素的量化（RTN）**：自己写几十行代码，把线性层权重按逐通道、按分组（group size 128）量化到 INT8 和 INT4，再反量化回去跑评测。这一步不用任何库，目的是看清误差从哪里来。
3. **GPTQ 与 AWQ**：用 llm-compressor 做 W4A16 的 GPTQ 和 AWQ，校准集用 128 条样本。比较困惑度和任务分数，再把校准集换成不同领域的文本，看结果怎么变。
4. **激活量化与 SmoothQuant**：做 W8A8，先不做平滑、再加 SmoothQuant。画出几层激活的逐通道最大值，直观看到离群通道。
5. **逐层敏感度**：每次只量化一层（其余保持 FP16），记录困惑度变化，画出敏感度图；把最敏感的几层保留为 FP16，再做一次混合精度。
6. **推理测速**：选两种量化结果用 vLLM 部署（有 Jetson 的话再测一次），分别测 prefill、decode 速度和显存；再开 KV Cache FP8 / INT8，看长上下文下省了多少显存。
7. **选做 · QAT**：用 torchao 的 QAT 流程在一个小下游任务上微调，比较 INT4 PTQ 与 QAT 的差距。

## 交付物

- 一张总表：方法 × 位宽 → 困惑度、任务分数、显存、prefill / decode 速度。
- 两张图：激活离群通道图、逐层敏感度图。
- 一段"掉点排查"说明：哪一层最敏感、为什么、怎么修。
- 代码仓库与 README，写明模型版本、校准集、硬件和重复次数；未做完的部分标"计划"。

## 精选资源

<div class="grid cards" markdown>

-   :material-tools:{ .lg .middle } __主工具 · llm-compressor__

    ---

    vLLM 官方的量化库，GPTQ、AWQ、SmoothQuant、FP8 都有现成示例，结果可直接用 vLLM 加载。

    [:octicons-arrow-right-24: 打开仓库](https://github.com/vllm-project/llm-compressor)

-   :material-fire:{ .lg .middle } __QAT 与低比特 · torchao__

    ---

    PyTorch 官方的量化与稀疏库，含 INT4 / INT8 权重量化和 QAT 流程。

    [:octicons-arrow-right-24: 打开仓库](https://github.com/pytorch/ao)

-   :material-clipboard-check:{ .lg .middle } __评测 · lm-evaluation-harness__

    ---

    统一跑困惑度和常见基准，保证量化前后用的是同一把尺子。

    [:octicons-arrow-right-24: 打开仓库](https://github.com/EleutherAI/lm-evaluation-harness)

</div>

**三篇必读论文**（只读方法部分和主要实验图）：

<div class="grid cards" markdown>

-   __GPTQ__

    ---

    用二阶信息逐列量化、补偿误差。

    [:octicons-arrow-right-24: arXiv 2210.17323](https://arxiv.org/abs/2210.17323)

-   __AWQ__

    ---

    按激活大小找出重要权重通道并放大保护。

    [:octicons-arrow-right-24: arXiv 2306.00978](https://arxiv.org/abs/2306.00978)

-   __SmoothQuant__

    ---

    把激活的离群值"转移"到权重上，使 W8A8 可行。

    [:octicons-arrow-right-24: arXiv 2211.10438](https://arxiv.org/abs/2211.10438)

</div>

## 自测

??? question "为什么大模型常用仅权重量化（W4A16），而不是权重和激活都量化？"
    LLM 的激活里有少数通道数值特别大（离群值），按张量量化会把其他通道的精度挤没；而 decode 阶段主要受显存带宽限制，只压缩权重就能拿到大部分速度收益。

??? question "GPTQ 和 RTN 的本质区别？"
    RTN 每个权重独立四舍五入；GPTQ 量化一列后，用 Hessian 近似把这一列的误差分摊到还没量化的列上，整体输出误差更小。

??? question "校准集换成另一个领域的文本，为什么结果会变？"
    GPTQ / AWQ / SmoothQuant 都依赖校准数据统计激活分布；分布和实际使用不一致时，挑出的"重要通道"和缩放因子就不准。

??? question "KV Cache 量化主要省的是什么？"
    长上下文和大批量下 KV Cache 会超过权重本身的显存占用；量化它主要省显存、放大可用上下文和批量，对短提示的单条速度影响不大。

## 相关

<div class="grid cards" markdown>

-   __推理与部署__

    ---

    量化、KV Cache 与服务系统的原理。

    [:octicons-arrow-right-24: 进入](../INFERENCE/inference.md)

-   __F5 · 推理对比__

    ---

    运行时层面的对比，可作为本项目第 6 步的延伸。

    [:octicons-arrow-right-24: 进入](project-based-learning.md#f5)

-   __AI · 高效训练与工程__

    ---

    测速方法、计算受限与内存受限。

    [:octicons-arrow-right-24: 进入](../../AI/SYSTEMS/systems.md)

</div>
