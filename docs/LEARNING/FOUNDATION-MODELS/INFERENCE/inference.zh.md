# 推理与部署

!!! abstract "这一页解决什么"
    模型训练好之后，怎样在给定硬件上跑得快、占用少、质量不掉太多。这是你的主场：MCU 上的量化和内存优化经验可以直接迁移到 LLM，只是规模不同。**必修基础，深度按方向**；前置：能加载并运行一个预训练模型。

## 学到什么程度

| 层次 | 学什么 | 完成标准 |
|---|---|---|
| 扫盲 | 预填充与解码、首 token 时延、吞吐、显存 | 解释首 token 时延与生成速度为什么不同 |
| 基础 | 固定条件下的性能测量、量化对比 | 完成 F5，质量、时延、内存一起报告 |
| 进阶 | KV cache 管理、连续批处理、投机解码、内核优化 | 用测量定位瓶颈并解释收益与限制 |

## 核心知识

**两个阶段**：预填充（prefill）一次性并行处理整段提示，主要受算力限制，决定首 token 时延（TTFT）；解码（decode）逐 token 生成，每步都要读一遍全部权重，主要受内存带宽限制，决定生成速度（tokens/s）。

**粗估解码速度**：单批次时，tokens/s 上限约等于"内存带宽 ÷ 模型权重字节数"。7B 模型 FP16 约 14 GB，在 100 GB/s 的设备上上限约 7 tokens/s；量化到 4-bit 约 3.5 GB，上限约 28 tokens/s。这就是量化提速的根本原因。

**KV cache**：缓存历史 token 的 K、V，避免每步重算；大小 = 2 × 层数 × KV 头数 × 头维度 × 序列长度 × 批大小 × 字节数。长上下文和高并发时它比权重还大。GQA / MQA 通过共享 KV 头缩小它，PagedAttention 解决碎片化。

**解码策略**：贪心、温度采样、top-k、top-p（核采样）、束搜索。温度越高越多样也越容易出错；事实性任务用低温，创作用高温。

**量化**：仅权重量化（GPTQ、AWQ、GGUF 的 k-quant）最常用，4-bit 通常质量损失可接受；权重加激活量化（如 W8A8）更难，因为激活有离群值，SmoothQuant 把难度从激活迁移到权重。

**蒸馏与剪枝**：蒸馏用大模型的输出训练小模型；剪枝去掉不重要的权重或结构（结构化剪枝才能真正加速）。

**投机解码**：小草稿模型先连猜几个 token，大模型一次并行验证，猜对就一次接受多个；输出分布与大模型完全一致，加速取决于接受率。

**连续批处理与服务**：服务端把不同请求在 token 级别动态拼批，充分利用 GPU；vLLM、TensorRT-LLM、SGLang 是代表。

**端侧运行时**：llama.cpp（CPU / GPU / Apple Silicon，GGUF 格式）、ONNX Runtime、TensorRT、MLC-LLM。端侧更关心内存峰值、功耗和首 token 时延。

## 路线之争 · 通用加速 vs. 专用硬件优化

通用派依赖成熟的框架和 GPU 生态，一次开发到处部署；专用派针对 NPU、特定芯片做算子融合和定制量化，在功耗和成本上赢得数量级优势，但移植成本高。端侧和车载场景越来越多走专用路线，这也是嵌入式背景的价值所在。

## 精选资源

<div class="grid cards" markdown>

-   :material-github:{ .lg .middle } __项目主入口 · llama.cpp__

    ---

    在 PC 和 Jetson 上跑小模型，做同一模型不同量化的对比（F5）。

    [:octicons-arrow-right-24: 打开仓库](https://github.com/ggml-org/llama.cpp)

-   :material-school:{ .lg .middle } __方向课 · MIT 6.5940 TinyML__

    ---

    量化、剪枝、蒸馏、高效推理的系统课程；确定高效推理方向后学。

    [:octicons-arrow-right-24: 打开课程](https://hanlab.mit.edu/courses/2024-fall-65940)

</div>

## 跟随实验

<div class="grid cards" markdown>

-   :material-hammer-wrench:{ .lg .middle } __F5 · 固定任务下的推理对比__

    ---

    同一模型的不同量化版本，记录首 token 时延、tokens/s、峰值内存和质量退化样例。

    [:octicons-arrow-right-24: 查看步骤](../PROJECT-BASED-LEARNING/project-based-learning.md)

</div>

## 自测

??? question "为什么批大小增加时，解码的 tokens/s 总量会上升？"
    单批次解码受内存带宽限制，读一次权重只算一个 token；批量越大，读一次权重能服务更多 token，算力利用率上升，直到变成算力受限。

??? question "一个 7B 模型 4-bit 量化后，在 50 GB/s 带宽的设备上，单批次解码速度大约多少？"
    权重约 3.5 GB，上限约 50 / 3.5 ≈ 14 tokens/s，实际会更低一些。

??? question "投机解码会改变输出质量吗？"
    不会。拒绝采样保证最终输出与大模型单独生成的分布一致，只影响速度。

??? question "为什么激活量化比权重量化难？"
    激活的分布随输入变化，且某些通道有极大的离群值，统一缩放会让大多数数值精度不足。

## 相关专题

<div class="grid cards" markdown>

-   __高效训练与工程（AI）__

    ---

    量化、蒸馏、剪枝在视觉与通用模型上的做法。

    [:octicons-arrow-right-24: 进入](../../AI/SYSTEMS/systems.md)

-   __仿真与虚实迁移（具身）__

    ---

    机载算力与实时性约束。

    [:octicons-arrow-right-24: 进入](../../EMBODIED-INTELLIGENCE/SIM-TO-REAL/sim-to-real.md)

</div>
