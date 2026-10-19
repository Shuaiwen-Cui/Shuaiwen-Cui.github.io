# 高效训练、推理与 AI 工程

!!! abstract "这一页解决什么"
    把算法变成可重复、可测量、能部署的系统。这是你的主场：嵌入式和边缘部署经验可以直接迁移到资源预算、性能分析和量化；GPU 训练和推理运行时按实际缺口补齐。**P0 工程与复现必修，P1 压缩与性能**。

## 知识地图

| 层次 | 主题 | 需要理解什么 |
|---|---|---|
| P0 | 工程与复现 | 数据版本、随机种子、配置、依赖、checkpoint、实验记录 |
| P1 | 性能与训练效率 | profiling、混合精度、batch、梯度累积、数据加载与内存 |
| P1 | 模型压缩 | 量化、蒸馏、剪枝与稀疏；质量与实际速度分开验证 |
| P2 | 大模型系统 | 高效注意力、KV cache、批处理、并行与通信 |
| P2 | 部署与生命周期 | ONNX 与运行时、边缘设备、监控、漂移、回归评估 |

## 核心机制

**可复现**：固定随机种子、记录代码版本、依赖、配置和数据版本。没有这些，"改了一个超参效果变好"无法确认。

**性能分析**：先测再优化。区分数据加载、前向、反向和优化器步；GPU 计时要同步，前几次迭代是预热不计入。很多训练其实瓶颈在数据加载。

**混合精度**：前向和反向用 FP16 / BF16，主权重保留 FP32。BF16 指数位与 FP32 相同，不容易溢出，现代 GPU 首选；FP16 需要损失缩放。

**量化**：把权重（和激活）从浮点换成 INT8 或更低位宽，模型变小、带宽需求降低。训练后量化（PTQ）用校准数据定缩放因子；量化感知训练（QAT）在训练中模拟量化误差。逐通道量化比逐张量更准。

**蒸馏**：让小模型模仿大模型的软输出（概率分布），比只学硬标签信息更多。

**剪枝**：去掉不重要的权重。非结构化剪枝压缩率高但硬件难加速；结构化剪枝（去掉通道、注意力头）才能真正提速。

**计算受限 vs. 内存受限**：判断一个算子快慢要看算术强度（每读一字节做多少运算）。很多层受内存带宽限制，FLOPs 降低不一定带来时延下降——这和你在 MCU 上优化算子是同一个道理。

**部署**：导出 ONNX，再交给 TensorRT、ONNX Runtime 等运行时做算子融合和量化；端侧还要看内存峰值、算子支持和功耗。

## 精选资源

<div class="grid cards" markdown>

-   :material-chip:{ .lg .middle } __方向课程 · MIT 6.5940 TinyML__

    ---

    高效深度学习：量化、剪枝、蒸馏、神经架构搜索与硬件部署。

    [:octicons-arrow-right-24: 打开课程主页](https://hanlab.mit.edu/courses/2024-fall-65940)

-   :material-timer-outline:{ .lg .middle } __测量 · PyTorch Profiler__

    ---

    在真实模型上定位计算瓶颈、内存开销与算子时延分布。

    [:octicons-arrow-right-24: 打开指南](https://docs.pytorch.org/tutorials/recipes/recipes/profiler_recipe.html)

-   :material-speedometer:{ .lg .middle } __量化实践 · ONNX Runtime__

    ---

    动态与静态量化、校准方法，以及量化后的精度调试。

    [:octicons-arrow-right-24: 打开文档](https://onnxruntime.ai/docs/performance/model-optimizations/quantization.html)

</div>

## 跟随实验 · AI-6 · 一个可解释的性能实验

1. 复用 F1、AI-2 或 F4 中已跑通的模型，不新增任务。
2. 固定设备、输入与精度；区分加载、预热、前向和数据处理。
3. 按设备要求同步计时，重复测量，记录时延分布与峰值内存。
4. 只改 batch、精度或一个优化，重新检查质量与实际资源。

**验收**：报告环境、测量范围、重复次数与质量限制；理论 FLOPs 或参数量不能直接当成部署时延。

## 自测

??? question "GPU 上用 time.time() 测前向时延有什么问题？"
    GPU 是异步执行的，不同步就计时只测到了"提交任务"的时间；要先 torch.cuda.synchronize()，并排除预热迭代。

??? question "BF16 和 FP16 怎么选？"
    BF16 动态范围与 FP32 相同，训练稳定、不需要损失缩放，硬件支持时优先；FP16 精度位多但容易溢出，需要损失缩放。

??? question "非结构化剪枝掉 50% 权重，为什么推理不一定变快？"
    稀疏的权重仍以稠密矩阵存储和计算，除非硬件或运行时支持特定的稀疏格式；结构化剪枝才能直接减少计算。

## 相关专题

<div class="grid cards" markdown>

-   __推理与部署（基础模型）__

    ---

    LLM 的量化、KV cache 与服务系统。

    [:octicons-arrow-right-24: 进入](../../FOUNDATION-MODELS/INFERENCE/inference.md)

-   __软件系统与仿真（具身）__

    ---

    机器人回路中的时序与执行。

    [:octicons-arrow-right-24: 进入](../../EMBODIED-INTELLIGENCE/SYSTEMS/systems.md)

</div>
