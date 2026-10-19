# AI · 两周速成

!!! abstract "这一页和系统学习的分工"
    系统学习（全貌、路线、各专题）解决"知识是否完整"；这一页只解决一件事：**两周内把深度学习的工程能力补到工业界算法岗的常见要求，并留下一个能展示的作品**。每天约 2–3 小时。做完再回到系统学习按兴趣深入。

## 目标能力（从常见岗位要求倒推）

深度学习算法、端侧 AI、模型优化类岗位几乎都会考：

| 能力 | 常见说法 | 你已有的基础 | 两周内要补到 |
|---|---|---|---|
| 训练与调试 | PyTorch、训练循环、过拟合、调参 | 论文项目里训练过 Transformer / Mamba | 能脱稿写训练循环，能解释每个超参的作用 |
| 主流架构 | CNN、Transformer、注意力 | 了解原理（checklist 已学） | 能讲清残差、归一化、注意力复杂度 |
| 评估 | 指标、数据泄漏、验证集 | 时间序列划分、AUROC | 补检测与分割指标（mAP、IoU） |
| 部署优化 | ONNX、量化、TensorRT、时延 | **TinyAuton：MCU 上的算子与 FP8 量化** | 补 GPU / Jetson 侧的 INT8 与测速方法 |

**你的差异化**：大多数候选人只会在 GPU 上训练；你能把同一个模型一路部署到 MCU。这两周的作品就围绕这一点。

## 14 天计划

| 天 | 学什么 | 动手 | 当天产出 |
|---|---|---|---|
| 1 | 训练循环复习 | 不看教程手写 PyTorch Quickstart 的训练与验证 | 一份能跑的脚本 + 自己写的注释 |
| 2 | 数据划分、过拟合、指标 | 故意制造过拟合与数据泄漏各一次 | 两张对比曲线 |
| 3 | CNN 与迁移学习 | 在小数据集上微调 ResNet18（冻结 vs 全微调） | 准确率与训练时间对比表 |
| 4 | Transformer 机制 | 读本站 Transformer；看 Karpathy "Let's build GPT" 前半 | 手画一次注意力的张量形状 |
| 5 | 训练技巧 | 学习率预热、归一化、混合精度各试一次 | 记录哪一项影响最大 |
| 6–7 | 检测 / 分割 | 用 Ultralytics 在一个小数据集上微调 YOLO11n | mAP、PR 曲线；能解释 IoU 与 mAP |
| 8–9 | 量化原理 | ONNX Runtime 量化：导出 ONNX，做静态 INT8 | FP32 vs INT8 的精度与 CPU 时延 |
| 10–11 | GPU / Jetson 部署 | Ultralytics Jetson 指南：TensorRT FP16 与 INT8 | 同一模型在 PC GPU 与 Jetson 上的时延表 |
| 12 | 接上自己的强项 | 用 TinyAuton 在 ESP32-S3 上跑一个同任务的小网络 | 三层设备（GPU / Jetson / MCU）的精度–时延–内存对比 |
| 13 | 整理作品 | README：问题、方法、结果表、局限 | GitHub 仓库首页 |
| 14 | 面试复盘 | 逐条回答下方问题，答不出的回到对应天 | 自己的答题笔记 |

## 这两周用到的资源

<div class="grid cards" markdown>

-   :material-open-in-new:{ .lg .middle } __PyTorch Quickstart__

    ---

    第 1 天 · 训练循环复习

    [:octicons-arrow-right-24: 打开](https://docs.pytorch.org/tutorials/beginner/basics/quickstart_tutorial.html)

-   :material-book-open-variant:{ .lg .middle } __本站 · Transformer__

    ---

    第 4 天 · Transformer 机制

    [:octicons-arrow-right-24: 打开](../TRANSFORMER/transformer.md)

-   :material-open-in-new:{ .lg .middle } __Karpathy · Zero to Hero__

    ---

    第 4 天 · Transformer 机制

    [:octicons-arrow-right-24: 打开](https://karpathy.ai/zero-to-hero.html)

-   :material-open-in-new:{ .lg .middle } __ONNX Runtime 量化__

    ---

    第 8–9 天 · 量化原理

    [:octicons-arrow-right-24: 打开](https://onnxruntime.ai/docs/performance/model-optimizations/quantization.html)

-   :material-open-in-new:{ .lg .middle } __Ultralytics Jetson 指南__

    ---

    第 10–11 天 · GPU / Jetson 部署

    [:octicons-arrow-right-24: 打开](https://docs.ultralytics.com/guides/nvidia-jetson/)

</div>


## 作品：同一模型，从 GPU 到 MCU

一张表说清楚：同一个任务，模型在 GPU、Jetson、MCU 上分别是什么精度、多快、占多少内存，量化在每一层损失了多少、为什么。结果按实测写，未做的部分标"计划"。

## 面试必会（答题要点）

| 问题 | 要点 |
|---|---|
| 为什么要残差连接？ | 梯度可以绕过非线性层直接回传，深层网络更容易优化；也让每层学"增量" |
| BatchNorm 和 LayerNorm 的区别？ | BN 在批次维度上统计，依赖批大小，推理时用滑动均值；LN 在特征维度上统计，适合序列和小批量，Transformer 用 LN |
| 自注意力的复杂度？怎么降？ | 序列长度 n 时时间和内存都是 O(n²)；稀疏 / 线性注意力、FlashAttention（不改复杂度，减少显存读写）、状态空间模型（Mamba，线性） |
| Adam 和 SGD 怎么选？ | Adam 收敛快、对学习率不敏感；SGD + 动量在视觉任务上常泛化更好；大模型常用 AdamW（解耦权重衰减） |
| 什么是数据泄漏？ | 测试信息进入训练或选模过程，例如随机划分时间序列、在全数据上做归一化；后果是指标虚高 |
| mAP 怎么算？ | 按置信度排序得到 PR 曲线，对曲线积分得到 AP，再对类别（和 IoU 阈值）求平均 |
| PTQ 和 QAT 的区别？ | PTQ：训练后用少量校准数据定缩放因子，快但可能掉点；QAT：训练时模拟量化误差，精度更好但要重训 |
| INT8 掉点怎么排查？ | 先逐层比较量化前后的输出，找敏感层（常见是首尾层和有离群值的层）；用逐通道量化、改校准方法，或保留个别层为 FP16 |
| FLOPs 低为什么不一定快？ | 很多层受内存带宽限制而不是算力限制；算子是否被硬件和推理引擎融合、是否支持也决定实际时延 |
| 端侧部署要先看什么？ | 内存（权重 + 激活峰值）、目标时延、算子支持、功耗；再选模型尺寸和量化方案 |

## 做完之后

回到 AI 学习路线 的第三轮，按兴趣选一个方向深入；部署方向继续看 高效训练与工程。


<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __全貌（系统学习总览）__

    ---

    系统学习入口

    [:octicons-arrow-right-24: 打开](../ai.md)

-   :material-book-open-variant:{ .lg .middle } __AI 学习路线__

    ---

    系统学习入口

    [:octicons-arrow-right-24: 打开](../ROADMAP/roadmap.md)

-   :material-book-open-variant:{ .lg .middle } __高效训练与工程__

    ---

    系统学习入口

    [:octicons-arrow-right-24: 打开](../SYSTEMS/systems.md)

</div>
