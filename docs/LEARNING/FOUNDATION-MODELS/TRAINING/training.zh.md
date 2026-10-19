# 训练与适配

!!! abstract "这一页解决什么"
    一个大模型从"读了很多书"到"会听指令、符合偏好"经历了哪些阶段，以及你在有限算力下能做什么（通常是适配，而不是从头预训练）。**必修概念，实践做到 LoRA**；前置：F1 的训练循环、数据划分。

## 学到什么程度

| 层次 | 学什么 | 完成标准 |
|---|---|---|
| 扫盲 | 预训练、SFT、偏好优化、参数高效微调 | 解释每个阶段改变了模型的什么 |
| 基础 | 数据划分、超参、checkpoint、LoRA | 完成一次小模型 LoRA 微调，用独立样例评估 |
| 进阶 | 扩展律、并行训练、偏好优化、可验证奖励 RL | 按方向选一个问题深入 |

## 核心知识

**预训练目标**：在海量文本上做下一个 token 预测（交叉熵损失）。模型的知识和大部分能力来自这一阶段，成本也集中在这里。

**扩展律（scaling laws）**：损失随参数量、数据量和算力按幂律下降。Chinchilla 的结论是给定算力时参数和数据应同比增长（约每参数 20 个 token）；现在的小模型常"过度训练"以换取更便宜的推理。

**数据**：去重、过滤低质量、配比（网页、代码、数学、多语言）对效果的影响常大于架构改动。数据泄漏（评测题进了训练集）会让指标失真。

**分布式训练**：数据并行（每卡一份模型、分数据）、张量并行（切分矩阵）、流水线并行（切分层）；ZeRO / FSDP 把优化器状态、梯度和参数切到多卡以省显存。

**指令微调（SFT）**：用"指令-回答"对训练，让模型学会对话格式和听指令。只计算回答部分的损失。

**偏好优化**：RLHF 先训练奖励模型再用 PPO 优化；DPO 直接用"好回答 vs. 差回答"对训练，不需要单独的奖励模型和强化学习循环，更简单稳定。

**可验证奖励的强化学习**：数学、代码这类答案能自动判对错的任务，用正确性当奖励做 RL（如 GRPO），是近期推理模型能力提升的主要来源。

**参数高效微调**：LoRA 冻结原权重 W，只训练低秩增量 ΔW = BA，可训练参数通常不到 1%；QLoRA 把底座量化到 4-bit 再训 LoRA，单卡就能微调数十亿参数的模型。

## 路线之争 · 继续扩展 vs. 方法创新

扩展派认为更多数据和算力仍是最可靠的进步来源；方法派认为高质量数据、合成数据、强化学习和测试时计算（推理时多想几步）能以更低成本获得同等收益。近两年的趋势是两者结合：预训练扩展放缓，后训练和推理时计算的比重上升。

## 精选资源

<div class="grid cards" markdown>

-   :material-code-braces:{ .lg .middle } __基础实现 · Karpathy__

    ---

    沿同一个小模型理解损失、训练 / 验证与采样（F4）。

    [:octicons-arrow-right-24: 打开课程](https://karpathy.ai/zero-to-hero.html)

-   :hugging:{ .lg .middle } __适配 · Hugging Face PEFT__

    ---

    LoRA / QLoRA 的原理与用法，跑一次小模型适配并独立评估。

    [:octicons-arrow-right-24: 打开文档](https://huggingface.co/docs/peft/index)

-   :material-school:{ .lg .middle } __深入（选修）· Stanford CS336__

    ---

    从数据、架构到训练系统从头实现语言模型；确定做模型训练方向后再学。

    [:octicons-arrow-right-24: 打开课程](https://cs336.stanford.edu/)

</div>

## 跟随实验

<div class="grid cards" markdown>

-   :material-hammer-wrench:{ .lg .middle } __F4 · 小型 GPT 与一次适配__

    ---

    先预训练小模型，再在现成小模型上做一次 LoRA 适配并评估。

    [:octicons-arrow-right-24: 查看步骤](../PROJECT-BASED-LEARNING/project-based-learning.md)

</div>

## 自测

??? question "SFT 时为什么只对回答部分计算损失？"
    指令部分是给定的输入，不需要模型学会生成它；对它计算损失会浪费容量，并让模型倾向于复述提示。

??? question "LoRA 的秩 r 设得太小或太大分别会怎样？"
    太小：表达能力不足，学不会任务；太大：参数增加、更易过拟合，收益递减。常见取 8–64，配合 alpha 控制缩放。

??? question "DPO 相比 RLHF 省掉了什么？"
    省掉了单独训练的奖励模型和 PPO 强化学习循环，直接在偏好对上做一个分类式损失，训练更稳定、资源需求更低。

??? question "为什么小模型常被训练远超 Chinchilla 最优的 token 数？"
    Chinchilla 只考虑训练算力最优；模型部署后推理成本按次累积，用更多数据训练一个更小的模型，推理更便宜、总成本更低。

## 相关专题

<div class="grid cards" markdown>

-   __推理与部署__

    ---

    训练好的模型怎样跑得快、跑得省。

    [:octicons-arrow-right-24: 进入](../INFERENCE/inference.md)

-   __评估与可信 AI（AI）__

    ---

    怎样判断适配真的有效。

    [:octicons-arrow-right-24: 进入](../../AI/EVALUATION/evaluation.md)

</div>
