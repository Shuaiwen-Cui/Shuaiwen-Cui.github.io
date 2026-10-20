# 基础模型 · 快速入门

!!! abstract "这一页和系统学习的分工"
    系统学习（全貌、路线、F0–F5）按"扫盲 → 实现 → 使用 → 适配"走完整条线；这一页只做<strong>两周内能上手的最小闭环</strong>：读懂一个小 GPT、用 LoRA 微调一个小模型、做一个带来源的检索问答、把它量化后离线跑在边缘设备上。每天约 2–3 小时。

## 目标能力（从常见岗位要求倒推）

LLM 应用、端侧大模型、推理优化类岗位最常出现的要求：

| 能力 | 常见说法 | 你已有的基础 | 两周内要补到 |
|---|---|---|---|
| 模型机制 | Transformer、注意力、KV cache | 训练过 Transformer / Mamba 时序模型 | 能讲清解码器、预填充和逐 token 解码 |
| 使用与适配 | Hugging Face、SFT、LoRA | PyTorch 训练经验 | 独立完成一次 LoRA 微调并评估 |
| 应用 | RAG、评测、幻觉 | 你的网站就是现成语料 | 一个能显示来源、有固定评测集的问答 |
| 推理部署 | 量化、llama.cpp、吞吐与时延 | **MCU 上的量化与内存优化** | 4-bit 量化、首 token 时延、tokens/s 的测法 |

## 14 天计划

| 天 | 学什么 | 动手 | 当天产出 |
|---|---|---|---|
| 1 | 全貌 | HF LLM Course 第 1 章；本站 Overview | 自己画的 LLM 系统图（训练、推理、检索、工具） |
| 2–4 | 小 GPT | Karpathy "Let's build GPT" 全程跟写 | 训练曲线；每个张量的形状注释 |
| 5 | 使用现成模型 | 用 transformers 加载 Qwen2.5-0.5B-Instruct，试对话模板与采样参数 | temperature / top-p 对比样例 |
| 6–7 | LoRA 微调 | PEFT：用 100–200 条你网站内容改写的问答对微调 | 微调前后在固定 20 题上的对比 |
| 8–9 | RAG | 把网站 Markdown 分块、向量化、检索（可用 LlamaIndex，见 F2） | 带来源显示的问答；20 题评测，区分检索失败与生成失败 |
| 10–11 | 量化与推理 | llama.cpp：转 GGUF，比较 Q8 与 Q4 量化；PC 和 Jetson 各测一次 | 首 token 时延、tokens/s、峰值内存、质量退化样例 |
| 12 | 全离线闭环 | 让 RAG 调用本地量化模型 | Jetson 上离线运行的个人知识库问答 |
| 13 | 整理作品 | README：架构图、结果表、局限 | GitHub 仓库首页 |
| 14 | 面试复盘 | 逐条回答下方问题 | 答题笔记 |

## 这两周用到的资源

<div class="grid cards" markdown>

-   :material-open-in-new:{ .lg .middle } __HF LLM Course 第 1 章__

    ---

    第 1 天 · 全貌

    [:octicons-arrow-right-24: 打开](https://huggingface.co/learn/llm-course/en/chapter1/1)

-   :material-book-open-variant:{ .lg .middle } __Overview__

    ---

    第 1 天 · 全貌

    [:octicons-arrow-right-24: 打开](../OVERVIEW/overview.md)

-   :material-open-in-new:{ .lg .middle } __Karpathy · Zero to Hero__

    ---

    第 2–4 天 · 小 GPT

    [:octicons-arrow-right-24: 打开](https://karpathy.ai/zero-to-hero.html)

-   :material-open-in-new:{ .lg .middle } __Qwen2.5-0.5B-Instruct__

    ---

    第 5 天 · 使用现成模型

    [:octicons-arrow-right-24: 打开](https://huggingface.co/Qwen/Qwen2.5-0.5B-Instruct)

-   :material-open-in-new:{ .lg .middle } __Hugging Face PEFT__

    ---

    第 6–7 天 · LoRA 微调

    [:octicons-arrow-right-24: 打开](https://huggingface.co/docs/peft/index)

-   :material-book-open-variant:{ .lg .middle } __F2 · 检索问答项目__

    ---

    第 8–9 天 · RAG

    [:octicons-arrow-right-24: 打开](../PROJECT-BASED-LEARNING/project-based-learning.md)

</div>


## 作品：跑在边缘设备上的个人知识库问答

LoRA 微调的小模型 + 检索增强 + 4-bit 量化，在 Jetson 上离线运行，回答关于你研究和项目的问题并给出来源。它同时覆盖"适配、应用、部署"三类要求，也能顺带展示你的网站内容。

## 面试必会（答题要点）

| 问题 | 要点 |
|---|---|
| Q、K、V 各是什么？ | 每个 token 通过三个线性映射得到；Q 与 K 的点积决定关注谁，权重再加权 V；除以 √d 防止 softmax 饱和 |
| KV cache 为什么能加速？ | 解码时历史 token 的 K、V 不变，缓存后每步只算新 token；代价是显存随序列长度线性增长 |
| 预填充和解码有什么不同？ | 预填充并行处理整段提示，偏算力受限，决定首 token 时延；解码逐 token 生成，偏内存带宽受限，决定 tokens/s |
| RoPE 是什么？ | 旋转位置编码：按位置旋转 Q、K，使点积只依赖相对位置，外推性较好 |
| 预训练、SFT、偏好优化分别改变什么？ | 预训练学语言与知识；SFT 学指令格式与行为；RLHF / DPO 按人类偏好调整输出倾向 |
| LoRA 的原理？ | 冻结原权重 W，只学低秩增量 ΔW = BA（秩 r 远小于维度），可训练参数少、可合并回原权重；alpha / r 控制缩放 |
| 常见量化方法？ | 仅权重量化（GPTQ、AWQ、GGUF 的 k-quant）为主；激活量化更难，因为有离群值；4-bit 通常质量损失可接受 |
| RAG 失败怎么定位？ | 先看检索到的片段是否包含答案：没有是检索问题（分块、嵌入、召回数），有却答错是生成问题（提示、模型能力） |
| 怎么评估幻觉？ | 固定问题集，含无答案题；检查回答是否有来源支撑，统计无依据断言的比例 |
| 上下文越长越好吗？ | 成本和时延随长度增长，注意力在长上下文中会分散；检索 + 精选上下文通常比塞满更好 |

## 做完之后

回到 基础模型路线 的"根据实践再选深入方向"；推理方向继续看 推理与部署，机器人方向接 具身智能速成。


<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __全貌（系统学习总览）__

    ---

    系统学习入口

    [:octicons-arrow-right-24: 打开](../OVERVIEW/overview.md)

-   :material-book-open-variant:{ .lg .middle } __基础模型路线__

    ---

    系统学习入口

    [:octicons-arrow-right-24: 打开](../ROADMAP/roadmap.md)

-   :material-book-open-variant:{ .lg .middle } __推理与部署__

    ---

    系统学习入口

    [:octicons-arrow-right-24: 打开](../INFERENCE/inference.md)

-   :material-book-open-variant:{ .lg .middle } __具身智能速成__

    ---

    系统学习入口

    [:octicons-arrow-right-24: 打开](../../EMBODIED-INTELLIGENCE/GET-STARTED/get-started.md)

</div>
