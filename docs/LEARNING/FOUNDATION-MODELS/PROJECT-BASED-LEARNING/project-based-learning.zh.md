# 基础模型 · 跟随式项目

!!! abstract "先复现，再修改，再解释"
    先做 F0；F1 按需补齐，F4 的小型 GPT 是实现主项目，任务适配部分以后选做。F2/F3 二选一体验，F5 按方向选修。每个项目先跑教程原例，再做下面的个人实验。这里提供学习与验收步骤，具体安装和命令以匹配版本的上游教程为准。

| 项目 | 层次 | 前置 | 资源与设备 |
|---|---|---|---|
| F0 · 模型概念体验 | 扫盲 | 导论 | 浏览器；本地 tokenizer 可选 |
| F1 · 小模型训练 | 基础 | Python、数组与函数 | 普通电脑 CPU 或在线 Notebook |
| F2 · 个人资料问答 | 应用选做 | Python、模型调用 | 文本资料；按教程选择 API 或本地模型，API 费用取决于服务 |
| F3 · 图片问答 | 体验选做 | 模型加载与推理 | 按所选 VLM 配置内存/GPU，或使用在线示例 |
| F4 · 小型 GPT 与适配 | 实现主线；适配选修 | F1、分词与 Attention | 小规模 CPU 实验起步；较大训练使用 GPU |
| F5 · 推理对比 | 进阶选修 | 能运行一个模型 | 固定模型、测试集与硬件 |
| F6 · 大模型量化实战 | 进阶选修 | F5 或会加载 HF 模型 | 8 GB 以上显存的 NVIDIA 显卡 |

## F0

### token、上下文与模型边界

**跟随：**[HF LLM Course 第一章](https://huggingface.co/learn/llm-course/en/chapter1/1)，再按第二章 tokenizer 内容做实验。

1. 选同一意思的中英文句子，比较分词与 token 数。
2. 对同一问题分别给少量资料、完整资料与无关资料。
3. 加入资料中没有答案的问题，检查是否出现无依据推断。
4. 固定问题、记录上下文与回答，标注模型和日期。

**个人修改：**只改变上下文，观察回答质量变化。

**验收：**一张结果表和系统图；能解释分词、上下文、生成与知识检索的关系。

## F1

### PyTorch 图片分类器

**先自查：**能解释训练/验证、梯度更新并独立保存加载模型时，直接进入 F4；否则只补 Quickstart 暴露的缺口。

**跟随：**[Learn the Basics](https://docs.pytorch.org/tutorials/beginner/basics/intro.html)，使用教程的 FashionMNIST 数据与 Notebook。

1. 跑 Quickstart，理解完整训练流程。
2. 依次学习 Tensor、DataLoader、模型、自动微分与优化。
3. 保存模型，在独立脚本中加载并预测。
4. 改学习率或模型宽度，保留其他条件不变。
5. 从训练数据中另划验证集做选择，最后才使用测试集报告结果；展示预测错误的图片。

**个人修改：**至少比较两种配置，并解释过拟合或训练不足。

**验收：**代码、依赖、训练/验证曲线、测试结果、错误图片与简短分析。F1 同时作为具身策略学习的准备。

## F2

### 个人网站资料问答

**跟随：**[LlamaIndex 官方 Starter Tutorial](https://developers.llamaindex.ai/python/framework/getting_started/starter_example/) / [官方 GitHub](https://github.com/run-llama/llama_index)。只完成最小文档问答流程；默认云端模型需要配置 API 并承担相应费用。

1. 运行原例，记录教程和依赖版本。
2. 选自己的 5–10 篇 Markdown，检查解析与分块结果。
3. 建索引，先查看检索片段，再检查生成回答。
4. 显示来源文件和片段，标记无答案情况。
5. 准备至少 20 个问题，含直接事实、跨段整合和无答案问题；预先记录应有依据。
6. 比较直接问答与 RAG；调整一个分块或检索设置，再用留出的样例评估。

**个人修改：**替换为自己的资料；增加来源展示和固定评测集。

**验收：**能运行的应用、流程图、问题集、来源支持检查与失败案例。至少区分检索失败和生成失败。

**后续：**只在当前代码中改进检索与评估；工具扩展见[应用与智能体](../APPLICATIONS/applications.md)，不另开一套 RAG 课程。

## F3

### 视觉语言模型的能力与限制

**跟随：**[HF Image-text-to-text](https://huggingface.co/docs/transformers/main/tasks/image_text_to_text)。文档为开发版，安装时匹配源码或选择对应稳定版页面。

1. 按教程运行一个图像问答原例。
2. 固定 10 张图片，分别询问物体、数量、空间关系与任务相关信息。
3. 保留提示、输出和人工判断，记录至少三类错误。
4. 改一处提示或输入尺寸，比较效果与推理时间。

**验收：**一份有成功和失败样例的报告，说明看图问答距可靠机器人动作还缺什么。进入动作学习请转[具身 E3–E4](../../EMBODIED-INTELLIGENCE/PROJECT-BASED-LEARNING/project-based-learning.md#e3)。

## F4

### 从小型 GPT 到一次适配

**主教程：**[Karpathy Zero to Hero](https://karpathy.ai/zero-to-hero.html) / [配套代码仓库](https://github.com/karpathy/nn-zero-to-hero)。先 micrograd（梯度基础不足时），再 makemore 前两讲，最后跟 GPT 视频及其代码。用小模型、小数据验证，不要求训练成通用助手。

**文字替代：**更喜欢逐章 Notebook 时，用 [LLMs-from-scratch](https://github.com/rasbt/LLMs-from-scratch) 第 2–5 章替代上述实现。完整叙述见配套书；一次选择一套实现。之后需要任务适配时，再选该书第 6 或 7 章，或 HF 官方微调说明。

**个人修改：**比较上下文长度、训练数据量或模型宽度中的一个因素。预训练和适配可使用不同规模模型，不要求把玩具模型训练成通用聊天助手。

**验收：**模型结构图、训练曲线、独立样例及能力限制；能解释主要 Tensor 的形状。LoRA/SFT 的进一步练习在[训练与适配](../TRAINING/training.md)中按需选择。

## F5

### 固定任务下的推理对比

**跟随：**[llama.cpp](https://github.com/ggml-org/llama.cpp) 的 README 与运行说明，结合[推理专题](../INFERENCE/inference.md)。

先运行受支持的小模型，再比较同一基座模型的不同量化版本；固定硬件、提示集、输出长度和运行设置，记录首 token 时间、生成速度、峰值内存与回答质量。重复测量并区分模型加载和实际推理。

**验收：**配置与结果表、质量退化样例、瓶颈解释；如使用不同模型或设备，明确说明比较限制。


<div class="grid cards" markdown>

-   :material-memory:{ .lg .middle } __F6 · 实战：小语言模型量化__

    ---

    从 RTN 到 GPTQ / AWQ / SmoothQuant：自己做校准和量化，比较精度、显存、速度，并做逐层误差分析。

    [:octicons-arrow-right-24: 进入项目](llm-quant.md)

</div>

## 统一项目记录

保存原教程链接、运行日期、代码 commit、依赖、模型/数据版本、运行方式、结果和个人修改。课程复现、个人实验和独立项目分别标注，不把计划成果写成已经完成。

[返回路线](../ROADMAP/roadmap.md)
