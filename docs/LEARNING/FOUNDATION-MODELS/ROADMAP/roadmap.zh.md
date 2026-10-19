# 基础模型 · 学习路线

!!! tip "时间紧时"
    先做 [两周速成](../GET-STARTED/get-started.md)：它把阶段 0、2、3 与 F2、F5 压缩成一个端到端作品；之后按本页选择深入方向。

## 两门主课，一条连续实现路线

**中文全貌：李宏毅导论 → 代码理解：Karpathy Zero to Hero。**先建立全貌，再运行和修改模型。时间按每周 8–10 小时估计，以实验结果决定节奏。

| 阶段 | 主资源 | 针对你的学习范围 | 实践与完成标准 |
|---|---|---|---|
| 0 · 全貌，1–2 周 | [李宏毅导论](https://speech.ee.ntu.edu.tw/~hylee/genai/2024-spring.php) | 第 1、6–10、12 讲，先理解系统与能力边界 | F0：解释训练、推理、检索和工具，画出流程 |
| 1 · 基础检查，1–3 次学习时段 | [PyTorch Quickstart](https://docs.pytorch.org/tutorials/beginner/basics/quickstart_tutorial.html) | 先试着解释训练循环、验证、保存加载；有缺口再补 Basics | F1 按需；已经熟练即可通过，不重学 Python 和基础数学 |
| 2 · 实现，2–4 周 | [Karpathy Zero to Hero](https://karpathy.ai/zero-to-hero.html) | 不懂反向传播先看 micrograd；再做 makemore 前两讲，最后 GPT；其余讲次按需回访 | F4 前半：小型 GPT、训练曲线、Tensor 形状与一个变量实验 |
| 3 · 使用现成模型，约 1 周 | [HF LLM Course](https://huggingface.co/learn/llm-course/en/chapter1/1) 的第 2 章作接口参考 | tokenizer、模型加载、上下文与生成；不另开一套完整课程 | 能用现成模型，并解释它与自己的小模型有什么区别 |
| 4 · 应用体验，1–2 周 | F2 或 F3 的单个官方教程 | 资料问答与图片问答二选一 | 固定测试样例，保留成功和失败结果 |

阶段 0 与[具身 Overview](../../EMBODIED-INTELLIGENCE/OVERVIEW/overview.md)交错阅读。F1 是共享准备，不必在两个栏目各做一次。

## 根据实践再选深入方向

| 主线 | 下一步 | 保留的资源 |
|---|---|---|
| 模型算法 | 完成 F4，再做一次任务适配 | [HF 微调文档](https://huggingface.co/docs/transformers/training)；需要系统训练研究时再选 CS336 |
| 应用与 Agent | 扩展 F2，增加一个可检查的只读工具 | [应用与智能体](../APPLICATIONS/applications.md)中的单个项目入口 |
| 推理系统 | F5：质量、时延、内存比较 | [推理与部署](../INFERENCE/inference.md)中的 llama.cpp |
| 多模态与机器人 | E3 后选择 E4 | [具身路线](../../EMBODIED-INTELLIGENCE/ROADMAP/roadmap.md)中的 LeRobot |

**文字学习替代路线：**[LLMs-from-scratch](https://github.com/rasbt/LLMs-from-scratch) 第 2–5 章可替代 Karpathy 实现部分；选定一种媒介后沿同一套代码推进。

MoE、完整偏好优化、大规模并行、世界模型和多智能体暂作知识地图中的选修主题。不要将每个专题页都理解为一门必修课。

## 选课与求职判断

先完成两个 Overview 和一个实现项目，再对照 3–5 个目标岗位。按“能解释、能复现、能独立改进”记录能力，使用项目结果决定下一步。

每个阶段先跑原例，再改一个因素，最后用未参与调试的样例评估。记录依赖、代码版本、结果和自己的取舍。

[F0–F5 项目步骤](../PROJECT-BASED-LEARNING/project-based-learning.md)
