# 基础模型

先理解 LLM 与多模态模型在解决什么问题，再通过小实验认识模型、数据和应用，最后选择方向深入。

<div class="grid cards" markdown>

-   :material-rocket-launch:{ .lg .middle } __Get Started · 两周速成__

    ---

    14 天闭环：小 GPT、LoRA 微调、带来源的检索问答，量化后离线跑在 Jetson 上。

    [开始速成](GET-STARTED/get-started.md)

-   :material-compass-outline:{ .lg .middle } __Overview · 全貌导读__

    ---

    从零认识术语、系统组成与方向关系，配合中文导论和可视化视频。

    [从这里开始](OVERVIEW/overview.md)

-   :material-map-marker-path:{ .lg .middle } __Roadmap · 学习路线__

    ---

    按“扫盲 → 基础 → 使用 → 适配 → 深入”推进，明确必学、选学和完成标准。

    [查看路线](ROADMAP/roadmap.md)

-   :material-hammer-wrench:{ .lg .middle } __Guided Projects · 跟随式项目__

    ---

    从 PyTorch 小模型到个人资料问答，每个项目都有教程、步骤、个人修改和验收要求。

    [开始实践](PROJECT-BASED-LEARNING/project-based-learning.md)

</div>

## 知识地图

| 单元 | 核心问题 | 学习层次 |
|---|---|---|
| [架构与表示](ARCHITECTURE/architecture.md) | token、嵌入和 Attention 如何组成模型？ | 基础；架构变体选学 |
| [训练与适配](TRAINING/training.md) | 预训练、SFT、LoRA 和偏好优化分别改变什么？ | 基础；大规模训练选学 |
| [应用与智能体](APPLICATIONS/applications.md) | 如何连接资料、工具与任务，并评估结果？ | 基础与实践 |
| [多模态与行动接口](MULTIMODAL/multimodal.md) | 图像、语言、状态和动作怎样连接？ | 导论与体验；VLA 适配进阶 |
| [推理与部署](INFERENCE/inference.md) | 怎样权衡质量、延迟、吞吐和内存？ | 基础；运行时优化进阶 |
| [世界模型](WORLD-MODELS/world-models.md) | 怎样预测环境变化，并用于规划或学习？ | 进阶选修 |

!!! tip "如何使用这个栏目"
    主课只保留李宏毅导论与 Karpathy 实现两条连续路线。先做 F0；F1 用于按需补 PyTorch，F4 的小型 GPT 是实现主项目；F2/F3 二选一体验，F5 以后按方向选择。文档在遇到具体接口时查阅，不需要逐个学完。

## 与具身智能的连接

- **技能调用**：语言模型选择目标或技能，机器人控制器负责执行。
- **VLA**：视觉、语言和机器人状态共同条件化动作策略，需要机器人数据与闭环评估。
- **世界模型**：预测动作带来的环境变化，可支持规划或策略学习。

在[具身智能全貌导读](../EMBODIED-INTELLIGENCE/OVERVIEW/overview.md)中查看这些模块如何进入机器人系统。
