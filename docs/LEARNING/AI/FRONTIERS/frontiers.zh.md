# 前沿技术索引 · 放回知识体系

!!! note "范围与依据"
    截至 2026-10-18 核对的代表技术族，不是热度排行榜或所有研究方向的穷举。用课程、原作者代码和文档确认范围；模型名作实例，知识位置、训练目标和限制优先。

## 架构、表示与生成

| 技术族 | 一句话定位 | 在哪里学 | 第一轮边界 |
|---|---|---|---|
| Transformer、ViT、DiT | 相似组件用于语言、视觉或扩散任务，架构不决定唯一训练目标 | [架构](../ARCHITECTURES/architectures.md) / [生成](../GENERATIVE-MODELING/generative-modeling.md) | 基本注意力与数据流 |
| 高效注意力、长上下文 | 结构、近似、位置与实现共同影响成本和有效上下文 | [Transformer](../TRANSFORMER/transformer.md) / [系统](../SYSTEMS/systems.md) | 区分线性注意力与 FlashAttention |
| MoE 与条件计算 | 路由到部分专家，需平衡质量、激活计算与通信 | [架构](../ARCHITECTURES/architectures.md) | 理解总参数与激活参数 |
| SSM、Mamba、混合序列层 | 用状态更新与选择性机制处理序列 | [Mamba](../MAMBA/mamba.md) | 先与 RNN/Attention 做概念对照 |
| 自监督、对比、CLIP、DINO、MAE | 从多视图、教师或掩码构造不同监督目标 | [表示学习](../REPRESENTATION/representation.md) | 先冻结特征与微调实验 |
| JEPA 与预测式表示 | 在表征空间预测，和像素生成目标不同 | [表示](../REPRESENTATION/representation.md) / [世界模型](../../FOUNDATION-MODELS/WORLD-MODELS/world-models.md) | 阅读 [V-JEPA 2 原作者代码](https://github.com/facebookresearch/vjepa2)，不要求重训练 |
| VAE、GAN、normalizing flow | 潜变量、对抗与可逆变换等生成路线 | [生成地图](../GENERATIVE-MODELING/generative-modeling.md) | 分清目标与采样方式 |
| 扩散、flow matching、rectified flow | 在噪声/概率路径上学习生成过程，具体目标与路径不同 | [扩散](../GENERATIVE-MODELING/DIFFUSION/diffusion.md) / [流匹配](../GENERATIVE-MODELING/FLOW-MATCHING/flow-matching.md) | 同一套低维 Labs |
| 一致性/蒸馏与少步采样 | 以额外训练或目标设计压缩生成步骤 | [生成](../GENERATIVE-MODELING/generative-modeling.md) / [系统](../SYSTEMS/systems.md) | 暂了解用途，不宣称无代价加速 |
| 离散扩散与非自回归生成 | 在离散状态上定义不同于逐 token 自回归的生成过程 | [生成](../GENERATIVE-MODELING/generative-modeling.md) | MIT 6.S184 第 5 讲选修 |

## 感知、推理、决策与系统

| 技术族 | 一句话定位 | 在哪里学 | 第一轮边界 |
|---|---|---|---|
| 提示式分割、视觉语言与视频模型 | 从特定标签任务走向条件化与多模态处理 | [视觉](../CV/cv.md) / [多模态](../../FOUNDATION-MODELS/MULTIMODAL/multimodal.md) | 明确任务/指标与空间局限 |
| NeRF、3DGS、隐式几何与点云 | 不同的三维形状或外观表示 | [视觉与三维](../CV/cv.md) | 新视角质量与几何正确分开验证 |
| 图 Transformer、等变网络 | 显式利用关系和对称性 | [图与几何](../GNN/gnn.md) | 小图 GCN 与无边基线 |
| SFT、LoRA/PEFT、偏好优化 | 任务适配与行为优化的方法，不能只比较缩写 | [训练与适配](../../FOUNDATION-MODELS/TRAINING/training.md) | 一个适配任务与独立评估 |
| RLVR、推理预算与测试时计算 | 用可验证反馈或额外推理计算改善特定任务 | [训练](../../FOUNDATION-MODELS/TRAINING/training.md) / [决策](../DECISION-LEARNING/decision-learning.md) | 认识奖励、验证与成本；GRPO 等算法细节后深入 |
| RAG、工具、记忆、Agent、多智能体 | 组织模型与外部资料/接口形成任务系统 | [应用与智能体](../../FOUNDATION-MODELS/APPLICATIONS/applications.md) | 先一个可检查的工具工作流 |
| 模仿、离线/模型式 RL、VLA | 从演示或交互学习动作并接入物理系统 | [决策](../DECISION-LEARNING/decision-learning.md) / [机器人策略](../../EMBODIED-INTELLIGENCE/POLICY-LEARNING/policy-learning.md) | 先 ACT 与闭环评估 |
| 世界模型、潜空间规划 | 预测环境变化，并检验对决策的帮助 | [世界模型](../../FOUNDATION-MODELS/WORLD-MODELS/world-models.md) | 预测能力与任务成功分开 |
| 时序基础模型、PINN、神经算子 | 时间数据迁移、物理约束、函数间映射 | [科学 AI](../SCIENTIFIC-AI/scientific-ai.md) | 一个基线与一个受控实验 |
| 量化、蒸馏、稀疏、并行与运行时 | 从算法效率走向硬件实际效率 | [系统](../SYSTEMS/systems.md) | 一台设备上测质量、内存与时延 |
| 持续/主动/联邦学习 | 分别处理数据/任务变化、标注预算或分布式数据约束 | [表示与范式](../REPRESENTATION/representation.md) | 先理解问题，不同时重建训练栈 |
| 校准、OOD、解释、数据污染与可靠性 | 检查信心、分布变化、证据与失效 | [评估与可信 AI](../EVALUATION/evaluation.md) | 所有项目共用检查 |

## 校验覆盖的少量来源

- [CS231n 2026](https://cs231n.stanford.edu/schedule)：将视觉、表示、生成、三维与视觉语言放进同一课程结构；用于检查知识位置。
- [CS336 2026](https://cs336.stanford.edu/)：数据、训练系统、对齐与推理 RL；要求较强 Python、深度学习与系统基础，作为后续方向课。
- [MIT 6.S184 2026](https://diffusion.csail.mit.edu/2026/index.html)：连续与离散生成、采样、潜空间及可跟随 Labs。

这些来源是分类与范围依据，不构成新增必修清单。新模型出现时，先判断其任务、架构、训练目标、数据、系统代价与证据，再归入已有模块；只有形成独立学习问题时才新增页面。

[AI 全貌](../ai.md) · [学习路线](../ROADMAP/roadmap.md)
