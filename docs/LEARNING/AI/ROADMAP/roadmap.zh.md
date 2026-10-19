# AI · 学习路线

!!! tip "时间紧时"
    先做 [两周速成](../GET-STARTED/get-started.md)，再按本页三轮补全。速成覆盖第二轮的共享基础和一个部署作品。

知识地图负责覆盖，路线负责取舍。结合你已有的数学、信号处理、几何与嵌入式背景，先查缺口，再做一次连续实现，随后选择一到两个跨方向体验。

## 分三轮推进

| 轮次 | 内容 | 主资源与范围 | 完成标准 |
|---|---|---|---|
| 第一轮 · 全貌 | AI 首页 + 各模块概念表 + 前沿索引 | 只读输入输出、目标、关系与局限；约 1–2 周的选读 | 看见一个术语能判断它是架构、学习目标、任务还是系统方法 |
| 第二轮 · 共享基础与一次实现 | P0 数据/优化/评估，按需补 PyTorch | 已有 Karpathy/F4 主线；基础够用时跳过重学；约 2–4 周参考 | 独立跑训练/验证，改一个因素，解释模型与结果 |
| 第三轮 · 选一到两个体验 | 视觉表示、生成、图、RL、科学 AI 中选择 | 每个方向只跟一个官方例子；不并行完整上多门课 | 用结果和兴趣选择求职深入方向 |

时间是规划参考，按实际缺口和可用时间调整。第一轮不要求先学完概率图、拓扑、SDE、完整 SLAM 或大规模预训练。

## 体验项目池：选择，不是逐项打卡

| 项目 | 入口 | 前置与设备 | 应得到的经验 |
|---|---|---|---|
| AI-1 · 简单基线 | [基础](../FOUNDATIONS/foundations.md) | Python、数据划分；CPU | 过拟合、指标与数据泄漏 |
| AI-2 · 表示迁移 | [表示](../REPRESENTATION/representation.md) | F1 等价训练能力；GPU 更方便 | 冻结/微调、标签量与域偏移 |
| G1–G3 · 几何与生成 | [生成专题](../GENERATIVE-MODELING/generative-modeling.md) | G1 可 CPU；G2/G3 需概率与 PyTorch | 数据几何、去噪、向量场与采样 |
| AI-3 · 小图 GCN | [图学习](../GNN/gnn.md) | 图概念、PyTorch；小图 CPU | 关系结构、图划分与基线 |
| AI-4 · CartPole | [决策学习](../DECISION-LEARNING/decision-learning.md) | MDP、训练循环；CPU 可起步 | 交互、探索与重复评估 |
| AI-5 · 信号/物理实验 | [科学 AI](../SCIENTIFIC-AI/scientific-ai.md) | 时间划分；算子/PINN 后续 | 衔接已有研究，验证物理与泛化 |
| AI-6 / AI-7 · 系统与可靠性 | [工程](../SYSTEMS/systems.md) / [评估](../EVALUATION/evaluation.md) | 复用以上任意模型 | 测性能、可靠性与实际局限 |

**对你的默认建议**：先完成共享基础自查与 F4 小实现；然后在 AI-2 和 G2/G3 中选一个，在 AI-5 或具身 E3 中选一个。图或 RL 若更感兴趣可替换。AI-6/AI-7 是已有项目上的检查，不再新增一个大项目。

## 核心资源与学习预算

<div class="grid cards" markdown>

-   :material-school:{ .lg .middle } __主线基石 · 概念与纯手写实现__

    ---

    李宏毅导论打通全局视野；Karpathy Zero to Hero 亲手写出小型 GPT。已会内容自测跳过。

    [:octicons-arrow-right-24: 访问 Karpathy 系列](https://karpathy.ai/zero-to-hero.html)

-   :material-eye:{ .lg .middle } __视觉与表示 · CS231n 2026__

    ---

    只选读对应卷积、残差连接与自监督讲次，使用公开讲义补齐表征视角，不强制通刷作业。

    [:octicons-arrow-right-24: 查阅 Stanford 课表](https://cs231n.stanford.edu/schedule)

-   :material-chart-bell-curve-cumulative:{ .lg .middle } __几何与生成 · MIT 6.S184__

    ---

    一套讲义与 Labs 串通扩散、流匹配与潜空间 ODE 演化，避免在多套符号与代码间切换。

    [:octicons-arrow-right-24: 访问主课与 Labs](https://diffusion.csail.mit.edu/2026/index.html)

-   :material-compass-rose:{ .lg .middle } __专精方向 · 图 / 强化学习 / 科学 AI__

    ---

    选定深入方向后再进入对应主课（CS224W / David Silver / NeuralOperator），其他保留知识定位。

    [查看体验项目池](#体验项目池选择不是逐项打卡)

</div>

## 何时结合岗位深入

全貌和至少一个实现完成后，对照 3–5 个真实岗位，比较兴趣、可迁移能力、算力/设备与可展示成果。选择一条主线做持续改进，而不是把“认识术语”写成“掌握技术”。

LLM 后训练、RAG、Agent、推理系统沿[基础模型路线](../../FOUNDATION-MODELS/ROADMAP/roadmap.md)；VLA、控制与迁移沿[具身路线](../../EMBODIED-INTELLIGENCE/ROADMAP/roadmap.md)。

[返回 AI 全貌](../ai.md) · [前沿技术索引](../FRONTIERS/frontiers.md)
