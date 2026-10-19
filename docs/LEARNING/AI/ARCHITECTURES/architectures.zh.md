# 神经网络架构与计算结构

用输入结构、信息交互方式和计算代价理解架构。架构、训练目标与应用任务是三个维度：同一 Transformer 可以用于语言、视觉或生成动作。

## 知识地图与优先级

| 层次 | 主题 | 需要理解什么 |
|---|---|---|
| P0 | MLP、CNN 与残差 | 前馈、局部共享、感受野、残差与归一化 |
| P1 | 序列架构 | RNN/LSTM/GRU、Transformer、SSM/Mamba；状态压缩与上下文访问 |
| P1 | 图与几何结构 | GNN 消息传递、图 Transformer、对称性与等变性 |
| P2 | 高效注意力与混合结构 | 稀疏/线性注意力、混合序列层；FlashAttention 属于实现优化 |
| P2 | MoE 与条件计算 | 路由、专家、负载均衡；总参数与激活参数的差别 |

## 核心计算范式对比

| 架构族 | 归纳偏置（Inductive Bias） | 上下文交互方式 | 复杂度（序列长度 $L$） | 端侧硬件亲和性 |
|---|---|---|---|---|
| **CNN** | 局部性与平移不变性 | 局部卷积核滑动，感受野随层数线性增加 | 时间 $O(L)$，显存 $O(L)$ | 极高（权重复用、规则访存、NPU友好） |
| **Transformer** | 弱偏置（置换等变，依赖位置编码） | 全对全内容寻址（动态软路由） | 时间 $O(L^2)$，显存 $O(L^2)$ | 需硬件加速（推理 KV Cache 成为带宽与显存瓶颈） |
| **SSM / Mamba** | 连续时间递推与马尔可夫演化 | 选择性连续状态压缩与记忆更新 | 训练 $O(L)$，推理单步 $O(1)$ | 良好（自回归推理无无限增长的缓存换入开销） |
| **GNN** | 拓扑不变性与置换等变性 | 沿显式图邻接矩阵的消息传递与聚合 | 依赖图稀疏度与采样策略 | 中等（稀疏矩阵与不规则内存访存挑战） |

## 精选入口

<div class="grid cards" markdown>

-   :material-code-braces:{ .lg .middle } __实现主线 · Karpathy Zero to Hero__

    ---

    沿小型语言模型亲手实现自注意力、Multi-Head Attention 与前向传播计算图。

    [:octicons-arrow-right-24: 访问实现教程](https://karpathy.ai/zero-to-hero.html)

-   :material-eye:{ .lg .middle } __视觉结构查阅 · CS231n 2026__

    ---

    经典视觉网络架构参考，查阅卷积（CNN）、残差连接（ResNet）与注意力的对应讲次。

    [:octicons-arrow-right-24: 查阅官方讲义](https://cs231n.stanford.edu/schedule)

</div>

## 跟随实验 · 结构对照表，而非同时复现所有架构

1. 先选择已完成的 F1 或 F4 模型，标出输入、输出和参数。
2. 计算主要层的 Tensor 形状、参数量与一次前向过程。
3. 只替换一个模块或改变上下文长度，记录质量与资源变化。
4. 再阅读一条替代路线，列出适用条件，不仅比较参数量。

**验收**：一张结构图、一个受控修改、输入输出与复杂度说明；不要求训练 MoE 或安装 Mamba CUDA 扩展。

## 认知自测 · 关键思考

- [ ] 残差连接（$x + F(x)$）在数学上为反向传播梯度流（Gradient Highway）提供了怎样的恒等项保护？
- [ ] 为什么在许多边缘传感器实时监测或视觉检测任务中，轻量现代 CNN（如 ConvNeXt/MobileNet）依然比 ViT 更易部署且延迟更稳定？
- [ ] 当输入序列长度扩展至 32k 以上时，全注意力机制与状态空间模型各自的显存和计算瓶颈分别发生在哪里？

Transformer · Mamba / SSM · 图与几何学习；大模型训练与部署在基础模型。

## 相关专题

<div class="grid cards" markdown>

-   __Transformer__

    ---

    自注意力的张量骨架与数学细节

    [:octicons-arrow-right-24: 进入](../TRANSFORMER/transformer.md)

-   __Mamba / SSM__

    ---

    状态空间序列模型与选择性机制

    [:octicons-arrow-right-24: 进入](../MAMBA/mamba.md)

-   __图与几何学习__

    ---

    消息传递、等变性与图结构

    [:octicons-arrow-right-24: 进入](../GNN/gnn.md)

-   __基础模型__

    ---

    大模型的架构、训练、应用与部署

    [:octicons-arrow-right-24: 进入](../../FOUNDATION-MODELS/foundation-models.md)

-   __AI 全貌__

    ---

    AI 知识地图与模块关系

    [:octicons-arrow-right-24: 进入](../ai.md)

-   __学习路线__

    ---

    三轮学习路线与体验项目池

    [:octicons-arrow-right-24: 进入](../ROADMAP/roadmap.md)

</div>
