# Mamba 与状态空间序列模型

!!! abstract "这一页解决什么"
    用"状态如何更新、历史如何压缩、计算如何随长度变化"来理解另一类序列模型。你熟悉的状态空间方程和系统辨识是很好的起点，你也已经在盾构项目里用 Mamba 替换过 Transformer。**P1 认识，P2 深入**。

## 核心机制

**线性状态空间模型**：连续形式 h'(t) = A h(t) + B x(t)，y(t) = C h(t)。离散化（如零阶保持）后变成递推 h_t = Ā h_{t−1} + B̄ x_t，和控制理论里的状态方程一致。

**两种计算方式**：A、B、C 不随输入变化（LTI，如 S4）时，可以把递推展开成一个长卷积并行训练；推理时用递推，每步 O(1)、不需要随长度增长的缓存。

**选择性（Mamba 的关键）**：让 B、C 和步长 Δ 依赖当前输入，模型可以决定"记住什么、忽略什么"，弥补 LTI 模型不能按内容选择信息的弱点。代价是不能再写成卷积，改用硬件感知的并行扫描实现。

**和注意力的取舍**：注意力保留全部历史、可以精确回看任意位置，但 KV cache 随长度增长；SSM 把历史压缩成固定大小的状态，长序列省内存、推理快，但精确回忆远处细节较弱。所以混合架构（少量注意力层 + 大量 SSM 层）越来越常见。

**Mamba-2 与 SSD**：揭示了结构化 SSM 与某类线性注意力的对偶关系，训练效率更高。

**和物理状态的区别**：SSM 的"状态"是学到的抽象向量，不一定对应真实的物理状态；不要把它当作可解释的系统辨识结果。

| 层次 | 主题 | 需要区分 |
|---|---|---|
| P1 | 线性状态空间、离散化与递推 | 数学状态与真实物理状态不一定相同 |
| P1 | 选择性 SSM / Mamba | 根据输入选择信息，不等于固定线性滤波器 |
| P2 | Mamba-2、SSD 与混合序列模型 | 结构联系、训练与推理代价、精度 |
| P2 | 时序与长序列应用 | 数据、采样率、长度与任务决定比较条件 |

## 精选资源

<div class="grid cards" markdown>

-   :material-github:{ .lg .middle } __官方主入口 · Mamba 仓库__

    ---

    原作者论文、架构图与硬件感知并行扫描的 CUDA 实现。

    [:octicons-arrow-right-24: 打开仓库](https://github.com/state-spaces/mamba)

-   :material-video:{ .lg .middle } __直觉解析 · Mamba 原理解读__

    ---

    中文视频，直观对比连续状态空间与注意力的结构差异。

    [:octicons-arrow-right-24: 观看视频](https://www.bilibili.com/video/BV1QrBbBhEyN/)

</div>

## 小实验：先建立可比较的基线

在时序与科学 AI 专题的数据上（见下方相关专题），先比较持续值 / 线性预测与一个小网络；再加入 Mamba，固定时间划分、输入窗口、训练预算和设备，记录预测质量、时延与内存。

**验收**：解释状态压缩与注意力访问历史的区别；不根据渐近复杂度直接宣称实际更快或已经替代 Transformer。

## 自测

??? question "为什么 LTI 的 SSM 可以写成卷积，而 Mamba 不行？"
    LTI 模型的 Ā、B̄、C 固定，输出是输入与一个固定核的卷积；Mamba 的参数随输入变化，每一步的"核"都不同，只能用扫描计算。

??? question "长序列推理时，Mamba 相对 Transformer 省在哪里？"
    Transformer 的 KV cache 随长度线性增长，每生成一个 token 都要读全部缓存；Mamba 只维护固定大小的状态，每步计算和内存都是常数。

??? question "什么任务上 Mamba 可能不如注意力？"
    需要精确复制或检索远处具体内容的任务（如长文中查找某个电话号码），固定大小的状态难以无损保存所有细节。

## 相关专题

<div class="grid cards" markdown>

-   __网络架构__

    ---

    CNN、Transformer、SSM、GNN 的计算范式对比。

    [:octicons-arrow-right-24: 进入](../ARCHITECTURES/architectures.md)

-   __时序与科学 AI__

    ---

    在你熟悉的信号数据上做对比实验。

    [:octicons-arrow-right-24: 进入](../SCIENTIFIC-AI/scientific-ai.md)

</div>
