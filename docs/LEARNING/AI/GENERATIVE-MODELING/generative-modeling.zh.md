# 几何表示与生成建模

!!! abstract "这一页解决什么"
    模型怎样"生成"新数据，以及数据本身有怎样的几何结构。流形和潜空间描述数据与表示，扩散和流匹配是生成方法；它们也是机器人扩散策略、π0 类 VLA 动作头的数学基础。**专题选修**：不作为 LLM 或机器人入门的统一前置。

## 核心专题

<div class="grid cards" markdown>

-   :material-vector-curve:{ .lg .middle } __流形与潜空间__

    ---

    高维数据蕴含怎样的低维结构？环境空间、内在维数与潜变量。

    [:octicons-arrow-right-24: 进入专题](MANIFOLDS/manifolds.md)

-   :material-chart-scatter-plot:{ .lg .middle } __扩散模型__

    ---

    在多个噪声尺度上学习去噪与 score，再逐步采样生成。

    [:octicons-arrow-right-24: 进入专题](DIFFUSION/diffusion.md)

-   :material-arrow-decision:{ .lg .middle } __流匹配__

    ---

    学习一条从噪声到数据的速度场，用 ODE 数值积分生成。

    [:octicons-arrow-right-24: 进入专题](FLOW-MATCHING/flow-matching.md)

</div>

## 核心机制 · 生成方法全貌

所有生成模型都在回答同一个问题：怎样从一个简单分布（通常是高斯噪声）得到符合数据分布的样本。区别在于"怎么搬运"：

| 方法族 | 核心思路 | 优点 | 代价 |
|---|---|---|---|
| 自回归 | 把联合分布拆成逐个条件分布，按顺序生成 | 训练简单、似然精确 | 生成慢、只能顺序 |
| VAE | 编码到潜变量、再解码，优化证据下界 ELBO | 有显式潜空间 | 样本偏模糊 |
| GAN | 生成器骗过判别器 | 样本清晰、生成快 | 训练不稳定、模式坍塌 |
| Normalizing flow | 可逆变换，精确计算密度变化 | 精确似然 | 架构受可逆约束 |
| 扩散 / score | 学习逐步去噪（或对数密度梯度） | 质量高、训练稳定、能表示多峰 | 采样步数多 |
| 流匹配 | 直接回归从噪声到数据的速度场 | 训练更简单、路径更直、步数可更少 | 同样需要数值积分 |

Transformer 是网络架构，扩散和流匹配是生成方法——Transformer 可以作为其中的去噪网络（如 DiT）。流形假设是理解数据的视角，不是生成算法。

## 精选资源

<div class="grid cards" markdown>

-   :material-school:{ .lg .middle } __主课 · MIT 6.S184（2026）__

    ---

    同一套讲义、视频和 Labs 串起 ODE/SDE、score matching、流匹配与条件生成，避免在多套符号间切换。

    [:octicons-arrow-right-24: 打开课程](https://diffusion.csail.mit.edu/2026/index.html)

-   :material-github:{ .lg .middle } __Labs · MIT 6.S184 代码__

    ---

    Lab 1：ODE/SDE；Lab 2：flow / score matching（G2、G3）；Lab 3：图像与 VAE。

    [:octicons-arrow-right-24: 打开代码](https://github.com/eje24/iap-diffusion-labs/tree/2026)

-   :material-eye:{ .lg .middle } __补充 · CS231n 生成模型讲次__

    ---

    VAE / GAN 的基础概念，只作同一知识地图的补充。

    [:octicons-arrow-right-24: 打开课表](https://cs231n.stanford.edu/schedule)

</div>

**针对你的背景**：线性代数、信号、系统与几何经验有助于理解向量场和时间演化；概率、随机过程与 PyTorch 仍需通过例子核对。课程要求线性代数、多元微积分、基础概率和 PyTorch，不是零门槛入门课。

| 阶段 | 只做这些 | 完成标准 |
|---|---|---|
| 扫盲 | 本页的方法对比、课程第 1 讲图示 | 区分表示、架构与生成方法 |
| 几何体验 | G1：Swiss roll 降维 | 区分环境维数、内在维数与嵌入 |
| 数值基础 | Lab 1：ODE / SDE，先看轨迹再看推导 | 理解步长、确定性与随机性 |
| 学习与采样 | Lab 2：G2 与 G3 | 区分训练目标与采样过程，比较二维样本 |
| 方向深入 | 第 3-B / 4 讲与 Lab 3：条件、DiT、VAE | 在同一代码体系中完成一个图像实验 |

**概念阅读可以先做；编码放在 F1 或等价 PyTorch 能力之后。**

## 三个跟随实验

<div class="grid cards" markdown>

-   :material-cube-scan:{ .lg .middle } __G1 · 数据几何与流形展开__

    ---

    scikit-learn Swiss roll，CPU 即可对比线性 PCA 与非线性局部邻域方法。

    [:octicons-arrow-right-24: 打开官方示例](https://scikit-learn.org/stable/auto_examples/manifold/plot_swissroll.html)

-   :material-wave-square:{ .lg .middle } __G2 · 二维去噪与 score matching__

    ---

    MIT Lab 2 的 score matching 部分，观察不同噪声水平下的梯度场。

    [:octicons-arrow-right-24: 打开 Lab 代码](https://github.com/eje24/iap-diffusion-labs/tree/2026)

-   :material-ray-start-arrow:{ .lg .middle } __G3 · 向量场生成与 ODE 采样__

    ---

    MIT Lab 2 的 flow matching 部分，测试 Euler 步数对生成质量的影响。

    [:octicons-arrow-right-24: 打开 Lab 代码](https://github.com/eje24/iap-diffusion-labs/tree/2026)

</div>

先完成原例，每个实验只改一个因素。

## 自测

??? question "扩散模型和 GAN 各自的主要问题是什么？"
    GAN 训练不稳定、容易模式坍塌（只生成少数几种样本）；扩散模型训练稳定、覆盖全面，但采样需要几十到上千步，生成慢。

??? question "为什么机器人动作生成会用扩散或流匹配，而不是直接回归？"
    同一状态下可能有多种合理动作（多峰分布），直接回归会输出它们的平均；扩散 / 流匹配能从多峰分布中采样出其中一个合理动作。

??? question "潜空间扩散（如 Stable Diffusion）为什么要先编码到潜空间？"
    像素空间维度高、计算昂贵；先用 VAE 压缩到低维潜空间，在那里做扩散，再解码回图像，计算量大幅下降。

## 相关专题

<div class="grid cards" markdown>

-   __基础模型__

    ---

    模型的使用、训练与部署。

    [:octicons-arrow-right-24: 进入](../../FOUNDATION-MODELS/foundation-models.md)

-   __学习型策略（具身）__

    ---

    扩散策略与流匹配动作头的闭环使用。

    [:octicons-arrow-right-24: 进入](../../EMBODIED-INTELLIGENCE/POLICY-LEARNING/policy-learning.md)

</div>
