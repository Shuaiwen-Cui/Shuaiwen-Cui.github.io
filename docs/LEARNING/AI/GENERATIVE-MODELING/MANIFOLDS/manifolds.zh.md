# 流形与潜空间

!!! abstract "这一页解决什么"
    高维数据为什么往往"其实没那么高维"，以及模型学到的潜空间意味着什么。你的几何与线性代数背景可以直接复用，你在 Manifold Dynamics 中学到的潜在状态也属于这个视角。**专题选修**；第一轮不要求学微分几何或拓扑。

## 核心机制

| 概念 | 含义 | 注意 |
|---|---|---|
| 环境空间 | 数据数值所在的空间 | 一张 224×224 彩色图片是 15 万维的向量 |
| 流形 | 局部可以用低维欧氏坐标描述的空间 | 复杂数据不一定严格构成单一光滑流形 |
| 流形假设 | 真实数据集中在高维空间中的低维结构附近 | 是有用的建模假设，需要结合数据检查 |
| 潜空间 | 模型使用或学到的内部变量空间 | 不保证低维、可解释，或忠实保留全部几何 |

**Swiss roll 的启示**：点位于三维空间，但沿卷曲曲面的变化只需两个参数描述。三维直线距离会把不同卷层的点误判为邻居，沿曲面的距离（测地距离）才反映真实相似度。

**线性 vs. 非线性降维**：PCA 找方差最大的线性方向，展不开卷曲的流形；Isomap、LLE、UMAP 等利用局部邻域关系做非线性展开。t-SNE 和 UMAP 的可视化图中，簇之间的距离和簇的大小通常没有意义。

**神经网络与表示**：深度网络可以看作一连串坐标变换，把纠缠在一起的类别逐层"拉开"，最后变成线性可分。好的表示让下游任务变简单。

## 精选资源

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __直觉主读 · Chris Olah__

    ---

    Neural Networks, Manifolds, and Topology：用图示看网络怎样变换表示。把结论当作特定网络和数据的例子，不当作普遍定理。

    [:octicons-arrow-right-24: 打开文章](https://colah.github.io/posts/2014-03-NN-Manifolds-Topology/)

-   :material-chart-scatter-plot:{ .lg .middle } __跟随 G1 · scikit-learn Swiss roll__

    ---

    官方可下载代码与 Notebook，CPU 即可。

    [:octicons-arrow-right-24: 打开示例](https://scikit-learn.org/stable/auto_examples/manifold/plot_swissroll.html)

</div>

## G1：把卷曲的数据展开

1. 跑官方 Swiss roll 原例，查看三维点云与颜色对应关系。
2. 明确环境维数与生成数据的参数维数。
3. 在原例基础上加入 PCA 作为线性基线，与所选非线性方法比较。
4. 只改变邻居数或噪声，观察展开是否破坏邻域、产生短路或分裂。
5. 保存原始点云、二维嵌入和参数表；在留出的邻域检查中比较局部结构，不只看图形是否好看。

**验收**：至少三组设置，解释局部邻域与全局形状的差别；说明二维可视化不能证明发现了"真实流形"。

## 自测

??? question "为什么 PCA 展不开 Swiss roll？"
    PCA 只能做线性投影，找到的是方差最大的平面方向；卷曲结构需要非线性映射才能展平，线性投影会把不同卷层压在一起。

??? question "邻居数 k 设得太大会怎样？"
    邻域会跨越相邻卷层形成"短路"，把实际上很远的点连在一起，展开结果折叠或错乱。

??? question "为什么不能根据 t-SNE 图里两个簇的距离判断它们有多不同？"
    t-SNE 只保留局部邻域关系，簇间距离和簇大小受参数影响很大，不对应原空间的距离。

## 相关专题

<div class="grid cards" markdown>

-   __扩散模型__

    ---

    在潜空间里做生成。

    [:octicons-arrow-right-24: 进入](../DIFFUSION/diffusion.md)

-   __表示学习__

    ---

    表示从哪里来、怎样迁移。

    [:octicons-arrow-right-24: 进入](../../REPRESENTATION/representation.md)

</div>
