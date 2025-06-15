# 决策树

## 决策树的定义

决策树是一种用于分类和回归的监督学习算法。它通过将数据集分割成更小的子集来构建一个树状模型，从而实现对数据的预测。每个内部节点代表一个特征的测试，每个叶节点代表一个类别或数值。
决策树的构建过程通常包括以下几个步骤：

1. **选择最佳特征**：使用某种准则（如信息增益、基尼指数等）来选择最能区分数据的特征。
2. **分割数据集**：根据选定的特征将数据集分割成子集。
3. **递归构建子树**：对每个子集重复上述步骤，直到满足停止条件（如达到最大深度或叶节点样本数小于某个阈值）。
4. **生成决策树**：将所有子树组合成一棵完整的决策树。
5. **剪枝**：通过去除一些不必要的分支来简化决策树，防止过拟合。
6. **预测**：使用生成的决策树对新数据进行预测。

## 决策树的优点和缺点

### 优点
- **易于理解和解释**：决策树的结构直观，易于可视化和解释。
- **处理非线性关系**：能够处理特征之间的非线性关系。
- **无需特征缩放**：不需要对特征进行标准化或归一化处理。
- **适用于大数据集**：可以处理大规模数据集，且训练速度较快。

### 缺点
- **过拟合**：决策树容易对训练数据过拟合，尤其是在树的深度较大时。
- **不稳定性**：小的变化可能导致决策树结构的显著变化。
- **偏向于多值特征**：决策树可能偏向于具有更多类别的特征。
- **难以捕捉复杂模式**：对于某些复杂的模式，决策树可能无法很好地拟合。

## 学习链接

<div class="grid cards" markdown>

-   :fontawesome-brands-bilibili:{ .lg .middle } __决策树🏆__

    ---


    [:octicons-arrow-right-24: <a href="https://www.bilibili.com/video/BV1Xp4y1U7vW/?spm_id_from=333.788.recommend_more_video.0&vd_source=5a427660f0337fedc22d4803661d493f" target="_blank"> Portal </a>](#)

-   :material-file:{ .lg .middle } __这是我见过关于决策树最详细的分析🏆__

    ---


    [:octicons-arrow-right-24: <a href="https://zhuanlan.zhihu.com/p/21018652275" target="_blank"> Portal </a>](#)

-   :material-file:{ .lg .middle } __【机器学习】决策树（上）——ID3、C4.5、CART（非常详细）🏆__

    ---


    [:octicons-arrow-right-24: <a href="https://zhuanlan.zhihu.com/p/85731206" target="_blank"> Portal </a>](#)

-   :material-file:{ .lg .middle } __【机器学习】决策树（中）——Random Forest、Adaboost、GBDT （非常详细）🏆__

    ---


    [:octicons-arrow-right-24: <a href="https://zhuanlan.zhihu.com/p/86263786" target="_blank"> Portal </a>](#)

-   :material-file:{ .lg .middle } __【机器学习】决策树（下）——XGBoost、LightGBM（非常详细）🏆__

    ---


    [:octicons-arrow-right-24: <a href="https://zhuanlan.zhihu.com/p/87885678" target="_blank"> Portal </a>](#)

</div>

## 发展脉络

基本树：基本树（包括 ID3、C4.5、CART）
进阶树：进阶树（包括 Random Forest、Adaboost、GBDT）
高级树：高级树（包括 Xgboost、LightGBM）
