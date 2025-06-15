# DECISION TREE

## Definition of Decision Trees

A decision tree is a supervised learning algorithm used for classification and regression. It builds a tree-like model by splitting the dataset into smaller subsets, enabling predictions to be made based on the data. Each internal node represents a test on a feature, and each leaf node represents a class label or a numerical value.

The process of building a decision tree typically involves the following steps:

1. **Feature Selection**: Use a criterion (e.g., information gain, Gini index) to select the feature that best separates the data.
2. **Data Splitting**: Split the dataset into subsets based on the selected feature.
3. **Recursive Subtree Construction**: Repeat the above steps for each subset until a stopping condition is met (e.g., maximum depth or minimum number of samples per leaf).
4. **Tree Generation**: Combine all subtrees into a complete decision tree.
5. **Pruning**: Simplify the tree by removing some unnecessary branches to prevent overfitting.
6. **Prediction**: Use the resulting decision tree to make predictions on new data.

## Advantages and Disadvantages of Decision Trees

### Advantages
- **Easy to Understand and Interpret**: The structure of a decision tree is intuitive and easy to visualize and explain.
- **Handles Nonlinear Relationships**: Capable of modeling nonlinear interactions between features.
- **No Need for Feature Scaling**: Standardization or normalization of features is not required.
- **Suitable for Large Datasets**: Can handle large-scale datasets and has relatively fast training speed.

### Disadvantages
- **Overfitting**: Decision trees are prone to overfitting, especially when the tree is deep.
- **Instability**: Small changes in the data can lead to significant changes in the tree structure.
- **Bias Toward Multivalued Features**: Trees may favor features with more categories.
- **Limited in Capturing Complex Patterns**: May struggle to model more intricate patterns in data.

## Study Links

<div class="grid cards" markdown>

-   :fontawesome-brands-bilibili:{ .lg .middle } __Decision Tree🏆__

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

## Development History

Basic Tree Structures: ID3 C4.5 CART
Advanced Tree Structures: Random Forest Adaboost GBDT
More Advanced Tree Structures: XGBoost LightGBM