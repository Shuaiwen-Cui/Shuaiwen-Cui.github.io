# Learning Fundamentals, Optimization & Classical ML

This is shared vocabulary. Use small experiments to diagnose probability, optimization and evaluation gaps while reusing your mathematics/engineering background. P0 is shared foundation, P1 orientation, P2 specialization.

## Map and priority

| Depth | Topic | Understand |
|---|---|---|
| P0 | Mathematics and statistics | Vectors, matrices, derivatives, conditional probability, likelihood and moments; association vs. causation |
| P0 | Learning and data | Supervision, loss, train/validation/test, generalization, leakage and imbalance |
| P0 | Optimization and networks | Backprop, SGD/Adam, schedules, normalization, regularization and initialization |
| P1 | Classical methods | Linear/logistic models, SVM, trees/ensembles/boosting, clustering and PCA; baseline first |
| P2 | Probabilistic and causal modeling | Bayesian inference, graphical models, interventions and identifiability as needed |

## Core mechanistic intuition · From classical ML to neural nets
- **Feature space & inductive bias**: Linear models assume hyperplane separation; decision trees perform axis-aligned orthogonal partitions (invariant to monotonic feature scaling); kernel SVMs map to implicit reproducing Hilbert spaces; neural networks learn continuous coordinate warpings via layered function compositions.
- **Optimization & backprop**: Backpropagation is dynamic programming caching intermediate jacobians of the chain rule. Gradient descent stability depends on the curvature (Hessian conditioning) of the loss surface—the primary reason learning rate warmup, normalization layers, and adaptive moment estimation (Adam) are required.

## Selected resources

<div class="grid cards" markdown>

-   :material-book-open-page-variant:{ .lg .middle } __Backbone · Karpathy Zero to Hero__

    ---

    Code-first backbone to fill gaps in loss, gradient backprop, training loops and validation.

    [:octicons-arrow-right-24: Access course material](https://karpathy.ai/zero-to-hero.html)

-   :material-chart-scatter-plot:{ .lg .middle } __Official experiment · scikit-learn classifier comparison__

    ---

    Inspect decision boundaries and variance across linear, kernel and tree models on shared 2D datasets.

    [:octicons-arrow-right-24: View official example](https://scikit-learn.org/stable/auto_examples/classification/plot_classifier_comparison.html)

</div>

## Guided experiment · AI-1 · Baselines and overfitting

1. Run the classifier comparison; verify that preprocessing fits training data only.
2. Keep one linear and one tree model; vary only noise or model complexity.
3. Tune on validation and report test once; use grouped/temporal splits when required.
4. Show learning curves, errors and at least two appropriate metrics.

**Acceptance**: CPU is sufficient. Record three settings and explain bias, variance, leakage and metric choices.

## Check yourself · Critical thinking

- [ ] Why are tree-based models largely invariant to feature scaling, whereas gradient-based neural networks and distance metrics (KNN/SVM) require strict normalization?
- [ ] In highly imbalanced tasks (e.g., 99:1), why is accuracy misleading, and how do precision-recall AUC and cost-sensitive losses resolve this?
- [ ] When splitting data with temporal correlation or sensor-group clustering, why does naive random shuffling cause optimistic performance leakage?

See trees and ensembles; use evaluation throughout.

## Related topics

<div class="grid cards" markdown>

-   __trees and ensembles__

    ---

    Decision trees and ensembles in detail

    [:octicons-arrow-right-24: Open](../DECISION-TREE/decision-tree.md)

-   __evaluation__

    ---

    Data quality, calibration, uncertainty, reliability

    [:octicons-arrow-right-24: Open](../EVALUATION/evaluation.md)

-   __AI overview__

    ---

    The AI knowledge map and how modules relate

    [:octicons-arrow-right-24: Open](../ai.md)

-   __Roadmap__

    ---

    Three-round route and project pool

    [:octicons-arrow-right-24: Open](../ROADMAP/roadmap.md)

</div>
