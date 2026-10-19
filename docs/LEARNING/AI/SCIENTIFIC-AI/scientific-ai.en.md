# Time Series & Scientific Machine Learning

This directly connects signals, identification and structural research. Start temporal splits and classical baselines, then distinguish data prediction, physical constraints and operator learning.

## Map and priority

| Depth | Topic | Understand |
|---|---|---|
| P1 | Time-series tasks and baselines | Forecasting, classification, anomalies and imputation; persistence/linear models, frequency features and leakage |
| P1 | Temporal networks | TCN, RNN, Transformers and SSM; locality, long range and multivariate inputs |
| P2 | Time-series foundation models | Pretraining data, sample rates, normalization, transfer and uncertainty |
| P2 | PINNs and inverse problems | Equation/boundary residuals as losses; accuracy and superiority require validation |
| P2 | Neural operators and hybrids | FNO/DeepONet learn function maps; test resolution, parameters and conservation |

## Core mechanistic intuition · From temporal representations to physical operators

- **Evolution of deep time series**:
    1. *Recurrent models (RNN/LSTM)*: Hidden state recurrence, prone to vanishing gradients, struggles with long-range vibrational phase preservation.
    2. *Point-wise self-attention (Vanilla Time Transformer)*: Treating each timestamp as a token disrupts temporal continuity, risks overfitting noise, and incurs quadratic cost on long sensor records.
    3. *Patching & Channel Independence (PatchTST/iTransformer)*: Aggregating consecutive time steps into sub-series patches preserves local context. Treating channels independently or projecting variables as tokens significantly improves generalization across sensor arrays.
- **PINNs vs. Neural Operators (FNO/DeepONet)**:
    - **PINN (Loss-constrained)**: Parameterizes a specific trajectory $u(t, x)$ under given parameters, penalized by PDE/ODE residuals. Requires no vast dataset but needs re-optimization per condition, often vulnerable to loss landscape imbalances.
    - **Neural Operators (Operator mapping)**: Maps between infinite-dimensional function spaces (e.g., input excitation force field $\to$ structural response field). Fourier layers compute global convolutions in spectral space, achieving **mesh independence (zero-shot super-resolution)**.

## Selected resources

<div class="grid cards" markdown>

-   :material-function-variant:{ .lg .middle } __Operator entry · NeuralOperator__

    ---

    Begin with small Darcy flow data to understand spectral convolutions and mesh-independent mappings.

    [:octicons-arrow-right-24: Official Quickstart](https://neuraloperator.github.io/dev/index.html)

-   :material-atom:{ .lg .middle } __Alternative PINN elective · DeepXDE__

    ---

    Physics-informed library reference; run a 1D Poisson example to examine PDE residual penalties.

    [:octicons-arrow-right-24: Official Forward Demos](https://deepxde.readthedocs.io/en/stable/demos/pinn_forward.html)

</div>

## Guided experiment · AI-5 · Temporal baseline, then one physics model

1. Define a signal task; split by time/device so neighboring windows do not leak.
2. Compare persistence, linear/feature baselines and one small network under changed noise/conditions.
3. For scientific AI, follow one Darcy or Poisson example; verify units/boundaries.
4. Test held-out parameters/resolutions against numerical/analytic references, not training loss alone.

**Acceptance**: Keep splits, conditions, errors and cost; evaluate physics and prediction separately without assuming deep models win.

## Check yourself · Critical thinking

- [ ] In multi-sensor monitoring, why does a channel-independent design often outperform complicated cross-channel attention architectures?
- [ ] In PINN optimization, what gradient pathologies occur when PDE residual losses and boundary condition losses differ by multiple orders of magnitude?
- [ ] How does the discrete Fourier transform and weight truncation in FNO relate mathematically to evaluating a continuous Green's function integral?

SSM and graphs are later model choices; retain domain-specific detail in its existing section.

## Related topics

<div class="grid cards" markdown>

-   __SSM__

    ---

    State-space sequence models and selectivity

    [:octicons-arrow-right-24: Open](../MAMBA/mamba.md)

-   __graphs__

    ---

    Message passing, equivariance and graph structure

    [:octicons-arrow-right-24: Open](../GNN/gnn.md)

-   __AI overview__

    ---

    The AI knowledge map and how modules relate

    [:octicons-arrow-right-24: Open](../ai.md)

-   __Roadmap__

    ---

    Three-round route and project pool

    [:octicons-arrow-right-24: Open](../ROADMAP/roadmap.md)

</div>
