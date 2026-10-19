# Data, Evaluation & Trustworthy AI

Evaluation runs through every method: data suitability, task-aligned metrics and reliability under shift. Use this page as an experimental checklist, not another long course.

## Map and priority

| Depth | Topic | Understand |
|---|---|---|
| P0 | Data and design | Provenance, licensing, cleaning, duplicates, missingness, leakage, splits and held-out evaluation |
| P0 | Task metrics | Different metrics for classification, regression, retrieval, generation and control; variance and failures |
| P1 | Uncertainty and robustness | Calibration, confidence, OOD, noise, shift, anomalies and abstention |
| P1 | Interpretation and traceability | Attribution, visualization, ablation and errors; explanation is not causal proof |
| P2 | Frontier reliability | Contamination, hallucination, tool traces, reward gaming, safety/privacy and bias; test the scenario |

## Core mechanistic intuition · Reliability & evaluation

- **Confidence vs. correctness**: Modern deep neural networks with Softmax output are routinely overconfident. Calibration (e.g., Temperature Scaling) measures whether predicted probabilities align with observed empirical frequencies.
- **Out-of-Distribution (OOD) & domain shift**: Identical distribution (i.i.d.) is the fundamental premise of empirical risk minimization. When sensor drift or environmental conditions change, performance can degrade abruptly; density or reconstruction-based rejection filters are essential.
- **Hidden leakage pathways**: Leakage extends far beyond identical sample duplication. Preprocessing statistics computed across the entire corpus, random splits on autocorrelated time series, and cross-split sensor contamination yield illusory benchmark results.

## Selected resources

<div class="grid cards" markdown>

-   :material-table-split:{ .lg .middle } __Validation & Splits · scikit-learn CV__

    ---

    Understand temporal splits (TimeSeriesSplit), grouping (GroupKFold) and stratification to avoid leakage.

    [:octicons-arrow-right-24: View Official CV Guide](https://scikit-learn.org/stable/modules/cross_validation.html)

-   :material-gauge:{ .lg .middle } __Calibration · scikit-learn Probability Calibration__

    ---

    Diagnose deep neural network overconfidence; evaluate reliability diagrams and temperature scaling.

    [:octicons-arrow-right-24: View Official Calibration Guide](https://scikit-learn.org/stable/modules/calibration.html)

</div>

## Guided experiment · AI-7 · Reliability checks on an existing model

1. Select an existing project; specify success, splits and metrics without adding models.
2. Define representative noise/missingness/condition changes and distinguish their origins.
3. Report errors by subgroup/condition; check calibration or generation support as appropriate.
4. Make one improvement and check gains/regressions on reserved examples.

**Acceptance**: Report reproducible conditions, evidence, failures and limits; historical completion checkmarks are not independent project evidence.

## Check yourself · Critical thinking

- [ ] Why is Expected Calibration Error (ECE) frequently more informative than raw accuracy in mission-critical safety or damage detection?
- [ ] How can a lightweight selective classification or uncertainty-based abstention mechanism be deployed to prevent confident failures on anomalous inputs?

Reuse these checks throughout the roadmap; LLM and robot-specific details stay in their sections.

## Related topics

<div class="grid cards" markdown>

-   __roadmap__

    ---

    Three-round route and project pool

    [:octicons-arrow-right-24: Open](../ROADMAP/roadmap.md)

-   __AI overview__

    ---

    The AI knowledge map and how modules relate

    [:octicons-arrow-right-24: Open](../ai.md)

-   __Roadmap__

    ---

    Three-round route and project pool

    [:octicons-arrow-right-24: Open](../ROADMAP/roadmap.md)

</div>
