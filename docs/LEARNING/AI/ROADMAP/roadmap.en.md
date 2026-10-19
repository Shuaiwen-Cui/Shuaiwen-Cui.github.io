# AI · Roadmap

!!! tip "Short on time?"
    Do the [two-week Get Started](../GET-STARTED/get-started.md) first, then complete the three rounds below. It covers round two's shared foundations plus one deployment project.

The map provides coverage; this route sets priorities. Reuse your mathematics, signals, geometry and embedded experience, diagnose gaps, implement one model, then explore one or two directions.

## Three passes

| Pass | Content | Resource scope | Completion |
|---|---|---|---|
| 1 · Overview | AI map, topic tables and frontier index | Purpose, interfaces, relationships and limitations; roughly 1–2 weeks of selected reading | Classify terms as architectures, objectives, tasks or systems methods |
| 2 · Foundation and implementation | P0 data, optimization, evaluation and PyTorch gaps | Existing Karpathy/F4 track; skip mastered material; roughly 2–4 weeks | Run training/validation, change one factor and explain outcomes |
| 3 · One or two explorations | Vision, generation, graphs, RL or scientific AI | One official example per chosen direction, not parallel full courses | Use evidence and interest to choose specialization |

Durations are estimates. Full graphical models, topology, SDEs, SLAM and large-scale pretraining are not entry requirements.

## Project pool: choose rather than collect checkmarks

| Project | Entry | Prerequisites/resources | Experience |
|---|---|---|---|
| AI-1 baselines | [Foundations](../FOUNDATIONS/foundations.md) | Python/splits; CPU | Overfitting, metrics and leakage |
| AI-2 transfer | [Representation](../REPRESENTATION/representation.md) | F1-equivalent training; GPU convenient | Frozen/full fine-tuning, labels and shift |
| G1–G3 geometry/generation | [Generative topics](../GENERATIVE-MODELING/generative-modeling.md) | G1 CPU; G2/G3 probability and PyTorch | Geometry, denoising, fields and sampling |
| AI-3 small GCN | [Graphs](../GNN/gnn.md) | Graphs/PyTorch; small CPU | Relations, graph splits and baselines |
| AI-4 CartPole | [Decision learning](../DECISION-LEARNING/decision-learning.md) | MDPs/training; CPU entry | Interaction, exploration and repeated evaluation |
| AI-5 signals/physics | [Scientific AI](../SCIENTIFIC-AI/scientific-ai.md) | Temporal splits; operators/PINNs later | Research transfer, physical checks and generalization |
| AI-6/AI-7 systems/reliability | [Systems](../SYSTEMS/systems.md) / [evaluation](../EVALUATION/evaluation.md) | Reuse any existing model | Performance, reliability and limitations |

**Default suggestion**: shared-gap checks and F4 first; then one of AI-2 or G2/G3, plus one of AI-5 or robot E3. Substitute graphs/RL if preferred. AI-6/AI-7 audit existing work rather than add major projects.

## Core resources & learning budget

<div class="grid cards" markdown>

-   :material-school:{ .lg .middle } __Backbone · Intuition & Scratch Implementation__

    ---

    Lee's lectures establish the global taxonomy; Karpathy Zero to Hero builds a GPT from tensors.

    [:octicons-arrow-right-24: Karpathy Tutorial series](https://karpathy.ai/zero-to-hero.html)

-   :material-eye:{ .lg .middle } __Vision & Representations · CS231n 2026__

    ---

    Read lectures on convolutions, residual connections and self-supervision using public slides.

    [:octicons-arrow-right-24: Stanford CS231n schedule](https://cs231n.stanford.edu/schedule)

-   :material-chart-bell-curve-cumulative:{ .lg .middle } __Geometry & Generation · MIT 6.S184__

    ---

    A unified lecture and lab track covering diffusion, flow matching and latent vector fields without notation conflicts.

    [:octicons-arrow-right-24: Access Course & Labs](https://diffusion.csail.mit.edu/2026/index.html)

-   :material-compass-rose:{ .lg .middle } __Specializations · Graphs / Decision / Science AI__

    ---

    Enter primary tracks (CS224W / David Silver / NeuralOperator) only after choosing a specialization focus.

    [Explore Project Pool](#guided-project-pool-choose-rather-than-complete-every-item)

</div>

## Connect to jobs

After orientation and an implementation, compare 3–5 actual roles using interest, transferable skills, resources and demonstrable results. Improve one main project; vocabulary familiarity is distinct from independent capability.

Use the [foundation roadmap](../../FOUNDATION-MODELS/ROADMAP/roadmap.md) for LLM post-training, RAG, agents and inference, and the [embodied roadmap](../../EMBODIED-INTELLIGENCE/ROADMAP/roadmap.md) for policies/control/transfer.

[AI overview](../ai.md) · [Frontier index](../FRONTIERS/frontiers.md)
