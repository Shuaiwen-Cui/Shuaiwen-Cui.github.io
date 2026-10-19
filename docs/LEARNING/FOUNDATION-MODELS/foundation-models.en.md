# Foundation Models

Understand what LLMs and multimodal models do, explore models and data through small experiments, then choose a direction for deeper work.

<div class="grid cards" markdown>

-   :material-rocket-launch:{ .lg .middle } __Get Started · Two Weeks__

    ---

    A 14-day loop: a small GPT, LoRA tuning, retrieval QA with sources, quantised and running offline on a Jetson.

    [Start here](GET-STARTED/get-started.md)

-   :material-compass-outline:{ .lg .middle } __Overview__

    ---

    Start with terminology, system diagrams and introductory courses, including Chinese lectures and visual explanations.

    [Start here](OVERVIEW/overview.md)

-   :material-map-marker-path:{ .lg .middle } __Roadmap__

    ---

    Move from orientation to fundamentals, model use, adaptation and specialization, with explicit completion criteria.

    [Follow the roadmap](ROADMAP/roadmap.md)

-   :material-hammer-wrench:{ .lg .middle } __Guided Projects__

    ---

    Follow tutorials from a small PyTorch model to a personal knowledge assistant; make changes and evaluate them.

    [Start building](PROJECT-BASED-LEARNING/project-based-learning.md)

</div>

## Knowledge map

| Unit | Main question | Depth |
|---|---|---|
| [Architecture & representation](ARCHITECTURE/architecture.md) | How do tokens, embeddings and attention form a model? | Fundamentals; architecture variants optional |
| [Training & adaptation](TRAINING/training.md) | What do pretraining, SFT, LoRA and preference optimization change? | Fundamentals; large-scale training optional |
| [Applications & agents](APPLICATIONS/applications.md) | How do models connect to documents, tools and tasks? | Fundamentals and practice |
| [Multimodality & action interfaces](MULTIMODAL/multimodal.md) | How do images, language, state and actions connect? | Orientation and exploration; VLA adaptation advanced |
| [Inference & deployment](INFERENCE/inference.md) | How do quality, latency, throughput and memory interact? | Fundamentals; runtime optimization advanced |
| [World models](WORLD-MODELS/world-models.md) | How can predicted changes support planning or learning? | Advanced elective |

!!! tip "How to use this section"
    Follow two continuous courses: Lee's overview and Karpathy's implementation. Start with F0, use F1 only to fill PyTorch gaps, and use the small GPT in F4 as the implementation project. Choose F2 or F3 for exploration; F5 is a later elective. Consult documentation for specific interfaces.

## Connection to embodied intelligence

- **Skill calling**: a language model selects goals or skills; robot controllers execute them.
- **VLA**: vision, language and robot state condition an action policy, requiring robot data and closed-loop evaluation.
- **World models**: predict consequences of actions for planning or policy learning.

See the [embodied intelligence overview](../EMBODIED-INTELLIGENCE/OVERVIEW/overview.md) for their place in a robot system.
