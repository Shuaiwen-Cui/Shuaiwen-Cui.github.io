# Foundation Models · Roadmap

!!! tip "Short on time?"
    Do the [two-week Get Started](../GET-STARTED/get-started.md) first: it compresses stages 0, 2, 3 and F2, F5 into one end-to-end project; then choose a direction below.

## Two backbone courses, one implementation track

**Lee's Chinese overview → Karpathy's Zero to Hero.** Build the system picture before modifying models. Estimates assume 8–10 hours per week; use experimental results to set the pace.

| Stage | Main resource | Scope for your background | Practice and completion |
|---|---|---|---|
| 0 · Overview, 1–2 weeks | [Lee's introduction](https://speech.ee.ntu.edu.tw/~hylee/genai/2024-spring.php) | Lectures 1, 6–10 and 12: systems and limitations | F0: explain training, inference, retrieval and tools |
| 1 · Diagnostic, 1–3 sessions | [PyTorch Quickstart](https://docs.pytorch.org/tutorials/beginner/basics/quickstart_tutorial.html) | Explain training, validation and save/load; study Basics only for gaps | F1 as needed; skip already-mastered Python and mathematics |
| 2 · Implementation, 2–4 weeks | [Zero to Hero](https://karpathy.ai/zero-to-hero.html) | Micrograd if backprop is unclear; first two makemore lectures, then GPT; revisit other lectures as needed | F4 first half: small GPT, curves, tensor shapes and one controlled change |
| 3 · Pretrained models, about 1 week | [HF LLM Course](https://huggingface.co/learn/llm-course/en/chapter1/1) chapter 2 as interface reference | Tokens, loading, context and generation; not a second full course | Use a model and explain differences from your own implementation |
| 4 · Application exploration, 1–2 weeks | One official tutorial in F2 or F3 | Choose document Q&A or image Q&A | Fixed examples with successes and failures |

Read the [embodied overview](../../EMBODIED-INTELLIGENCE/OVERVIEW/overview.md) alongside stage 0. F1 is shared preparation and only needs completing once.

## Choose depth after practice

| Direction | Next step | Resource |
|---|---|---|
| Model algorithms | Finish F4, then adapt to one task | [HF fine-tuning](https://huggingface.co/docs/transformers/training); CS336 later for training research |
| Applications and agents | Extend F2 with one verifiable read-only tool | One project entry in [Applications](../APPLICATIONS/applications.md) |
| Inference systems | F5 quality/latency/memory comparison | llama.cpp in [Inference](../INFERENCE/inference.md) |
| Multimodality and robots | E3, then select E4 | LeRobot in the [embodied roadmap](../../EMBODIED-INTELLIGENCE/ROADMAP/roadmap.md) |

**Text-based alternative:** chapters 2–5 of [LLMs-from-scratch](https://github.com/rasbt/LLMs-from-scratch) can replace Karpathy's implementation track. Choose a medium and stay with its codebase.

MoE, full preference optimization, large-scale parallelism, world models and multi-agent systems remain electives. Topic pages are not additional required courses.

## Learning and career decisions

Complete both overviews and one implementation project, then compare 3–5 target roles. Classify skills as explain, reproduce or independently improve; use project evidence to select the next step.

Run the original example, change one factor and evaluate on examples not used for tuning. Keep versions, results and engineering decisions.

[F0–F5 project steps](../PROJECT-BASED-LEARNING/project-based-learning.md)
