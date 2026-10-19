# Foundation Models · Two-Week Get Started

!!! abstract "How this page differs from the systematic track"
    The systematic track (overview, roadmap, F0–F5) walks the full path from orientation to implementation, use and adaptation. This page is **the smallest loop you can close in two weeks**: understand a small GPT, LoRA-tune a small model, build a retrieval QA with sources, quantise it and run it offline on an edge device. About 2–3 hours a day.

## Target skills (worked back from common role requirements)

| Skill | Usually phrased as | What you already have | Close the gap to |
|---|---|---|---|
| Model mechanics | Transformer, attention, KV cache | Trained Transformer / Mamba time-series models | Explain decoders, prefill and token-by-token decoding |
| Use & adaptation | Hugging Face, SFT, LoRA | PyTorch training | One LoRA fine-tune, evaluated, on your own |
| Applications | RAG, evaluation, hallucination | Your website is a ready corpus | A QA that shows sources and has a fixed test set |
| Inference | Quantisation, llama.cpp, throughput / latency | **Quantisation and memory work on MCUs** | 4-bit quantisation; how to measure time-to-first-token and tokens/s |

## 14-day plan

| Day | Learn | Do | Output |
|---|---|---|---|
| 1 | The big picture | HF LLM Course ch. 1; the site's Overview | Your own LLM system diagram (training, inference, retrieval, tools) |
| 2–4 | A small GPT | Code along all of Karpathy "Let's build GPT" | Training curve; every tensor shape annotated |
| 5 | Using a pretrained model | Load Qwen2.5-0.5B-Instruct with transformers; chat template and sampling | temperature / top-p comparison |
| 6–7 | LoRA | PEFT: fine-tune on 100–200 QA pairs written from your site | Before / after on a fixed 20-question set |
| 8–9 | RAG | Chunk, embed and retrieve the site's Markdown (LlamaIndex works, see F2) | QA with sources; 20-question evaluation separating retrieval from generation failures |
| 10–11 | Quantisation & inference | llama.cpp: convert to GGUF, compare Q8 and Q4; measure on PC and Jetson | Time-to-first-token, tokens/s, peak memory, quality regressions |
| 12 | Fully offline loop | Point the RAG at the local quantised model | Personal knowledge-base QA running offline on a Jetson |
| 13 | Package the project | README: architecture, results, limits | Repository front page |
| 14 | Interview review | Answer every question below | Answer notes |

## Resources for these two weeks

<div class="grid cards" markdown>

-   :material-open-in-new:{ .lg .middle } __HF LLM Course ch. 1__

    ---

    Day 1 · The big picture

    [:octicons-arrow-right-24: Open](https://huggingface.co/learn/llm-course/en/chapter1/1)

-   :material-book-open-variant:{ .lg .middle } __Overview (systematic track)__

    ---

    Day 1 · The big picture

    [:octicons-arrow-right-24: Open](../OVERVIEW/overview.md)

-   :material-open-in-new:{ .lg .middle } __Karpathy · Zero to Hero__

    ---

    Day 2–4 · A small GPT

    [:octicons-arrow-right-24: Open](https://karpathy.ai/zero-to-hero.html)

-   :material-open-in-new:{ .lg .middle } __Qwen2.5-0.5B-Instruct__

    ---

    Day 5 · Using a pretrained model

    [:octicons-arrow-right-24: Open](https://huggingface.co/Qwen/Qwen2.5-0.5B-Instruct)

-   :material-open-in-new:{ .lg .middle } __Hugging Face PEFT__

    ---

    Day 6–7 · LoRA

    [:octicons-arrow-right-24: Open](https://huggingface.co/docs/peft/index)

-   :material-book-open-variant:{ .lg .middle } __F2 · RAG project__

    ---

    Day 8–9 · RAG

    [:octicons-arrow-right-24: Open](../PROJECT-BASED-LEARNING/project-based-learning.md)

</div>


## Project: a personal knowledge-base QA on an edge device

A LoRA-tuned small model + retrieval + 4-bit quantisation, running offline on a Jetson, answering questions about your research and projects with sources. It covers adaptation, application and deployment at once.

## Must-know interview questions (answer key)

| Question | Key points |
|---|---|
| What are Q, K and V? | Three linear projections of each token; Q·K decides what to attend to, the weights mix V; divide by √d to keep softmax from saturating |
| Why does a KV cache speed things up? | Past tokens' K and V do not change during decoding, so cache them and compute only the new token; memory grows linearly with length |
| Prefill vs decode? | Prefill processes the whole prompt in parallel, compute-bound, sets time-to-first-token; decode generates one token at a time, memory-bandwidth-bound, sets tokens/s |
| What is RoPE? | Rotary position embedding: rotate Q and K by position so their dot product depends only on relative position; extrapolates better |
| What do pretraining, SFT and preference optimisation each change? | Pretraining learns language and knowledge; SFT learns instruction format and behaviour; RLHF / DPO shift outputs toward human preferences |
| How does LoRA work? | Freeze W, learn a low-rank update ΔW = BA (rank r ≪ dimension); few trainable parameters, can be merged back; alpha / r scales it |
| Common quantisation methods? | Mostly weight-only (GPTQ, AWQ, GGUF k-quants); activations are harder because of outliers; 4-bit usually costs acceptable quality |
| How do you locate a RAG failure? | Check whether retrieved chunks contain the answer: if not, it is retrieval (chunking, embeddings, top-k); if yes but wrong, it is generation (prompt, model) |
| How do you evaluate hallucination? | A fixed question set including unanswerable ones; check whether each claim is supported by a source; count unsupported claims |
| Is a longer context always better? | Cost and latency grow with length and attention dilutes; retrieval plus a curated context usually beats filling the window |

## Afterwards

Return to "choose a direction from practice" in the foundation-model roadmap; for inference continue with Inference & Deployment, for robots continue with the embodied Get Started.


<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Overview (systematic track)__

    ---

    Systematic learning

    [:octicons-arrow-right-24: Open](../OVERVIEW/overview.md)

-   :material-book-open-variant:{ .lg .middle } __Foundation-model roadmap__

    ---

    Systematic learning

    [:octicons-arrow-right-24: Open](../ROADMAP/roadmap.md)

-   :material-book-open-variant:{ .lg .middle } __Inference & Deployment__

    ---

    Systematic learning

    [:octicons-arrow-right-24: Open](../INFERENCE/inference.md)

-   :material-book-open-variant:{ .lg .middle } __Embodied Get Started__

    ---

    Systematic learning

    [:octicons-arrow-right-24: Open](../../EMBODIED-INTELLIGENCE/GET-STARTED/get-started.md)

</div>
