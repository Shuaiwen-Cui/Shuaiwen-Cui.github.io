# Architecture & Representation

!!! abstract "What this page solves"
    What is inside a large model: how text becomes vectors, how attention lets tokens look at each other, and how blocks stack into a GPT. With this page, training, inference and multimodality share one skeleton. **Required**; prerequisites: matrix products, PyTorch tensors.

## How far to go

| Level | Learn | Done when |
|---|---|---|
| Orientation | Tokens, embeddings, context, attention, generation | Input-to-output drawn; F0 done |
| Foundation | Residuals, normalisation, positions, the training loop | Small GPT built in F4; every tensor shape explained |
| Advanced | Long context, MoE, state-space models, linear attention | A comparison experiment designed for one question |

## Core knowledge

**Tokenisation**: BPE and similar algorithms split text into sub-word units. Vocabulary size sets the embedding size and sequence length; how Chinese, code and numbers are split affects cost and capability (digits split badly make arithmetic harder).

**Embeddings**: each token looks up a vector; the output layer maps vectors back to vocabulary probabilities. Many models tie input embeddings and output weights.

**Self-attention**: each token produces Q, K, V; softmax(QKᵀ/√d)·V lets each position gather information from others by content similarity. A causal mask keeps generation from seeing the future. Multiple heads attend to different relations.

**Transformer block**: attention + feed-forward network (FFN, most of the parameters), each wrapped in a residual connection and normalisation. Modern LLMs use pre-norm (RMSNorm) and SwiGLU for stabler training.

**Positional encoding**: attention itself is order-blind. Absolute encodings add a vector per position; RoPE rotates Q and K by position so attention depends only on relative distance — today's default.

**Long context**: attention is O(n²) and the KV cache grows linearly. Remedies: sliding windows, sparse attention, FlashAttention (same complexity, far less memory traffic), RoPE extrapolation and interpolation.

**Mixture of experts (MoE)**: replace the FFN with many experts and let a router activate a few per token. Large total parameters, little compute per token; hard parts are load balancing and memory (all experts must fit).

**Alternatives**: state-space models (Mamba) and linear attention compress history into a fixed-size state, O(1) per inference step and cheaper on long sequences, but weaker at exact long-range recall; hybrids are spreading.

**Unified multimodal representation**: images become patches, audio becomes frames, encoders turn them into tokens of the same width as text, and everything enters one Transformer.

## Debate · autoregressive vs. diffusion

Autoregressive generation, token by token, is simple to train and natural for language — the LLM default. Diffusion denoises in parallel and dominates images, video and robot actions; diffusion language models are explored for parallel decoding. The trade-off is the simplicity of sequential generation against the efficiency and multi-modality of parallel generation.

## Curated resources

<div class="grid cards" markdown>

-   :material-code-braces:{ .lg .middle } __Main track · Karpathy Zero to Hero__

    ---

    First two makemore lectures → "Let's build GPT"; micrograd first if gradients are unclear.

    [:octicons-arrow-right-24: Open course](https://karpathy.ai/zero-to-hero.html)

-   :material-book-open-page-variant:{ .lg .middle } __Text alternative · LLMs-from-scratch__

    ---

    Chapters 2–5 as notebooks, from tokenisation to GPT pretraining; choose this or the previous one.

    [:octicons-arrow-right-24: Open repository](https://github.com/rasbt/LLMs-from-scratch)

-   :hugging:{ .lg .middle } __Interface reference · HF LLM Course__

    ---

    Chapter 2 on tokenizers and model loading, for when an interface question comes up.

    [:octicons-arrow-right-24: Open course](https://huggingface.co/learn/llm-course/en/chapter1/1)

</div>

## Follow-along experiments

<div class="grid cards" markdown>

-   :material-hammer-wrench:{ .lg .middle } __F0 · Tokens and context__

    ---

    Compare Chinese and English tokenisation; change only the context and watch the answer.

    [:octicons-arrow-right-24: See steps](../PROJECT-BASED-LEARNING/project-based-learning.md)

-   :material-hammer-wrench:{ .lg .middle } __F4 · A small GPT__

    ---

    Code along a GPT, then change context length or width.

    [:octicons-arrow-right-24: See steps](../PROJECT-BASED-LEARNING/project-based-learning.md)

</div>

## Self-check

??? question "Why divide attention scores by √d?"
    The variance of Q·K grows with d; large values push softmax towards one-hot and kill gradients. Dividing by √d keeps the variance near 1.

??? question "How is the causal mask implemented?"
    Set the upper triangle (future positions) of the score matrix to −∞ before softmax so their weights become 0.

??? question "Where are most parameters in a Transformer block?"
    In the FFN: two linear layers with an inner width around 4× the hidden size (SwiGLU ≈ 8/3× with three matrices) — about twice the attention projections.

??? question "Why is MoE large in parameters but cheap at inference? And memory?"
    Each token visits only a few experts, so compute follows active parameters; but all expert weights must be loaded, so memory follows total parameters.

## Related topics

<div class="grid cards" markdown>

-   __Transformer (AI)__

    ---

    Tensor skeleton and maths.

    [:octicons-arrow-right-24: Open](../../AI/TRANSFORMER/transformer.md)

-   __Training & Adaptation__

    ---

    How the same architecture is trained and adapted.

    [:octicons-arrow-right-24: Open](../TRAINING/training.md)

</div>
