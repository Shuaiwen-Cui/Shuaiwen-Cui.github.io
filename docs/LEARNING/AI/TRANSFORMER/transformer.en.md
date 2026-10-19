# Transformer

Start tokens, vectors, attention, positions and residuals, then one small implementation. A research paper is a later reference rather than the full prerequisite course.

| Depth | Topic | Completion |
|---|---|---|
| P1 | Attention, Q/K/V, positions, residuals and normalization | Trace data and tensor shapes |
| P2 | ViT, cross-modal, long-context, sparse/linear attention | Define tasks/resources before comparison |
| P2 | MoE and efficient implementation | Separate conditional compute, architecture and kernels |

## Selected resources

<div class="grid cards" markdown>

-   :material-code-braces:{ .lg .middle } __Implementation · Karpathy Zero to Hero__

    ---

    Code-first backbone from makemore to GPT; implement embeddings, causal attention and residual streams from scratch.

    [:octicons-arrow-right-24: Access Tutorial Series](https://karpathy.ai/zero-to-hero.html)

-   :material-video:{ .lg .middle } __Visual Intuition · 3Blue1Brown Attention__

    ---

    Premier geometric visualization explaining Query, Key, Value vectors and multi-head projection spaces.

    [:octicons-arrow-right-24: Watch Video Guide](https://www.youtube.com/watch?v=wjZofJX0v4M)

</div>

## Core mechanistic intuition · Soft routing & tensor flow

- **Intuition**: Self-attention is fundamentally **content-based dynamic weighted information retrieval (soft routing)**. Each token projects a Query (what information am I seeking?) and a Key (what features do I contain?), computing attention affinities via inner products to aggregate Values (what information payload do I carry?).
- **Minimal tensor pipeline**:
    - Input embedding: $X \in \mathbb{R}^{B \times T \times D}$
    - Linear projections: $Q = X W_Q,\; K = X W_K,\; V = X W_V \in \mathbb{R}^{B \times T \times D}$
    - Attention map: $A = \text{Softmax}\left(\frac{Q K^T}{\sqrt{d_k}}\right) \in \mathbb{R}^{B \times T \times T}$ (explicit all-to-all dependency graph)
    - Value aggregation & output projection: $O = A V W_O \in \mathbb{R}^{B \times T \times D}$
- **Trade-offs**: Capturing unrestricted global dependencies requires $O(T^2)$ pairwise interactions. In autoregressive inference, the memory bottleneck shifts from compute to **KV Cache** memory footprint and memory bandwidth limits.

## Check yourself · Critical thinking

- [ ] What is the mathematical purpose of the scaling factor $\sqrt{d_k}$? What happens to the Softmax gradient when $d_k$ is large without scaling?
- [ ] How does a causal mask operate at the matrix arithmetic level, and why is an architectural mask necessary rather than relying solely on the loss function?
- [ ] Why does Pre-LN (normalizing before the sub-layer) offer significantly smoother gradient flow in deep networks compared to original Post-LN?

Revisit papers and historical bookmarks after implementing the small model.

## Related topics

<div class="grid cards" markdown>

-   __Architecture map__

    ---

    CNN, Transformer, SSM and GNN compared

    [:octicons-arrow-right-24: Open](../ARCHITECTURES/architectures.md)

-   __AI overview__

    ---

    The AI knowledge map and how modules relate

    [:octicons-arrow-right-24: Open](../ai.md)

</div>

<details markdown="1">
<summary>Earlier reading bookmarks (optional)</summary>

![](transformer.png)

## __GET STARTED__

<div class="grid cards" markdown>

-   :material-file:{ .lg .middle } __Attention is all you need (Paper)__ 

    ---


    [:octicons-arrow-right-24: <a href="https://arxiv.org/abs/1706.03762" target="_blank"> Portal </a>](#)

-   :material-youtube:{ .lg .middle } __Attention is all you need (Video)✅__ 

    ---


    [:octicons-arrow-right-24: <a href="https://www.youtube.com/watch?v=nzqlFIcCSWQ" target="_blank"> Portal </a>](#)

</div>


## __DIVING DEEPER__

<div class="grid cards" markdown>

-   :material-youtube:{ .lg .middle } __Illustrated Guide to Transformers Neural Network: A step by step explanation🏆__ 

    ---

    A really good one, concise and clear!

    [:octicons-arrow-right-24: <a href="https://www.youtube.com/watch?v=4Bdc55j80l8" target="_blank"> Portal </a>](#)

-   :material-youtube:{ .lg .middle } __Transformers (how LLMs work) explained visually | DL5🏆__ 

    ---

    

    [:octicons-arrow-right-24: <a href="https://www.youtube.com/watch?v=wjZofJX0v4M" target="_blank"> Portal </a>](#)

-   :material-youtube:{ .lg .middle } __Attention in transformers, step-by-step | DL6🏆__ 

    ---

    

    [:octicons-arrow-right-24: <a href="https://www.youtube.com/watch?v=eMlx5fFNoYc" target="_blank"> Portal </a>](#)

-   :material-youtube:{ .lg .middle } __How might LLMs store facts | DL7🏆__ 

    ---

    

    [:octicons-arrow-right-24: <a href="https://www.youtube.com/watch?v=9-Jl0dxWQs8" target="_blank"> Portal </a>](#)


-   :fontawesome-brands-bilibili:{ .lg .middle } __Attention is all you need (Explanation CN)__ 

    ---


    [:octicons-arrow-right-24: <a href="https://www.bilibili.com/video/BV14m421u7EM/?spm_id_from=333.337.search-card.all.click&vd_source=5a427660f0337fedc22d4803661d493f" target="_blank"> Portal </a>](#)


</div>

## __PRACTICE__

<div class="grid cards" markdown>

-   :material-file:{ .lg .middle } __annotated-transformer__ 

    ---


    [:octicons-arrow-right-24: <a href="http://nlp.seas.harvard.edu/annotated-transformer" target="_blank"> Notes </a>](#)

    [:octicons-arrow-right-24: <a href="https://github.com/harvardnlp/annotated-transformer" target="_blank"> Codes </a>](#)

-   :material-file:{ .lg .middle } __transformer tutorial🎯🏆__ 

    ---

    [:octicons-arrow-right-24: <a href="https://medium.com/towards-data-science/all-you-need-to-know-about-attention-and-transformers-in-depth-understanding-part-1-552f0b41d021" target="_blank"> Theory 1 </a>](#)

    [:octicons-arrow-right-24: <a href="https://medium.com/towards-data-science/all-you-need-to-know-about-attention-and-transformers-in-depth-understanding-part-2-bf2403804ada" target="_blank"> Theory 2 </a>](#)

    [:octicons-arrow-right-24: <a href="https://medium.com/towards-data-science/build-your-own-transformer-from-scratch-using-pytorch-84c850470dcb" target="_blank"> Code Tutorial </a>](#)

</div>

</details>
