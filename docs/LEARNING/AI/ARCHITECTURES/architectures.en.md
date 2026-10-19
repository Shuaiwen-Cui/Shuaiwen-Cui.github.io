# Neural Architectures & Computation

Organize architectures by input structure, information interaction and computation. Architecture, objective and task are separate axes; Transformers can serve language, vision or action generation.

## Map and priority

| Depth | Topic | Understand |
|---|---|---|
| P0 | MLP, CNN and residuals | Feedforward computation, locality, receptive fields, residuals and normalization |
| P1 | Sequence architectures | RNN/LSTM/GRU, Transformers and SSM/Mamba; state compression and context access |
| P1 | Graphs and geometry | Message passing, graph Transformers, symmetries and equivariance |
| P2 | Efficient attention and hybrids | Sparse/linear attention and hybrid layers; FlashAttention is an implementation optimization |
| P2 | MoE and conditional compute | Routing, experts and balancing; total vs. active parameters |

## Core computation motifs comparison

| Family | Inductive Bias | Interaction Mechanism | Complexity (length $L$) | Edge Hardware Affinity |
|---|---|---|---|---|
| **CNN** | Locality & translation equivariance | Sliding window, receptive field grows linearly | Time $O(L)$, Memory $O(L)$ | Very high (weight reuse, structured memory access, NPU friendly) |
| **Transformer** | Minimal bias (permutation equivariant, needs PE) | All-to-all content addressing (dynamic soft routing) | Time $O(L^2)$, Memory $O(L^2)$ | Demands acceleration (KV Cache bandwidth bound during generation) |
| **SSM / Mamba** | Continuous recurrence & state dynamics | Selective continuous compression & hidden update | Training $O(L)$, Step inference $O(1)$ | Good (no expanding KV cache footprint during generation) |
| **GNN** | Topological invariance & permutation equivariance | Message passing along explicit graph edges | Depends on sparsity & sampling | Moderate (irregular memory access and sparse matrix operations) |

## Selected resources

<div class="grid cards" markdown>

-   :material-code-braces:{ .lg .middle } __Implementation · Karpathy Zero to Hero__

    ---

    Follow along a small language model to understand attention, multi-head projections and forward computational graphs.

    [:octicons-arrow-right-24: Access implementation series](https://karpathy.ai/zero-to-hero.html)

-   :material-eye:{ .lg .middle } __Vision architectures · CS231n 2026__

    ---

    Read corresponding lecture slides on convolutions (CNN), residual highways (ResNet) and self-attention.

    [:octicons-arrow-right-24: View official schedule](https://cs231n.stanford.edu/schedule)

</div>

## Guided experiment · Compare structures before implementations

1. Choose your completed F1/F4 model; label inputs, outputs and parameters.
2. Compute key shapes, parameter counts and a forward pass.
3. Replace one component or change context length; record quality and resources.
4. Read one alternative and document assumptions, not parameter count alone.

**Acceptance**: An architecture diagram, controlled change and shape/complexity explanation; MoE training and Mamba CUDA installation are optional.

## Check yourself · Critical thinking

- [ ] Mathematically, how does a residual shortcut ($x + F(x)$) provide an identity gradient highway to preserve backpropagation stability?
- [ ] Why do lightweight modern CNNs (e.g., MobileNet/ConvNeXt) often outperform ViTs on real-time embedded sensor or edge inspection tasks?
- [ ] When context length scales beyond 32k tokens, where are the primary computational and memory bottlenecks for full attention vs. state space models?

Transformer · Mamba / SSM · graphs and geometry; see Foundation Models for large-model training/deployment.

## Related topics

<div class="grid cards" markdown>

-   __Transformer__

    ---

    Tensor skeleton and maths of self-attention

    [:octicons-arrow-right-24: Open](../TRANSFORMER/transformer.md)

-   __Mamba / SSM__

    ---

    State-space sequence models and selectivity

    [:octicons-arrow-right-24: Open](../MAMBA/mamba.md)

-   __graphs and geometry__

    ---

    Message passing, equivariance and graph structure

    [:octicons-arrow-right-24: Open](../GNN/gnn.md)

-   __Foundation Models__

    ---

    Architecture, training, applications and deployment of large models

    [:octicons-arrow-right-24: Open](../../FOUNDATION-MODELS/foundation-models.md)

-   __AI overview__

    ---

    The AI knowledge map and how modules relate

    [:octicons-arrow-right-24: Open](../ai.md)

-   __Roadmap__

    ---

    Three-round route and project pool

    [:octicons-arrow-right-24: Open](../ROADMAP/roadmap.md)

</div>
