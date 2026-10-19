# Graph Learning & Geometric Deep Learning

!!! abstract "What this page solves"
    Put relations, topology and symmetry explicitly into the model. It meets your sensor networks (spatial relations between nodes), point clouds and 3D geometry, but graph structure must be justified by the task or data — not used for its own sake. **P1: know it; P2: go deep by direction.**

## Knowledge map

| Level | Topic | What to understand |
|---|---|---|
| P1 | Graph basics and tasks | Nodes, edges, adjacency; node-, edge- and graph-level prediction; permutation invariance and equivariance |
| P1 | Message passing | GCN, GraphSAGE, GAT; neighbourhood aggregation, sampling and features |
| P2 | Modern geometric methods | Graph Transformers, positional encodings, point clouds and SE(3)-equivariant models |
| P2 | Limits and scientific use | Over-smoothing, over-squashing, graph-construction leakage; molecules, sensor graphs, physical fields |

## Core mechanics

**Message passing**: in each layer every node collects its neighbours' features (messages), aggregates them (sum, mean or attention-weighted) and combines with its own to update. After k layers a node "sees" information within k hops.

**Permutation invariance and equivariance**: graphs have no natural node order, so aggregation must be order-independent. Graph-level outputs should be invariant; node-level outputs should reorder along with the nodes (equivariant).

**GCN, GraphSAGE, GAT**: GCN averages with degree normalisation; GraphSAGE samples neighbours and works on nodes unseen in training (inductive); GAT learns per-neighbour weights with attention.

**Over-smoothing and over-squashing**: with many layers all node representations converge (over-smoothing); distant information must squeeze through few paths and is lost (over-squashing). So GNNs are usually shallow.

**Geometric equivariant networks**: properties of point clouds and molecules should not change under global rotation and translation. SE(3)-equivariant networks build the symmetry into the structure, needing less data than learning it through augmentation.

**Transductive vs inductive**: transductive trains and tests on the same graph (only labels are hidden); inductive must generalise to new graphs or nodes. State which one you evaluate, and never use test labels when building the graph.

## Curated resources

<div class="grid cards" markdown>

-   :material-graph-outline:{ .lg .middle } __Main course · Stanford CS224W__

    ---

    Node embeddings, the message-passing paradigm and higher-order graph structure.

    [:octicons-arrow-right-24: Open course page](https://web.stanford.edu/class/cs224w/)

-   :material-code-tags:{ .lg .middle } __Official follow-along · PyG quick start__

    ---

    Cora node classification; reproduce GCN on a small graph (AI-3).

    [:octicons-arrow-right-24: Open docs](https://pytorch-geometric.readthedocs.io/en/latest/get_started/introduction.html)

</div>

## Follow-along · AI-3 · What does graph structure actually add?

1. Run the PyG example; understand x, edge_index, labels and data masks.
2. On the same split, compare an MLP that ignores edges with a GCN.
3. Perturb some edges and record validation and test results; avoid leaking test labels when rebuilding the graph.
4. State whether the example is transductive or inductive, then discuss sensor applications.

**Done when**: starts on a CPU with a small graph; report the baseline, edge perturbation, masks and limits, without equating node classification with structural damage identification.

## Self-check

??? question "Why are GNNs usually only a few layers deep?"
    More layers over-smooth node representations until they are indistinguishable and over-squash distant information. For long-range information use graph Transformers, virtual nodes or skip connections.

??? question "How do you show the graph structure really helps?"
    Compare against an edge-free MLP on the same split, then shuffle or delete some edges and check that performance drops accordingly.

??? question "What to watch when modelling a sensor network as a graph?"
    Edge definitions (physical connection, distance threshold or correlation) need justification; building edges from correlations computed on all data leaks test information into the model.

## Related topics

<div class="grid cards" markdown>

-   __Architectures__

    ---

    GNNs compared with CNN, Transformer and SSM computation patterns.

    [:octicons-arrow-right-24: Open](../ARCHITECTURES/architectures.md)

-   __Time-series & scientific AI__

    ---

    Sensor graphs and physical-field modelling.

    [:octicons-arrow-right-24: Open](../SCIENTIFIC-AI/scientific-ai.md)

</div>
