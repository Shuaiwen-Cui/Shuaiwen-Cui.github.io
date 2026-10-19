# Manifolds & Latent Spaces

!!! abstract "What this page solves"
    Why high-dimensional data is often "not really that high-dimensional", and what a learned latent space means. Your geometry and linear-algebra background carries over, and the latent state in your Manifold Dynamics work belongs to this view. **Topic elective**; no differential geometry or topology needed in round one.

## Core mechanics

| Concept | Meaning | Caution |
|---|---|---|
| Ambient space | Where the raw numbers live | A 224×224 colour image is a 150,000-dimensional vector |
| Manifold | A space locally described by low-dimensional Euclidean coordinates | Complex data need not form one smooth manifold |
| Manifold hypothesis | Real data concentrates near low-dimensional structure | A useful modelling assumption to check against data |
| Latent space | Internal variables a model uses or learns | Not guaranteed low-dimensional, interpretable or geometry-preserving |

**The Swiss-roll lesson**: points live in 3D, yet variation along the rolled surface needs only two parameters. Straight-line 3D distance wrongly makes neighbours of points on adjacent layers; distance along the surface (geodesic) reflects true similarity.

**Linear vs non-linear reduction**: PCA finds directions of maximum variance and cannot unroll a curled manifold; Isomap, LLE and UMAP unroll using local neighbourhoods. In t-SNE and UMAP plots, distances between clusters and cluster sizes usually mean nothing.

**Networks and representations**: a deep network can be seen as a chain of coordinate transforms that pulls tangled classes apart layer by layer until they are linearly separable. Good representations make downstream tasks easy.

## Curated resources

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Intuition · Chris Olah__

    ---

    Neural Networks, Manifolds, and Topology: figures of how networks transform representations. Treat conclusions as examples for particular networks and data, not general theorems.

    [:octicons-arrow-right-24: Open article](https://colah.github.io/posts/2014-03-NN-Manifolds-Topology/)

-   :material-chart-scatter-plot:{ .lg .middle } __Follow G1 · scikit-learn Swiss roll__

    ---

    Official downloadable code and notebook; CPU only.

    [:octicons-arrow-right-24: Open example](https://scikit-learn.org/stable/auto_examples/manifold/plot_swissroll.html)

</div>

## G1: unrolling curled data

1. Run the official Swiss-roll example and check the colour-to-point correspondence.
2. State the ambient dimension and the generating parameter dimension.
3. Add PCA as a linear baseline and compare with the chosen non-linear method.
4. Change only the neighbour count or noise; see whether unrolling breaks neighbourhoods, short-circuits or splits.
5. Save the point cloud, 2D embedding and parameters; compare local structure on held-out neighbourhood checks, not just how nice the plot looks.

**Done when**: at least three settings, local neighbourhoods vs global shape explained; a 2D plot does not prove a "true manifold" was found.

## Self-check

??? question "Why can't PCA unroll the Swiss roll?"
    PCA is a linear projection onto the plane of maximum variance; a curled structure needs a non-linear map, and projecting presses layers together.

??? question "What happens when k neighbours is too large?"
    Neighbourhoods span adjacent layers and short-circuit, linking far-apart points, so the unrolled result folds or scrambles.

??? question "Why not judge how different two clusters are from their distance in a t-SNE plot?"
    t-SNE preserves only local neighbourhoods; inter-cluster distance and cluster size depend heavily on parameters and do not reflect the original space.

## Related topics

<div class="grid cards" markdown>

-   __Diffusion Models__

    ---

    Generation in latent space.

    [:octicons-arrow-right-24: Open](../DIFFUSION/diffusion.md)

-   __Representation Learning__

    ---

    Where representations come from and how they transfer.

    [:octicons-arrow-right-24: Open](../../REPRESENTATION/representation.md)

</div>
