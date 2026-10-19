# Geometric Representation & Generative Modelling

!!! abstract "What this page solves"
    How models generate new data, and what geometric structure data itself has. Manifolds and latent spaces describe data and representations; diffusion and flow matching are generation methods — and the mathematical basis of robot diffusion policies and π0-style VLA action heads. **Topic elective**: not a shared prerequisite for LLM or robotics entry.

## Topics

<div class="grid cards" markdown>

-   :material-vector-curve:{ .lg .middle } __Manifolds & latent spaces__

    ---

    What low-dimensional structure hides in high-dimensional data? Ambient space, intrinsic dimension, latent variables.

    [:octicons-arrow-right-24: Open topic](MANIFOLDS/manifolds.md)

-   :material-chart-scatter-plot:{ .lg .middle } __Diffusion models__

    ---

    Learn denoising and the score across noise levels, then sample step by step.

    [:octicons-arrow-right-24: Open topic](DIFFUSION/diffusion.md)

-   :material-arrow-decision:{ .lg .middle } __Flow matching__

    ---

    Learn a velocity field from noise to data and generate by ODE integration.

    [:octicons-arrow-right-24: Open topic](FLOW-MATCHING/flow-matching.md)

</div>

## Core mechanics · the landscape of generative methods

Every generative model answers one question: how to turn samples from a simple distribution (usually Gaussian noise) into samples from the data distribution. They differ in how they transport:

| Family | Core idea | Strength | Cost |
|---|---|---|---|
| Autoregressive | Factor the joint into conditionals, generate in order | Simple training, exact likelihood | Slow, sequential |
| VAE | Encode to latents and decode, optimise the ELBO | Explicit latent space | Blurrier samples |
| GAN | A generator fools a discriminator | Sharp, fast samples | Unstable training, mode collapse |
| Normalizing flow | Invertible maps with exact density change | Exact likelihood | Architectures constrained by invertibility |
| Diffusion / score | Learn step-by-step denoising (the log-density gradient) | High quality, stable, multi-modal | Many sampling steps |
| Flow matching | Regress the velocity field from noise to data | Simpler training, straighter paths, fewer steps possible | Still needs numerical integration |

A Transformer is an architecture; diffusion and flow matching are generation methods — a Transformer can be their denoising network (e.g. DiT). The manifold hypothesis is a way of viewing data, not a generation algorithm.

## Curated resources

<div class="grid cards" markdown>

-   :material-school:{ .lg .middle } __Main course · MIT 6.S184 (2026)__

    ---

    One set of notes, videos and labs from ODE/SDE through score matching, flow matching and conditional generation, in one notation.

    [:octicons-arrow-right-24: Open course](https://diffusion.csail.mit.edu/2026/index.html)

-   :material-github:{ .lg .middle } __Labs · MIT 6.S184 code__

    ---

    Lab 1: ODE/SDE; Lab 2: flow / score matching (G2, G3); Lab 3: images and VAE.

    [:octicons-arrow-right-24: Open code](https://github.com/eje24/iap-diffusion-labs/tree/2026)

-   :material-eye:{ .lg .middle } __Supplement · CS231n generative lectures__

    ---

    VAE / GAN basics, only as a supplement to the same map.

    [:octicons-arrow-right-24: Open schedule](https://cs231n.stanford.edu/schedule)

</div>

**For your background**: linear algebra, signals, systems and geometry help with vector fields and time evolution; probability, stochastic processes and PyTorch still need checking through examples. The course expects linear algebra, multivariable calculus, basic probability and PyTorch — not a zero-prerequisite course.

| Stage | Do only this | Done when |
|---|---|---|
| Orientation | The comparison above, lecture 1 figures | You can separate representation, architecture and generation |
| Geometry | G1: Swiss roll | Ambient dimension, intrinsic dimension and embedding distinguished |
| Numerics | Lab 1: ODE / SDE, trajectories before derivations | Step size, determinism and randomness understood |
| Learning & sampling | Lab 2: G2 and G3 | Training objective vs sampling process; 2D samples compared |
| Deeper | Lectures 3-B / 4 and Lab 3: conditioning, DiT, VAE | One image experiment in the same codebase |

**Read concepts any time; code after F1 or equivalent PyTorch.**

## Three follow-along experiments

<div class="grid cards" markdown>

-   :material-cube-scan:{ .lg .middle } __G1 · Data geometry and unrolling__

    ---

    scikit-learn Swiss roll on a CPU: linear PCA vs non-linear neighbourhood methods.

    [:octicons-arrow-right-24: Open example](https://scikit-learn.org/stable/auto_examples/manifold/plot_swissroll.html)

-   :material-wave-square:{ .lg .middle } __G2 · 2D denoising and score matching__

    ---

    The score-matching part of MIT Lab 2; watch gradient fields at different noise levels.

    [:octicons-arrow-right-24: Open lab code](https://github.com/eje24/iap-diffusion-labs/tree/2026)

-   :material-ray-start-arrow:{ .lg .middle } __G3 · Vector-field generation and ODE sampling__

    ---

    The flow-matching part of MIT Lab 2; test how Euler step count affects quality.

    [:octicons-arrow-right-24: Open lab code](https://github.com/eje24/iap-diffusion-labs/tree/2026)

</div>

Run the original first; change one factor per experiment.

## Self-check

??? question "What are the main problems of diffusion models and GANs?"
    GANs train unstably and mode-collapse (generating only a few kinds of samples); diffusion trains stably and covers the distribution but needs tens to thousands of sampling steps.

??? question "Why do robot action generators use diffusion or flow matching instead of regression?"
    One state may allow several valid actions (a multi-modal distribution); regression outputs their average, while diffusion / flow matching samples one valid mode.

??? question "Why does latent diffusion (e.g. Stable Diffusion) encode to a latent space first?"
    Pixel space is high-dimensional and expensive; a VAE compresses to a small latent space, diffusion runs there, and the decoder maps back — far less compute.

## Related topics

<div class="grid cards" markdown>

-   __Foundation Models__

    ---

    Using, training and deploying models.

    [:octicons-arrow-right-24: Open](../../FOUNDATION-MODELS/foundation-models.md)

-   __Learned Policies (Embodied)__

    ---

    Diffusion policies and flow-matching action heads in closed loop.

    [:octicons-arrow-right-24: Open](../../EMBODIED-INTELLIGENCE/POLICY-LEARNING/policy-learning.md)

</div>
