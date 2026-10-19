# Diffusion Models

!!! abstract "What this page solves"
    Generate data by denoising pure noise step by step. Diffusion is the mainstream for image and video generation and the basis of robot diffusion policies. **Topic elective**; read the ideas of noising, denoising and sampling first; coding needs basic probability, gradients and PyTorch.

## Core mechanics

**Forward noising**: add Gaussian noise step by step until, after T steps, the data is nearly pure noise. The noised sample at any t has a closed form: x_t = √ᾱ_t · x_0 + √(1−ᾱ_t) · ε, so training never actually iterates.

**Training objective**: the network takes x_t and t and predicts the added noise ε (equivalently the score, the gradient of the log density), with mean-squared error. One network serves all noise levels, so t is a required input.

**Reverse sampling**: start from noise and remove part of the predicted noise each step, walking back to the data distribution. DDPM is stochastic per step; DDIM is deterministic and allows fewer steps.

**The score view**: the score ∇log p(x) points towards denser data. Diffusion learns the score at many noise scales and follows it (Langevin dynamics).

**Conditional generation**: add class, text, image or robot observations as inputs. Classifier-free guidance trains conditional and unconditional models together and amplifies their difference at sampling, strengthening the condition at the cost of diversity.

**Faster sampling**: deterministic solvers (DDIM, DPM-Solver), consistency models and distillation compress thousands of steps to a few or one.

Diffusion is not a particular U-Net or Transformer; and "stochastic vs deterministic" does not separate diffusion from flow matching — both can sample deterministically.

## Curated resources

<div class="grid cards" markdown>

-   :material-school:{ .lg .middle } __Main course · MIT 6.S184 (2026)__

    ---

    Lecture 1 for sampling intuition, lecture 3-A for score matching; 3-B guidance and 4 architectures when needed.

    [:octicons-arrow-right-24: Open course](https://diffusion.csail.mit.edu/2026/index.html)

-   :material-github:{ .lg .middle } __Labs · same year as the notes__

    ---

    Lab 1 ODE / SDE trajectories, Lab 2 score matching.

    [:octicons-arrow-right-24: Open code](https://github.com/eje24/iap-diffusion-labs/tree/2026)

</div>

## G2: denoising and generating 2D data

1. Do Lab 1's ODE / SDE trajectory exercises first to see step size and random perturbation.
2. Enter the score-matching part of Lab 2 and run the original 2D example.
3. Inspect data and predictions at different noise levels; identify the training input and supervision.
4. Record the training loss, then sample from noise and compare generated and training point clouds.
5. Change only the step count or noise setting, keeping data, training and seeds traceable.

**Done when**: training record, plots and timing for several sampling settings, failure modes. Discuss loss, visual similarity and sample quality separately; 2D results do not speak for images or robots.

**Compute**: start low-dimensional and small; image Lab 3 needs a GPU per the official notebook.

## Self-check

??? question "Why can training compute the noised sample at any t directly?"
    Repeatedly adding Gaussian noise stays Gaussian with closed-form mean and variance, so x_t follows from x_0 and one noise sample.

??? question "What happens with a larger classifier-free guidance scale?"
    Samples follow the condition (e.g. the text) more closely but lose diversity and may over-saturate or distort.

??? question "How do noise prediction and score prediction relate?"
    They differ only by a noise-level-dependent scale and sign: score ≈ −ε / √(1−ᾱ_t).

## Related topics

<div class="grid cards" markdown>

-   __Flow Matching__

    ---

    Another way to transport noise to data, with straighter paths.

    [:octicons-arrow-right-24: Open](../FLOW-MATCHING/flow-matching.md)

-   __Embodied E3 · closed-loop policies__

    ---

    Diffusion for actions still needs matching states, control rates and closed-loop validation.

    [:octicons-arrow-right-24: Open](../../../EMBODIED-INTELLIGENCE/PROJECT-BASED-LEARNING/project-based-learning.md)

</div>
