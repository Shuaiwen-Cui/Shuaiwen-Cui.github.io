# Flow Matching

!!! note "Prerequisites and role"
    Use vectors, probability basics, time-dependent functions, numerical ODEs and PyTorch. Begin with 2D distribution transport before large image or robot applications.

## Intuition: learn how points move over time

Choose a probability path from a simple distribution toward data, then learn a time-dependent vector field. Generate by drawing initial points and integrating the learned field with an ODE solver.

A common teaching setup draws a noise endpoint $x_0$ and data endpoint $x_1$, constructs $x_t=(1-t)x_0+t x_1$, and regresses conditional velocity $x_1-x_0$. Here $t=0$ is noise and $t=1$ is data; specify paths, endpoint coupling and time convention.

This is one conditional path choice. Sampling does not know a destination data point, and learned trajectories need not all be straight. **Training can regress a constructible target; generation still requires numerical integration.**

## Related concepts

| Concept | Relationship |
|---|---|
| Diffusion / score matching | Commonly use noise levels and denoising/score objectives; connections to flows exist under appropriate assumptions |
| Flow matching | Commonly fit a velocity field for a probability path; path choices affect learning and sampling |
| Normalizing flow | A family of invertible distribution transformations, not a synonym for the flow-matching objective |
| Rectified flow | Related to specific paths and velocity learning, not a generic name for every flow-matching method |

Speed and quality depend on models, data, solvers and steps. The name alone does not establish faster sampling or replacement of diffusion.

## Core mechanistic intuition · Simulation-free training vs. numerical ODE sampling

- **Training phase (Simulation-free regression)**:
    - Sample noise point $x_0 \sim \mathcal{N}(0, I)$ and clean data point $x_1 \sim q(x)$.
    - Interpolate along a straight conditional path: $x_t = (1-t)x_0 + t x_1$.
    - Target conditional velocity field is constant: $u_t(x|x_1, x_0) = \frac{d x_t}{d t} = x_1 - x_0$.
    - The neural network simply minimizes MSE: $\mathcal{L}_{\text{CFM}}(\theta) = \mathbb{E}_{t, x_0, x_1} \left[ \| v_\theta(x_t, t) - (x_1 - x_0) \|^2 \right]$. **No expensive numerical ODE simulation is required during training**.
- **Sampling phase (ODE integration)**:
    - Sample initial point $x_0 \sim \mathcal{N}(0, I)$.
    - Numerically integrate the ODE $\frac{d x_t}{d t} = v_\theta(x_t, t)$ using Euler, Midpoint, or RK4 solvers from $t=0$ to $t=1$ to synthesize clean sample $x_1$.
- **Contrast with diffusion**: Diffusion models typically evolve along curved, noisy stochastic trajectories. Flow matching with optimal-transport paths yields straighter probability trajectories, enabling high-fidelity generation in substantially fewer numerical solver steps (e.g., 10–20 steps).

## Same course, experiment G3

Use MIT 6.S184, 2026, lectures 1–2 and the flow-matching section of Lab 2, sharing code and notation with G2.

1. Identify noise, data, sampled time and conditional paths.
2. Plot paths, derive the target and train the original 2D example.
3. Visualize vector fields at multiple times and trajectories from fixed initial noise.
4. Vary solver steps; record distributions, trajectories and sampling time.
5. When comparing G2/G3, match data, capacity and evaluation where possible and report training/sampling costs separately.

**Acceptance**: explain information available during training vs. sampling, plus one step/path experiment and failure analysis. State limitations when comparison settings differ.

## Check yourself · Critical thinking

- [ ] Why can training regress on conditional straight-line velocity targets without simulation, while sampling must integrate through the marginalized nonlinear field?
- [ ] Why are the empirical inference trajectories of $v_\theta(x, t)$ not guaranteed to be strictly linear, even when trained on straight conditional paths?
- [ ] In embodied policy learning (Action Flow / Diffusion Policy), why is flow matching better suited to capturing multimodal demonstrations than direct MSE regression?

## Connections to your background

Vector fields and time evolution build on systems/signals intuition. Conditional action generation still needs compatible robot interfaces, frequency and feedback; see Learned policies. Start low-dimensional, then decide whether to pursue VLA.

## Related topics

<div class="grid cards" markdown>

-   __MIT 6.S184, 2026__

    ---

    

    [:octicons-arrow-right-24: Open](https://diffusion.csail.mit.edu/2026/index.html)

-   __Lab 2__

    ---

    MIT 6.S184 lab code

    [:octicons-arrow-right-24: Open](https://github.com/eje24/iap-diffusion-labs/tree/2026)

-   __Learned policies__

    ---

    Imitation, RL and VLA on robots

    [:octicons-arrow-right-24: Open](../../../EMBODIED-INTELLIGENCE/POLICY-LEARNING/policy-learning.md)

-   __Topic map__

    ---

    Generative methods overview and the MIT course route

    [:octicons-arrow-right-24: Open](../generative-modeling.md)

-   __Diffusion__

    ---

    Step-by-step denoising and the score

    [:octicons-arrow-right-24: Open](../DIFFUSION/diffusion.md)

-   __Manifolds and latent spaces__

    ---

    Low-dimensional structure and latent spaces

    [:octicons-arrow-right-24: Open](../MANIFOLDS/manifolds.md)

</div>
