# Simulation & Sim-to-Real

!!! abstract "What this page solves"
    Success in simulation is not success on the robot. This page covers narrowing the gap between simulation and reality, deploying on constrained hardware, and validating in stages. **Advanced**: go deep once E3 reproduces reliably; your digital-twin-to-field work (adaptive triggering) and edge deployment experience apply directly.

## How far to go

| Level | Learn | Done when |
|---|---|---|
| Orientation | Simulation, software-in-the-loop, hardware-in-the-loop, real robot | You can say what each kind of experiment validates |
| Foundation | Fixed evaluation conditions, control rate, latency | E3 repeated with offline error and closed-loop success reported separately |
| Advanced | Domain randomisation, system identification, edge deployment, safety and recovery | E5 before/after comparison with the domain gap reported |

## Core knowledge

**The reality gap**: simulation and reality differ in dynamics (friction, mass, delay), sensing (lighting, noise, camera parameters) and actuation (motor response, backlash); policies exploit loopholes that only exist in simulation.

**Domain randomisation**: randomise friction, mass, lighting, textures, camera pose and more during training so reality looks like one more sample. Simple and effective; costs are more conservative policies and longer training.

**System identification**: estimate simulator parameters from real data so the simulation behaves like the robot. Complements randomisation: identify first to shrink the gap, then randomise around the identified values. It is the same kind of problem as your structural modal identification.

**Differentiable simulation**: gradients through the simulator let you fit real trajectories or optimise policies directly; gradients are unreliable across contact discontinuities.

**Staged validation**: software-in-the-loop (algorithm + simulation) → hardware-in-the-loop (real controller / compute + simulation) → controlled real robot → open environment. Change one factor per stage so failures can be attributed.

**Real-time and resource budgets**: control rate, end-to-end latency, compute and power together bound model size. Quantisation, distillation and kernel optimisation decide whether a policy runs on time on the onboard computer.

**Safety and failure handling**: action limits, speed limits, emergency stop, watchdog timeouts, a safe posture after detected anomalies — all in place before the real robot.

**Why evaluation is hard**: real trials are slow, initial conditions are hard to reproduce and environments change. Fix a set of initial states, record every failure type, and report "successes / trials" rather than a bare percentage.

## Debate · real data vs. simulated data

Real-data camp: simulation is never real enough; teleoperated real data is most reliable and scale will bring the cost down. Simulation camp: only simulation provides the data volume and safety RL needs, and rendering and physics keep improving. In practice locomotion leans on simulation + RL, fine manipulation on real demonstrations, and mixed training is becoming the norm.

## Curated resources

<div class="grid cards" markdown>

-   :material-cube-scan:{ .lg .middle } __Simulation · MuJoCo docs__

    ---

    What time steps, actuators and contact parameters mean; read before randomising them.

    [:octicons-arrow-right-24: Open docs](https://mujoco.readthedocs.io/)

-   :material-chip:{ .lg .middle } __Edge deployment · Jetson tutorials__

    ---

    NVIDIA's Jetson getting-started and inference deployment tutorials; measure latency with TensorRT.

    [:octicons-arrow-right-24: Open tutorials](https://developer.nvidia.com/embedded/learn/tutorials)

</div>

## Follow-along experiment

<div class="grid cards" markdown>

-   :material-hammer-wrench:{ .lg .middle } __E5 · Deployment and transfer__

    ---

    Fix model and task; measure software-loop latency, then hardware-in-the-loop, then a controlled real robot.

    [:octicons-arrow-right-24: See steps](../PROJECT-BASED-LEARNING/project-based-learning.md)

</div>

## Self-check

??? question "What happens if the randomisation range is too wide?"
    The policy must work in every extreme, becomes overly conservative, performs worse and trains less easily. Cover the real parameter uncertainty, not as much as possible.

??? question "How do system identification and domain randomisation relate?"
    They complement each other: identification pulls simulator parameters to the real centre, randomisation covers the remaining uncertainty; identify first, then randomise around it.

??? question "What does hardware-in-the-loop validate, and what not?"
    It validates latency, compute and interfaces on the real computer; not real physical interaction (contact, friction) or real sensor noise.

??? question "Why report successes / trials on the real robot instead of one success rate?"
    Real trials are few; 16 of 20 and 160 of 200 carry different confidence. Also report initial conditions and failure types.

## Related topics

<div class="grid cards" markdown>

-   __Learned Policies__

    ---

    Where the policy to transfer comes from.

    [:octicons-arrow-right-24: Open](../POLICY-LEARNING/policy-learning.md)

-   __Inference & Deployment (Foundation Models)__

    ---

    Quantisation, latency measurement and edge runtimes.

    [:octicons-arrow-right-24: Open](../../FOUNDATION-MODELS/INFERENCE/inference.md)

</div>
