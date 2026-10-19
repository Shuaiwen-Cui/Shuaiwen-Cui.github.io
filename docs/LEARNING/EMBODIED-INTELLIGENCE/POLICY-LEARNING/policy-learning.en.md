# Learned Policies

!!! abstract "What this page solves"
    Teach a robot to act from data: a policy maps observations to actions. It can imitate demonstrations, be reinforced by rewards, or generalise through a pretrained vision-language model. **Advanced, required**: prerequisites are a PyTorch training loop and the notions of observation and action.

## How far to go

| Level | Learn | Done when |
|---|---|---|
| Orientation | Policies, behaviour cloning, imitation, RL, VLA | You can say where each method's learning signal comes from |
| Foundation | Demonstrations, action representations, ACT, closed-loop evaluation | E3 done; demonstration coverage compared against success rate |
| Advanced | Diffusion policies, RL, VLA and cross-embodiment adaptation | Two methods compared on one fixed task |

## Core knowledge

**Behaviour cloning (BC)**: treat demonstrations as supervised data and regress actions. Simple and effective, but with **compounding error**: small mistakes push the robot into states the demonstrations never covered, where the policy is worse, and errors grow. This is why "low offline loss, closed-loop failure" happens.

**Fighting distribution shift**: DAgger asks the expert to label states the policy actually visits; wider data coverage (perturbations, varied initial states); action chunking (predicting a short sequence at once) reduces the number of decisions.

**ACT (Action Chunking with Transformers)**: predicts the next k actions at once, encoding images and joint states with a Transformer and smoothing execution with temporal ensembling. Data-efficient; a common baseline for bimanual manipulation.

**Diffusion policies**: generate action sequences with a diffusion model and can represent multi-modal distributions (several valid actions in one state, e.g. going around an obstacle left or right); BC with mean-squared error averages the modes into a wrong action.

**Reinforcement learning (RL)**: maximise cumulative reward through interaction. Suited to tasks with clear rewards and abundant simulation (legged locomotion, in-hand rotation); hard parts are reward design, exploration and sample efficiency, so training usually happens in simulation and transfers.

**VLA (vision-language-action models)**: add an action output to a pretrained vision-language model and train on large robot datasets (OpenVLA, π0, SmolVLA). Strength: language instructions and cross-task generalisation. Hard parts: action representation (discrete tokens or continuous generation), data scale and inference latency.

**Hierarchical policies**: a high level (language model or planner) chooses sub-goals or skills; a low level (learned policy or classical controller) executes at high frequency. The most deployed structure today.

**Data sources**: teleoperation (high quality, expensive), cross-embodiment datasets (e.g. Open X-Embodiment), human videos (abundant but without action labels), simulation. Data often matters more than the algorithm.

## Debate · modular vs. end-to-end · imitation vs. reinforcement

End-to-end: one large model from pixels to actions learns behaviours nobody could design, and scale will bring generalisation. Modular: only interpretable, debuggable, separately verifiable modules reach production. Imitation depends on demonstration quality and coverage; RL depends on rewards and simulation fidelity. The common recipe is imitation for an initial policy, then RL or online data to fine-tune.

## Curated resources

<div class="grid cards" markdown>

-   :material-school:{ .lg .middle } __Concepts · LeRobot Robot Learning__

    ---

    One tutorial from imitation learning and data formats to policy interfaces; E3 in the same ecosystem.

    [:octicons-arrow-right-24: Open tutorial](https://huggingface.co/spaces/lerobot/robot-learning-tutorial)

-   :material-robot:{ .lg .middle } __Practice · LeRobot ACT__

    ---

    ACT principles, configuration and training commands, used directly in E3.

    [:octicons-arrow-right-24: Open docs](https://huggingface.co/docs/lerobot/act)

</div>

## Follow-along experiment

<div class="grid cards" markdown>

-   :material-hammer-wrench:{ .lg .middle } __E3 · Imitation learning in simulation__

    ---

    Teleoperate demonstrations → train ACT → fixed-condition closed-loop evaluation → change coverage and compare.

    [:octicons-arrow-right-24: See steps](../PROJECT-BASED-LEARNING/project-based-learning.md)

</div>

## Self-check

??? question "Why can behaviour cloning have low training loss but fail in closed loop?"
    Compounding error causes distribution shift: execution errors lead to states the data never covered, where the policy never learned, and errors grow. Offline loss only measures demonstration states.

??? question "Why are multi-modal action distributions a problem for MSE regression?"
    If one state allows two valid actions (go left, go right), MSE learns their average — straight into the obstacle. Diffusion policies or mixture outputs can represent several modes.

??? question "Why does action chunking help?"
    Predicting several steps at once means fewer decisions and fewer chances for compounding error, and actions are more coherent over time.

??? question "Why is VLA inference latency a deployment bottleneck?"
    One large-model forward pass can take tens to hundreds of milliseconds while control runs at tens to hundreds of hertz. Common fix: a low-frequency large model outputs action chunks or sub-goals, a high-frequency small controller executes.

## Related topics

<div class="grid cards" markdown>

-   __Sim-to-Real__

    ---

    Getting a simulation-trained policy reliably onto the robot.

    [:octicons-arrow-right-24: Open](../SIM-TO-REAL/sim-to-real.md)

-   __Multimodality & Action (Foundation Models)__

    ---

    The model side of VLA: how vision-language models output actions.

    [:octicons-arrow-right-24: Open](../../FOUNDATION-MODELS/MULTIMODAL/multimodal.md)

</div>
