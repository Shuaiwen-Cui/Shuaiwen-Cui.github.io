# Multimodality & Action Interfaces

!!! abstract "What this page solves"
    How models see images and hear sound, and how they go from describing to acting: how vision-language models (VLMs) connect to robot actions (VLAs). This is where foundation models meet embodied intelligence. **Introduction required, VLA adaptation elective**.

## How far to go

| Level | Learn | Done when |
|---|---|---|
| Orientation | Image–text alignment, VLM, VLA inputs and outputs | You can state the difference between a VLM and a VLA |
| Foundation | Run a VLM, probe spatial and counting ability | F3 done with at least three error types recorded |
| Advanced | Action representations, robot data adaptation, closed-loop evaluation | One route of embodied E4 done |

## Core knowledge

**Contrastive image–text alignment (CLIP)**: an image encoder and a text encoder map matching pairs into one space, pulling matches together and pushing mismatches apart. Result: zero-shot classification and retrieval, and the vision encoder of many VLMs.

**VLM structure**: a vision encoder (ViT) turns the image into patch vectors, a projection maps them into the language model's embedding space, and they enter the LLM as "image tokens" alongside text.

**VLM limits**: good at description and common-sense QA, often wrong on exact counting, fine spatial relations, metric quantities (distance, size) and 3D reasoning. Verify spatial outputs separately before using them on a robot.

**From VLM to VLA**: add an action output and train on robot data (image + instruction + action). Actions can be discretised into tokens generated like words (RT-2, OpenVLA) or produced by a continuous head (diffusion or flow matching, e.g. π0) as smooth action chunks.

**Action representations**: joint angles or end-effector pose, absolute or delta, single step or chunk. The choice decides cross-robot data reuse and control smoothness.

**Cross-embodiment and data**: robots differ in joints, camera placement and control rate. Datasets like Open X-Embodiment train across bodies, but adapting to your robot still needs its own data.

**Two ways to connect a robot**: an LLM chooses skills and parameters for existing controllers (interpretable, easy to verify), or a VLA outputs actions directly (more generalisation potential, harder to debug).

## Debate · discrete action tokens vs. continuous action generation

Discrete: treat actions as vocabulary tokens and reuse language-model training and inference unchanged. Continuous: discretisation loses precision; high-frequency fine manipulation needs smooth continuous chunks. High-performing VLAs are moving to continuous heads; discrete is the simplest to implement.

## Curated resources

<div class="grid cards" markdown>

-   :hugging:{ .lg .middle } __VLM practice · Image-text-to-text__

    ---

    Hugging Face's official image–text QA task guide, used in F3.

    [:octicons-arrow-right-24: Open guide](https://huggingface.co/docs/transformers/main/tasks/image_text_to_text)

-   :material-robot:{ .lg .middle } __VLA practice · SmolVLA__

    ---

    A small VLA in the LeRobot ecosystem, from data to fine-tuning and evaluation.

    [:octicons-arrow-right-24: Open tutorial](https://huggingface.co/docs/lerobot/smolvla)

</div>

## Follow-along experiments

<div class="grid cards" markdown>

-   :material-hammer-wrench:{ .lg .middle } __F3 · VLM abilities and limits__

    ---

    Ten fixed images; ask about objects, counts and spatial relations; record error types.

    [:octicons-arrow-right-24: See steps](../PROJECT-BASED-LEARNING/project-based-learning.md)

-   :material-hammer-wrench:{ .lg .middle } __E4 · Language and action__

    ---

    Skill calling or VLA adaptation, with closed-loop validation.

    [:octicons-arrow-right-24: See steps](../../EMBODIED-INTELLIGENCE/PROJECT-BASED-LEARNING/project-based-learning.md)

</div>

## Self-check

??? question "How does CLIP differ from a VLM?"
    CLIP only maps images and text into one space for matching and cannot generate text; a VLM feeds image features into a language model and answers open questions about the image.

??? question "Why can't a VLM's "the cup is left of the plate" drive a grasp directly?"
    The description is qualitative, may be wrong and gives no metric pose; grasping needs millimetre-level 3D position and orientation from, or verified by, geometric perception.

??? question "Why do VLAs often predict action chunks instead of single steps?"
    Large-model inference is slow; single-step prediction cannot keep up with control rates, while a chunk covers several control cycles per inference and is more coherent.

??? question "What does using an LLM as a skill selector buy you?"
    Actions run on verified controllers, so they stay safe and controllable; the language model only interprets instructions and composes skills, and failures are easy to locate.

## Related topics

<div class="grid cards" markdown>

-   __Learned Policies (Embodied)__

    ---

    Where VLAs sit among policy-learning methods.

    [:octicons-arrow-right-24: Open](../../EMBODIED-INTELLIGENCE/POLICY-LEARNING/policy-learning.md)

-   __World Models__

    ---

    Predicting how the environment changes, rather than outputting actions.

    [:octicons-arrow-right-24: Open](../WORLD-MODELS/world-models.md)

</div>
