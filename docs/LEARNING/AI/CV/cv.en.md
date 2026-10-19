# Vision, Multimodality & 3D Space

!!! abstract "What this page solves"
    Vision tasks extend from "what is it" (classification) to "where" (detection, segmentation), "how far and how big" (depth, 3D) and "how to describe it in language" (vision-language). This is the core of robot perception and many industrial AI roles. **P1: tasks and metrics required; P2: go deep into 3D by direction.**

## Knowledge map

| Level | Topic | What to understand |
|---|---|---|
| P1 | Classification, detection, segmentation | Labels, boxes, pixel masks and their metrics; CNN / ViT, YOLO / DETR, promptable segmentation |
| P1 | Vision-language and video | Image-text alignment, VLMs, temporal context |
| P2 | 3D geometry and representations | Depth, point clouds, pose, occupancy / SDF, NeRF, 3D Gaussian splatting |
| P2 | Space and interaction | Geometric reconstruction, semantics, affordances and action conditioning |

## Core mechanics

**From classification to detection to segmentation**: classification outputs one label; detection outputs boxes plus classes (YOLO: one-stage, fast; Faster R-CNN: two-stage, accurate; DETR: a Transformer predicting the box set directly, no NMS); semantic segmentation classifies every pixel, and instance segmentation also separates individuals of the same class (Mask R-CNN, YOLO-seg).

**Metrics**: IoU is intersection-over-union of prediction and ground truth; mAP integrates the precision-recall curve ranked by confidence, then averages over classes and IoU thresholds. Segmentation often uses mIoU.

**ViT**: cut the image into 16×16 patches as tokens for a Transformer. Weaker than CNNs on little data (no locality prior), stronger with lots; CNNs are still friendlier on edge devices.

**Promptable segmentation (SAM)**: given a point or box, segment the object, with no class-specific training. Suits labelling and new objects in robotics.

**Depth and 3D**: monocular depth (e.g. Depth Anything), stereo, RGB-D; 3D representations include point clouds, meshes, voxels, SDFs; NeRF represents a scene implicitly with a network and renders new views; 3D Gaussian splatting uses many coloured Gaussian ellipsoids explicitly and renders far faster.

**Vision-language**: CLIP aligns image and text spaces; VLMs feed visual tokens into a language model (see the multimodal topic in Foundation Models).

## Curated resources

<div class="grid cards" markdown>

-   :material-school:{ .lg .middle } __Main course (selected) · Stanford CS231n 2026__

    ---

    Classification and transfer first, then detection, segmentation, self-supervision and 3D vision.

    [:octicons-arrow-right-24: Open schedule](https://cs231n.stanford.edu/schedule)

-   :material-target:{ .lg .middle } __Practice · Ultralytics detection & segmentation__

    ---

    YOLO train, validate and export end to end — the fastest way to get detection and instance segmentation running.

    [:octicons-arrow-right-24: Open docs](https://docs.ultralytics.com/tasks/segment/)

-   :material-cube-scan:{ .lg .middle } __3D elective · Nerfstudio Splatfacto__

    ---

    Try 3D Gaussian splatting reconstruction; measure novel-view quality and render latency.

    [:octicons-arrow-right-24: Open docs](https://docs.nerf.studio/nerfology/methods/splat.html)

</div>

## Follow-along

1. Finish AI-2 on the representation page first; do not build a duplicate classification project.
2. Pick one of detection, segmentation or 3D and sketch its inputs, outputs and evaluation.
3. For 3D, run the original Splatfacto example on the official data, holding out novel views from training.
4. Measure novel-view quality, render time and failure regions; separate visual realism from geometric accuracy.

**Done when**: a 3D report states camera, data, compute and evaluation views; a pretty reconstruction alone does not demonstrate spatial understanding. The full robotics perception project is P1 below.

## Self-check

??? question "What separates mAP@0.5 from mAP@0.5:0.95 for detection?"
    The first counts a prediction correct at IoU ≥ 0.5; the second averages over ten thresholds from 0.5 to 0.95 and is much stricter about localisation.

??? question "Why does ViT often lose to CNNs on small datasets?"
    CNNs have built-in locality and translation-invariance priors, so little data suffices; ViT's priors are weak and need large data or strong pretraining to learn these regularities.

??? question "Main difference between NeRF and 3D Gaussian splatting?"
    NeRF is an implicit MLP; rendering queries the network many times per ray, which is slow. 3DGS uses explicit Gaussian ellipsoids with rasterisation, rendering in real time but using more storage.

## Related topics

<div class="grid cards" markdown>

-   __P1 · 6DoF pose estimation (Embodied)__

    ---

    Detection, segmentation and depth to a graspable object pose.

    [:octicons-arrow-right-24: Open](../../EMBODIED-INTELLIGENCE/PROJECT-BASED-LEARNING/pose-6dof.md)

-   __Multimodality & action interfaces (Foundation Models)__

    ---

    Vision-language models and VLA.

    [:octicons-arrow-right-24: Open](../../FOUNDATION-MODELS/MULTIMODAL/multimodal.md)

-   __Representation learning__

    ---

    The AI-2 transfer-learning experiment.

    [:octicons-arrow-right-24: Open](../REPRESENTATION/representation.md)

</div>
