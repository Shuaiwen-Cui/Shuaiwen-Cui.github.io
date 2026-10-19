# Representation Learning & Learning Paradigms

!!! abstract "What this page solves"
    How features are learned and where the supervision signal comes from. Most of modern AI's capability comes from "learn a good representation on lots of unlabelled data first, then transfer to a specific task". **P0: transfer learning is required; P1: know self-supervision.**

## Knowledge map

| Level | Topic | What to understand |
|---|---|---|
| P0 | Representations, embeddings, transfer | Features, similarity, frozen encoder, linear probe, fine-tuning, domain shift |
| P1 | Self-supervised & contrastive learning | Augmented views, positives/negatives, representation collapse |
| P1 | Masked, teacher and predictive representations | How MAE, DINO and JEPA objectives differ |
| P2 | Continual, active, few-shot, federated learning | Forgetting, labelling budget, distributed privacy; pick by task |

## Core mechanics

**Representation and transfer**: a good representation puts similar inputs close in vector space, so a downstream task needs only a simple head on top. A linear probe (freeze the encoder, train one linear layer) measures representation quality; full fine-tuning usually does better but needs more data and compute and overfits more easily.

**Contrastive learning**: two augmented views of the same image are positives to pull together; different images are negatives to push apart (SimCLR). CLIP applies the same idea to image-text pairs. The key risk is representation collapse — every input maps to one point and the loss can still be low.

**Masked prediction**: hide part of the input and make the model fill it in. BERT masks words, MAE masks image patches (still learns well with 75% masked). No negatives needed.

**Teacher-student and self-distillation**: DINO makes a student match a momentum-updated teacher's outputs across views; no negatives, no collapse, and the features segment semantically well.

**Predictive representations (JEPA)**: predict the masked part in representation space rather than pixel space, so capacity is not wasted reconstructing irrelevant detail — close to the world-model idea.

**Domain shift**: when training and deployment distributions differ (lighting, devices, scenes), transfer breaks. Moving structural monitoring from lab to field is a textbook domain shift.

## Curated resources

<div class="grid cards" markdown>

-   :material-book-open-page-variant:{ .lg .middle } __Main read · CS231n 2026 lecture 12__

    ---

    Self-supervision, contrastive learning and masked-prediction representations, in the same course as the vision track.

    [:octicons-arrow-right-24: Open schedule](https://cs231n.stanford.edu/schedule)

-   :material-swap-horizontal-bold:{ .lg .middle } __Official lab · PyTorch transfer learning__

    ---

    Compare frozen feature extraction with end-to-end fine-tuning on cost and generalisation (AI-2).

    [:octicons-arrow-right-24: Open tutorial](https://docs.pytorch.org/tutorials/beginner/transfer_learning_tutorial.html)

</div>

## Follow-along · AI-2 · frozen features vs full fine-tuning

1. Run the original transfer-learning example; confirm pretrained weights and augmentation.
2. On the same data split, compare a frozen encoder with full fine-tuning.
3. Reduce training labels or change image conditions; observe accuracy and training cost.
4. Explain why better accuracy alone does not prove causal or complete semantics were learned.

**Done when**: training curves, validation results, time and failure cases; doable after F1, easier with a small GPU.

## Self-check

??? question "With very little labelled data, freeze the encoder or fine-tune fully?"
    Usually freeze first (or fine-tune only the last layers): fewer parameters, less overfitting; unfreeze gradually as data grows.

??? question "What is representation collapse in contrastive learning, and how is it avoided?"
    Every input maps to the same vector; the alignment loss is zero but there is no information. Avoid it with negatives (SimCLR), asymmetric structure with a momentum teacher (BYOL, DINO) or variance regularisation (VICReg).

??? question "Why can MAE mask 75% of image patches?"
    Images are highly redundant; at low mask ratios neighbouring pixels let the model interpolate. A high ratio forces higher-level structure, and the encoder sees only visible patches, so training is faster.

## Related topics

<div class="grid cards" markdown>

-   __Manifolds & latent spaces__

    ---

    The geometric view of representations.

    [:octicons-arrow-right-24: Open](../GENERATIVE-MODELING/MANIFOLDS/manifolds.md)

-   __World models (Foundation Models)__

    ---

    Where JEPA and predictive representations meet.

    [:octicons-arrow-right-24: Open](../../FOUNDATION-MODELS/WORLD-MODELS/world-models.md)

</div>
