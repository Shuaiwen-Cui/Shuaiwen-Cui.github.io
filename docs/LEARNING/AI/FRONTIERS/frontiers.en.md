# Frontier Index · Place Techniques in the Map

!!! note "Scope and evidence"
    Representative families checked on 2026-10-18, not a popularity ranking or exhaustive research inventory. Use courses, author code and documentation for scope; prioritize objectives, interfaces and limits over model names.

## Architectures, representations and generation

| Family | Role | Entry | First-pass scope |
|---|---|---|---|
| Transformer, ViT, DiT | Similar components across language, vision and diffusion; architecture does not fix the objective | [Architectures](../ARCHITECTURES/architectures.md) / [generation](../GENERATIVE-MODELING/generative-modeling.md) | Attention and data flow |
| Efficient attention and long context | Structure, approximations, positions and implementations affect cost/context | [Transformer](../TRANSFORMER/transformer.md) / [systems](../SYSTEMS/systems.md) | Linear attention vs. FlashAttention |
| MoE | Partial expert activation with routing/balance/communication trade-offs | [Architectures](../ARCHITECTURES/architectures.md) | Total vs. active parameters |
| SSM, Mamba and hybrids | State updates and selection for sequences | [Mamba](../MAMBA/mamba.md) | Compare recurrence and attention conceptually |
| Self-supervision, contrast, CLIP, DINO, MAE | Different objectives from views, teachers and masks | [Representation](../REPRESENTATION/representation.md) | Frozen vs. fine-tuned features |
| JEPA and predictive representations | Predict features rather than pixels | [Representation](../REPRESENTATION/representation.md) / [world models](../../FOUNDATION-MODELS/WORLD-MODELS/world-models.md) | Read [V-JEPA 2 author code](https://github.com/facebookresearch/vjepa2), not retraining |
| VAE, GAN and normalizing flows | Latent-variable, adversarial and invertible-transform approaches | [Generation](../GENERATIVE-MODELING/generative-modeling.md) | Objectives and sampling |
| Diffusion, flow matching, rectified flow | Generation along noise/probability paths with different targets | [Diffusion](../GENERATIVE-MODELING/DIFFUSION/diffusion.md) / [flow matching](../GENERATIVE-MODELING/FLOW-MATCHING/flow-matching.md) | Same low-dimensional labs |
| Consistency/distillation and few-step sampling | Extra training/objective design to reduce sampling steps | [Generation](../GENERATIVE-MODELING/generative-modeling.md) / [systems](../SYSTEMS/systems.md) | Purpose; acceleration has trade-offs |
| Discrete diffusion/non-autoregressive generation | Discrete-state generation beyond token-wise autoregression | [Generation](../GENERATIVE-MODELING/generative-modeling.md) | Optional MIT lecture 5 |

## Perception, reasoning, decisions and systems

| Family | Role | Entry | First-pass scope |
|---|---|---|---|
| Prompted segmentation, VLM and video | Conditional/multimodal processing beyond fixed-label tasks | [Vision](../CV/cv.md) / [multimodality](../../FOUNDATION-MODELS/MULTIMODAL/multimodal.md) | Tasks, metrics and spatial limits |
| NeRF, 3DGS, implicit geometry and points | Different shape/appearance representations | [3D vision](../CV/cv.md) | Appearance vs. geometry |
| Graph Transformers and equivariance | Explicit relations and symmetry | [Graphs](../GNN/gnn.md) | Small GCN and edge-free baseline |
| SFT, LoRA/PEFT and preferences | Task and behavior adaptation | [Training](../../FOUNDATION-MODELS/TRAINING/training.md) | One task with held-out checks |
| RLVR, reasoning budgets and test-time compute | Verifiable feedback or added computation for specific tasks | [Training](../../FOUNDATION-MODELS/TRAINING/training.md) / [decisions](../DECISION-LEARNING/decision-learning.md) | Rewards, verification and costs; GRPO details later |
| RAG, tools, memory, agents and multi-agent systems | Compose models with material/interfaces | [Applications](../../FOUNDATION-MODELS/APPLICATIONS/applications.md) | One traceable tool workflow |
| Imitation, offline/model-based RL and VLA | Learn actions from demonstrations/interaction | [Decisions](../DECISION-LEARNING/decision-learning.md) / [robot policies](../../EMBODIED-INTELLIGENCE/POLICY-LEARNING/policy-learning.md) | ACT and closed-loop evaluation |
| World models and latent planning | Predict environment changes and evaluate decision benefit | [World models](../../FOUNDATION-MODELS/WORLD-MODELS/world-models.md) | Prediction vs. task success |
| Time-series foundation models, PINNs and operators | Temporal transfer, physics constraints and function maps | [Scientific AI](../SCIENTIFIC-AI/scientific-ai.md) | Baselines and controlled experiments |
| Quantization, distillation, sparsity and runtimes | Algorithmic vs. hardware efficiency | [Systems](../SYSTEMS/systems.md) | Quality, memory and latency on one device |
| Continual/active/federated learning | Changing tasks/data, annotation budgets or distributed data constraints | [Paradigms](../REPRESENTATION/representation.md) | Understand problems before adding training stacks |
| Calibration, OOD, explanations and contamination | Confidence, shifts, evidence and failure | [Evaluation](../EVALUATION/evaluation.md) | Checks shared by all projects |

## A few coverage references

- [CS231n 2026](https://cs231n.stanford.edu/schedule): vision, representations, generation, 3D and vision-language in one curriculum.
- [CS336 2026](https://cs336.stanford.edu/): data, training systems, alignment and reasoning RL; substantial Python/deep-learning/systems prerequisites make it a later specialization.
- [MIT 6.S184 2026](https://diffusion.csail.mit.edu/2026/index.html): continuous/discrete generation, sampling, latent spaces and guided labs.

These inform scope and classification, not new required courses. Place new models by task, architecture, objective, data, cost and evidence; add pages only for distinct learning problems.

[AI overview](../ai.md) · [Roadmap](../ROADMAP/roadmap.md)
