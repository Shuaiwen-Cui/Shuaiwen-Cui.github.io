# Training & Adaptation

!!! abstract "What this page solves"
    The stages that take a large model from "read a lot" to "follows instructions and preferences", and what you can do with limited compute (usually adaptation, not pretraining from scratch). **Concepts required, practice up to LoRA**; prerequisites: F1's training loop and data splits.

## How far to go

| Level | Learn | Done when |
|---|---|---|
| Orientation | Pretraining, SFT, preference optimisation, parameter-efficient tuning | You can say what each stage changes in the model |
| Foundation | Splits, hyper-parameters, checkpoints, LoRA | One small-model LoRA fine-tune evaluated on held-out examples |
| Advanced | Scaling laws, parallel training, preference optimisation, RL with verifiable rewards | One question explored in depth |

## Core knowledge

**Pretraining objective**: next-token prediction (cross-entropy) on vast text. Most knowledge and capability — and most cost — come from here.

**Scaling laws**: loss falls as a power law in parameters, data and compute. Chinchilla: for a fixed compute budget scale parameters and data together (~20 tokens per parameter); small models today are often "over-trained" to make inference cheaper.

**Data**: deduplication, quality filtering and mixture (web, code, maths, languages) often matter more than architecture changes. Benchmark contamination inflates metrics.

**Distributed training**: data parallel (a copy per GPU, split data), tensor parallel (split matrices), pipeline parallel (split layers); ZeRO / FSDP shard optimiser states, gradients and parameters across GPUs to save memory.

**Instruction tuning (SFT)**: train on instruction–response pairs so the model learns the chat format and to follow instructions; loss only on the response.

**Preference optimisation**: RLHF trains a reward model then optimises with PPO; DPO trains directly on "better vs worse" pairs, with no separate reward model or RL loop — simpler and more stable.

**RL with verifiable rewards**: for maths and code, where correctness can be checked automatically, use it as the reward (e.g. GRPO); the main source of recent reasoning gains.

**Parameter-efficient tuning**: LoRA freezes W and trains a low-rank update ΔW = BA, usually under 1% of parameters; QLoRA quantises the base to 4-bit and trains LoRA on top, so billion-parameter models fit on one GPU.

**Distributed training systems**: data parallelism (DDP) keeps a full model per GPU and synchronises gradients with NCCL all-reduce; when the model does not fit, ZeRO / FSDP shard optimiser states, gradients and parameters across GPUs (DeepSpeed, PyTorch FSDP); larger models add tensor parallelism (splitting matrices within a layer) and pipeline parallelism (splitting layers into stages), as in Megatron-LM. Activation recomputation (gradient checkpointing) trades time for memory, and mixed precision (BF16) is close to the default.

## Debate · keep scaling vs. new methods

Scaling: more data and compute remain the most reliable source of progress. Methods: better data, synthetic data, RL and test-time compute (thinking longer at inference) deliver similar gains more cheaply. The recent trend combines both: pretraining scale slows while post-training and inference-time compute grow.

## Curated resources

<div class="grid cards" markdown>

-   :material-code-braces:{ .lg .middle } __Foundation · Karpathy__

    ---

    Loss, train / validation and sampling on one small model (F4).

    [:octicons-arrow-right-24: Open course](https://karpathy.ai/zero-to-hero.html)

-   :hugging:{ .lg .middle } __Adaptation · Hugging Face PEFT__

    ---

    LoRA / QLoRA principles and usage; one small adaptation with independent evaluation.

    [:octicons-arrow-right-24: Open docs](https://huggingface.co/docs/peft/index)

-   :material-school:{ .lg .middle } __Deeper (elective) · Stanford CS336__

    ---

    Language models from scratch: data, architecture, training systems; only if model training is your direction.

    [:octicons-arrow-right-24: Open course](https://cs336.stanford.edu/)

</div>

## Follow-along experiment

<div class="grid cards" markdown>

-   :material-hammer-wrench:{ .lg .middle } __F4 · A small GPT and one adaptation__

    ---

    Pretrain a small model, then LoRA-adapt an existing small model and evaluate.

    [:octicons-arrow-right-24: See steps](../PROJECT-BASED-LEARNING/project-based-learning.md)

</div>

## Self-check

??? question "Why compute SFT loss only on the response?"
    The instruction is given input the model need not learn to generate; including it wastes capacity and encourages parroting prompts.

??? question "What if LoRA rank r is too small or too large?"
    Too small: not expressive enough. Too large: more parameters, more overfitting, diminishing returns. 8–64 is common, with alpha scaling.

??? question "What does DPO remove compared with RLHF?"
    The separately trained reward model and the PPO loop; it applies a classification-style loss directly to preference pairs — more stable and cheaper.

??? question "Why are small models trained far beyond the Chinchilla-optimal token count?"
    Chinchilla optimises training compute only; inference cost accumulates with every call, so a smaller model trained on more data is cheaper overall.

## Related topics

<div class="grid cards" markdown>

-   __Inference & Deployment__

    ---

    Making the trained model fast and lean.

    [:octicons-arrow-right-24: Open](../INFERENCE/inference.md)

-   __Evaluation & Trustworthy AI (AI)__

    ---

    Showing that an adaptation really works.

    [:octicons-arrow-right-24: Open](../../AI/EVALUATION/evaluation.md)

</div>
