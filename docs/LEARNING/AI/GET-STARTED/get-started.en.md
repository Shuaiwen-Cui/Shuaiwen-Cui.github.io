# AI · Two-Week Get Started

!!! abstract "How this page differs from the systematic track"
    The systematic track (overview, roadmap, topic pages) is about complete understanding. This page does one thing: **bring deep-learning engineering skills up to what industry algorithm roles commonly ask for, in two weeks, and leave one project you can show**. About 2–3 hours a day. Afterwards, return to the systematic track and go deep where you are interested.

## Target skills (worked back from common role requirements)

| Skill | Usually phrased as | What you already have | Close the gap to |
|---|---|---|---|
| Training & debugging | PyTorch, training loops, overfitting, tuning | Trained Transformer / Mamba models in research projects | Write a training loop from memory and explain every hyper-parameter |
| Mainstream architectures | CNN, Transformer, attention | Principles (see the checklist) | Explain residuals, normalisation, attention complexity |
| Evaluation | Metrics, leakage, validation | Time-ordered splits, AUROC | Detection and segmentation metrics (mAP, IoU) |
| Deployment | ONNX, quantisation, TensorRT, latency | **TinyAuton: kernels and FP8 quantisation on MCUs** | INT8 and latency measurement on GPU / Jetson |

**Your edge:** most candidates only train on GPUs; you can carry the same model all the way down to an MCU. The project below is built around that.

## 14-day plan

| Day | Learn | Do | Output |
|---|---|---|---|
| 1 | The training loop | Rewrite the PyTorch Quickstart training/validation without looking | A running script with your own comments |
| 2 | Splits, overfitting, metrics | Cause overfitting and data leakage on purpose, once each | Two comparison curves |
| 3 | CNNs and transfer learning | Fine-tune ResNet18 on a small dataset (frozen vs full) | Accuracy and training-time table |
| 4 | Transformer mechanics | Read the site's Transformer page; first half of Karpathy "Let's build GPT" | Tensor shapes of one attention block, drawn by hand |
| 5 | Training tricks | Try warm-up, normalisation and mixed precision once each | Which one mattered most |
| 6–7 | Detection / segmentation | Fine-tune YOLO11n on a small dataset with Ultralytics | mAP and PR curve; explain IoU and mAP |
| 8–9 | Quantisation | ONNX Runtime quantisation: export ONNX, static INT8 | FP32 vs INT8 accuracy and CPU latency |
| 10–11 | GPU / Jetson deployment | Ultralytics Jetson guide: TensorRT FP16 and INT8 | Latency table on desktop GPU and Jetson |
| 12 | Connect your strength | Run a small network for the same task on ESP32-S3 with TinyAuton | Accuracy–latency–memory across GPU / Jetson / MCU |
| 13 | Package the project | README: problem, method, results, limits | Repository front page |
| 14 | Interview review | Answer every question below; revisit the day behind any gap | Your own answer notes |

## Resources for these two weeks

<div class="grid cards" markdown>

-   :material-open-in-new:{ .lg .middle } __PyTorch Quickstart__

    ---

    Day 1 · The training loop

    [:octicons-arrow-right-24: Open](https://docs.pytorch.org/tutorials/beginner/basics/quickstart_tutorial.html)

-   :material-book-open-variant:{ .lg .middle } __Transformer (this site)__

    ---

    Day 4 · Transformer mechanics

    [:octicons-arrow-right-24: Open](../TRANSFORMER/transformer.md)

-   :material-open-in-new:{ .lg .middle } __Karpathy · Zero to Hero__

    ---

    Day 4 · Transformer mechanics

    [:octicons-arrow-right-24: Open](https://karpathy.ai/zero-to-hero.html)

-   :material-open-in-new:{ .lg .middle } __ONNX Runtime quantisation__

    ---

    Day 8–9 · Quantisation

    [:octicons-arrow-right-24: Open](https://onnxruntime.ai/docs/performance/model-optimizations/quantization.html)

-   :material-open-in-new:{ .lg .middle } __Ultralytics Jetson guide__

    ---

    Day 10–11 · GPU / Jetson deployment

    [:octicons-arrow-right-24: Open](https://docs.ultralytics.com/guides/nvidia-jetson/)

</div>


## Project: one model, from GPU to MCU

One table that answers: for the same task, what accuracy, speed and memory does the model reach on a GPU, a Jetson and an MCU, how much does quantisation cost at each level, and why. Report measured results only; mark unfinished parts as planned.

## Must-know interview questions (answer key)

| Question | Key points |
|---|---|
| Why residual connections? | Gradients bypass the non-linear layers, so deep networks optimise more easily; each layer learns an increment |
| BatchNorm vs LayerNorm? | BN normalises across the batch, depends on batch size and uses running statistics at inference; LN normalises across features, suits sequences and small batches; Transformers use LN |
| Complexity of self-attention, and how to reduce it? | O(n²) time and memory in sequence length; sparse / linear attention, FlashAttention (same complexity, far less memory traffic), state-space models (Mamba, linear) |
| Adam or SGD? | Adam converges fast and tolerates learning-rate choice; SGD + momentum often generalises better in vision; large models use AdamW (decoupled weight decay) |
| What is data leakage? | Test information reaching training or model selection, e.g. random splits of time series or normalising on all data; metrics become inflated |
| How is mAP computed? | Sort by confidence, build the PR curve, integrate to AP, average over classes (and IoU thresholds) |
| PTQ vs QAT? | PTQ sets scales from a small calibration set after training: fast, may lose accuracy. QAT simulates quantisation error during training: better accuracy, needs retraining |
| INT8 lost accuracy — how to debug? | Compare per-layer outputs before and after quantisation to find sensitive layers (often first/last and outlier-heavy ones); per-channel scales, different calibration, or keep a few layers in FP16 |
| Why are fewer FLOPs not always faster? | Many layers are memory-bound, not compute-bound; operator fusion and support in the runtime decide real latency |
| First things to check for edge deployment? | Memory (weights + peak activations), latency target, operator support, power; then choose model size and quantisation |

## Afterwards

Return to round three of the AI roadmap and pick one direction; for deployment continue with Efficient Training & Engineering.


<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Overview (systematic track)__

    ---

    Systematic learning

    [:octicons-arrow-right-24: Open](../ai.md)

-   :material-book-open-variant:{ .lg .middle } __AI roadmap__

    ---

    Systematic learning

    [:octicons-arrow-right-24: Open](../ROADMAP/roadmap.md)

-   :material-book-open-variant:{ .lg .middle } __Efficient Training & Engineering__

    ---

    Systematic learning

    [:octicons-arrow-right-24: Open](../SYSTEMS/systems.md)

</div>
