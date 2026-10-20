# Foundation Models · Guided Projects

!!! abstract "Reproduce, modify, explain"
    Start with F0, fill gaps through F1, and use the small GPT in F4 as the implementation project; its adaptation part is optional later. Choose F2 or F3 for exploration; F5 is a specialization elective. Run original examples before personal experiments. Use upstream installation instructions matching your selected version.

| Project | Depth | Prerequisite | Resources |
|---|---|---|---|
| F0 · Model exploration | Orientation | Introductory concepts | Browser; local tokenizer optional |
| F1 · Small model training | Fundamentals | Python, arrays and functions | CPU or online notebook |
| F2 · Personal document assistant | Optional application | Python and model calls | Documents; API or local model per tutorial; API costs depend on provider |
| F3 · Image Q&A | Optional exploration | Loading and inference | Model-dependent memory/GPU or online example |
| F4 · Small GPT and adaptation | Core implementation; adaptation elective | F1, tokens and attention | CPU for small experiments; GPU for larger training |
| F5 · Inference comparison | Advanced elective | A working model | Fixed model, test set and hardware |
| F6 · LLM quantization project | Advanced elective | F5, or loading HF models | NVIDIA GPU with 8 GB+ memory |

## F0

### Tokens, context and limits

**Follow:** [HF LLM Course chapter 1](https://huggingface.co/learn/llm-course/en/chapter1/1), then chapter 2 tokenization.

1. Compare tokenization and counts for English and Chinese equivalents.
2. Give a fixed question partial, complete and irrelevant reference material.
3. Include questions unsupported by the material.
4. Record questions, context, answers, model and date.

**Personal change:** vary only context and compare outcomes.

**Acceptance:** a results table and diagram explaining tokenization, context, generation and retrieval.

## F1

### PyTorch image classifier

**Diagnostic:** if you can explain validation/gradient updates and independently save/load a model, proceed to F4. Otherwise study only gaps exposed by Quickstart.

**Follow:** [Learn the Basics](https://docs.pytorch.org/tutorials/beginner/basics/intro.html), using FashionMNIST and the tutorial notebooks.

1. Run Quickstart to see the full workflow.
2. Study tensors, DataLoader, models, autograd and optimization.
3. Save and reload the model in a separate inference script.
4. Change learning rate or width while holding other factors fixed.
5. Split validation data from training data for selection; use the test set for final reporting. Inspect wrong predictions.

**Personal change:** compare two configurations and explain overfitting or undertraining.

**Acceptance:** code, dependencies, training/validation curves, test results and error analysis. F1 also prepares you for robot policies.

## F2

### Personal website document assistant

**Follow:** [LlamaIndex official starter](https://developers.llamaindex.ai/python/framework/getting_started/starter_example/) / [GitHub](https://github.com/run-llama/llama_index). Build the minimum document-Q&A workflow; its default cloud models require API setup and usage costs.

1. Run the original example and record tutorial/dependency versions.
2. Select 5–10 personal Markdown files and inspect parsing and chunks.
3. Build the index; inspect retrieved passages before generated answers.
4. Display source files/passages and handle absent information.
5. Create at least 20 questions: direct facts, cross-passage synthesis and unanswerable questions, with expected evidence.
6. Compare direct answers and RAG. Tune one chunking or retrieval setting, then evaluate on reserved examples.

**Personal change:** use your material, expose sources and build fixed evaluations.

**Acceptance:** runnable app, diagram, question set, support checks and failures; distinguish retrieval and generation failures.

**Next:** improve retrieval/evaluation in the same codebase; see [Applications](../APPLICATIONS/applications.md) for a tool extension rather than start another RAG course.

## F3

### Vision-language capabilities and limits

**Follow:** [HF image-text-to-text](https://huggingface.co/docs/transformers/main/tasks/image_text_to_text). Match development documentation to source, or choose the corresponding stable page.

1. Run one original image Q&A example.
2. Use 10 fixed images for object, counting, spatial and task-related questions.
3. Keep prompts, outputs and human checks, including at least three failure types.
4. Change one prompt or input-size setting and compare quality and inference time.

**Acceptance:** successful and failed examples explaining what remains between image Q&A and reliable actions. Continue to [robot E3–E4](../../EMBODIED-INTELLIGENCE/PROJECT-BASED-LEARNING/project-based-learning.md#e3) for action learning.

## F4

### Small GPT and one adaptation task

**Primary:** [Karpathy Zero to Hero](https://karpathy.ai/zero-to-hero.html) / [companion code](https://github.com/karpathy/nn-zero-to-hero). Micrograd if gradients are unclear, then the first two makemore lectures, followed by the GPT video/code. Use small models/data rather than target a general assistant.

**Text alternative:** [LLMs-from-scratch](https://github.com/rasbt/LLMs-from-scratch), chapters 2–5, can replace this track. The full narrative is in its book; choose one implementation. Later use chapter 6 or 7, or HF's fine-tuning reference, for adaptation.

**Personal change:** vary context length, training data or width. Pretraining and adaptation may use different model sizes; a toy model need not become a general chatbot.

**Acceptance:** architecture, curves, held-out examples, limitations and tensor-shape explanations. Continue to [training](../TRAINING/training.md) for LoRA/SFT practice as needed.

## F5

### Inference on a fixed task

**Follow:** [llama.cpp](https://github.com/ggml-org/llama.cpp) README and run instructions, with the [inference topic](../INFERENCE/inference.md).

Run a supported small model, then compare quantizations of the same base model. Fix hardware, prompts, output length and settings; measure time to first token, generation speed, peak memory and answer quality. Repeat measurements and separate loading from inference.

**Acceptance:** configurations, results, quality regressions and bottleneck analysis. Document limitations when models or hardware differ.


<div class="grid cards" markdown>

-   :material-memory:{ .lg .middle } __F6 · Project: quantizing a small language model__

    ---

    From RTN to GPTQ / AWQ / SmoothQuant: calibrate and quantize yourself, compare accuracy, memory and speed, and analyse per-layer error.

    [:octicons-arrow-right-24: Open project](llm-quant.md)

</div>

## Common record

Keep tutorial links, date, commit, dependencies, model/data revisions, run instructions, results and personal changes. Label tutorial reproductions, personal experiments and independent projects accurately; planned work is not completed work.

[Roadmap](../ROADMAP/roadmap.md)
