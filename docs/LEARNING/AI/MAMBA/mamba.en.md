# Mamba & State-Space Sequence Models

!!! abstract "What this page solves"
    Understand another family of sequence models through "how the state updates, how history is compressed, how compute scales with length". Your familiarity with state-space equations and system identification is a strong start, and you have already swapped Transformer for Mamba in a tunnelling project. **P1: know it; P2: go deep.**

## Core mechanics

**Linear state-space model**: continuous form h'(t) = A h(t) + B x(t), y(t) = C h(t). After discretisation (e.g. zero-order hold) it becomes the recurrence h_t = Ā h_{t−1} + B̄ x_t — the same state equation as in control theory.

**Two ways to compute**: when A, B, C do not depend on the input (LTI, e.g. S4), the recurrence unrolls into one long convolution for parallel training; inference uses the recurrence, O(1) per step with no cache growing with length.

**Selectivity (Mamba's key idea)**: make B, C and the step size Δ depend on the current input, so the model decides what to remember and what to ignore — fixing LTI models' inability to select by content. The cost: no convolution form, so a hardware-aware parallel scan is used instead.

**Trade-off against attention**: attention keeps all history and can look back exactly at any position, but its KV cache grows with length; an SSM compresses history into a fixed-size state, saving memory and speeding long-sequence inference, but is weaker at exact recall of distant detail. Hence hybrids (a few attention layers + many SSM layers) are increasingly common.

**Mamba-2 and SSD**: shows a duality between structured SSMs and a class of linear attention, with better training efficiency.

**Not a physical state**: the SSM "state" is a learned abstract vector, not necessarily a real physical state; do not read it as an interpretable system-identification result.

| Level | Topic | Distinguish |
|---|---|---|
| P1 | Linear state space, discretisation, recurrence | Mathematical state ≠ physical state |
| P1 | Selective SSM / Mamba | Input-dependent selection, not a fixed linear filter |
| P2 | Mamba-2, SSD, hybrid sequence models | Structural links, training and inference cost, accuracy |
| P2 | Time-series and long-sequence use | Data, sampling rate, length and task set the comparison |

## Curated resources

<div class="grid cards" markdown>

-   :material-github:{ .lg .middle } __Official entry · Mamba repository__

    ---

    The authors' papers, architecture diagrams and CUDA implementation of the hardware-aware scan.

    [:octicons-arrow-right-24: Open repository](https://github.com/state-spaces/mamba)

-   :material-video:{ .lg .middle } __Intuition · Mamba explained__

    ---

    Chinese video contrasting continuous state space with attention.

    [:octicons-arrow-right-24: Watch video](https://www.bilibili.com/video/BV1QrBbBhEyN/)

</div>

## Small experiment: build a comparable baseline first

On data from the time-series & scientific AI topic (see related topics below), first compare persistence / linear prediction with a small network; then add Mamba, fixing time split, input window, training budget and device, and record prediction quality, latency and memory.

**Done when**: you can explain state compression vs attention's access to history, without claiming Mamba is faster in practice or has replaced Transformers from asymptotic complexity alone.

## Self-check

??? question "Why can an LTI SSM be written as a convolution but Mamba cannot?"
    In an LTI model Ā, B̄, C are fixed, so the output is the input convolved with one fixed kernel; Mamba's parameters change with the input, so each step's "kernel" differs and only a scan works.

??? question "In long-sequence inference, where does Mamba save over a Transformer?"
    A Transformer's KV cache grows linearly with length and each new token reads the whole cache; Mamba keeps a fixed-size state, so per-step compute and memory are constant.

??? question "On what tasks might Mamba lose to attention?"
    Tasks that require exact copying or retrieval of specific distant content (e.g. finding a phone number in a long document): a fixed-size state struggles to keep every detail losslessly.

## Related topics

<div class="grid cards" markdown>

-   __Architectures__

    ---

    CNN, Transformer, SSM and GNN compared as computation patterns.

    [:octicons-arrow-right-24: Open](../ARCHITECTURES/architectures.md)

-   __Time-series & scientific AI__

    ---

    Run the comparison on signal data you know.

    [:octicons-arrow-right-24: Open](../SCIENTIFIC-AI/scientific-ai.md)

</div>
