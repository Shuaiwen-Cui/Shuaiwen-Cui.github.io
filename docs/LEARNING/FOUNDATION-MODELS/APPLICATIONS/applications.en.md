# Applications & Agents

!!! abstract "What this page solves"
    Turning a model into a system that solves real problems: connect your own documents (retrieval augmentation), connect tools (agents), and prove it works (evaluation). **Required practice**: much day-to-day LLM work in industry lives here.

## How far to go

| Level | Learn | Done when |
|---|---|---|
| Orientation | Prompts, context, RAG, tool calls, agents | One QA flow with retrieval and a tool drawn |
| Foundation | Chunking, embeddings, retrieval, source display, test sets | F2 done; retrieval and generation failures told apart |
| Advanced | Agent loops, structured output, multi-step planning, automated evaluation | One read-only tool added to F2 and evaluated |

## Core knowledge

**Prompt and context engineering**: the model only sees its context. System prompts set role and rules, few-shot examples set format, relevant material sets the source of facts. Longer context costs more and dilutes attention.

**Retrieval-augmented generation (RAG)**: chunk documents, embed and index them; at question time retrieve the most relevant chunks into the context and answer from them. Knowledge stays updatable and answers traceable.

**Chunking and retrieval quality**: chunks too large add noise, too small lose context; hybrid retrieval (BM25 keywords + vectors) and rerankers improve recall markedly. Most RAG problems are retrieval problems, not generation problems.

**Tool calling**: the model emits a structured call (name + arguments), the program executes it and returns the result. The model executes nothing itself; permissions and safety live in the program.

**Agents**: a "think → call tool → observe → continue" loop until the task is done. Powerful but brittle: more steps, more accumulated error; limit steps and validate each output.

**Structured output**: ask for JSON matching a schema, or use constrained decoding to guarantee the format, so downstream code parses reliably.

**Evaluation**: a fixed question set (including unanswerable and adversarial ones), human or model grading, statistics by failure type. Without a test set every "improvement" is a feeling.

## Debate · long context vs. retrieval

Long context: windows reach a million tokens, so put everything in — simpler, nothing missed. Retrieval: cheaper, faster, traceable, and information in the middle of long contexts gets overlooked. In practice both: retrieve to narrow down, then use a long-enough context for the relevant material.

## Curated resources

<div class="grid cards" markdown>

-   :material-database-search:{ .lg .middle } __Practice · LlamaIndex Starter__

    ---

    The smallest document QA flow: load, chunk, index, retrieve, generate.

    [:octicons-arrow-right-24: Open tutorial](https://developers.llamaindex.ai/python/framework/getting_started/starter_example/)

-   :hugging:{ .lg .middle } __Agent concepts · HF Agents Course__

    ---

    Tool calling, agent loops and evaluation, from concepts to small examples.

    [:octicons-arrow-right-24: Open course](https://huggingface.co/learn/agents-course/unit0/introduction)

</div>

## Follow-along experiment

<div class="grid cards" markdown>

-   :material-hammer-wrench:{ .lg .middle } __F2 · QA over your website__

    ---

    Index your own Markdown, show sources, evaluate on 20 questions, compare plain QA with RAG.

    [:octicons-arrow-right-24: See steps](../PROJECT-BASED-LEARNING/project-based-learning.md)

</div>

## Self-check

??? question "A RAG answer is wrong. How do you tell retrieval from generation?"
    Check whether the retrieved chunks contain the answer: if not, it is retrieval (chunking, embeddings, top-k); if yes but wrong, it is generation (prompt, model, context too long).

??? question "Why include questions the documents cannot answer?"
    To test whether the model admits it does not know instead of inventing an answer — the most common hallucination scenario.

??? question "Why do agents get less reliable with more steps?"
    Each step has some error rate and errors propagate; ten steps at 95% each gives about 60% overall.

??? question "Where is the safety boundary of tool calling?"
    The model only requests calls; the program executes them, validates arguments, enforces permissions and requires human confirmation for side effects.

## Related topics

<div class="grid cards" markdown>

-   __Inference & Deployment__

    ---

    Running applications offline and with low latency.

    [:octicons-arrow-right-24: Open](../INFERENCE/inference.md)

-   __Evaluation & Trustworthy AI (AI)__

    ---

    Evaluation methods, uncertainty and reliability.

    [:octicons-arrow-right-24: Open](../../AI/EVALUATION/evaluation.md)

</div>
