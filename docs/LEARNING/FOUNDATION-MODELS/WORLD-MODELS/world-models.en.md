# World Models

!!! abstract "What this page solves"
    A world model predicts "if I do this, what will the environment become". It can be used to plan in imagination, generate data for policy learning, or capture physical regularities in video. **Advanced elective**; your Manifold Dynamics research (a learned latent state + transition operator) is structurally a world model.

## How far to go

| Level | Learn | Done when |
|---|---|---|
| Orientation | Prediction, state, action conditioning, planning | You can distinguish world models, policies and video generators |
| Foundation | Latent dynamics models, planning in the model | You can follow one latent world model's training and use |
| Advanced | Training policies with world models; prediction vs task benefit | An experiment that separates "predicts well" from "decides well" |

## Core knowledge

**Definition**: given the current state (or observation history) and an action, predict the next state or observation: s_{t+1} = f(s_t, a_t). A video predictor without action input can only watch, not decide.

**Latent modelling**: predicting pixels is expensive and full of irrelevant detail. The mainstream encodes observations into a compact latent state and learns dynamics there (e.g. Dreamer's RSSM), decoding only when needed.

**Planning in the model**: simulate the outcomes of candidate action sequences and execute the best (the MPC idea), or train a policy on imagined trajectories (Dreamer), cutting real interaction drastically.

**Predictive representation learning (JEPA)**: predict the representation of masked or future parts instead of reconstructing pixels, learning more abstract and robust states — usable as a world-model encoder.

**Video generation as a world simulator**: large video generators show some physical common sense and are used to synthesise robot training data or long-horizon predictions; physical consistency and action controllability remain unreliable.

**Compounding error and evaluation**: multi-step prediction errors accumulate. Judge a world model not only by one-step error but by whether policies planned or trained with it do better in the real environment.

## Debate · generating pixels vs. predicting representations

Generation: realistically generating future video demonstrates understanding and yields training data directly. Representation (e.g. LeCun's JEPA): reconstructing every pixel is wasteful and unnecessary; intelligence needs prediction over abstract states. Both are advancing; the final criterion is the benefit to downstream decisions.

## Curated resources

<div class="grid cards" markdown>

-   :material-file-document:{ .lg .middle } __Classic entry · World Models (Ha & Schmidhuber)__

    ---

    A paper page with interactive demos: visual encoder, latent dynamics, training a controller "in the dream".

    [:octicons-arrow-right-24: Open paper page](https://worldmodels.github.io/)

-   :material-file-document:{ .lg .middle } __Mainstream method · DreamerV3__

    ---

    Latent world-model RL that works across many domains with one set of hyper-parameters.

    [:octicons-arrow-right-24: Open project page](https://danijar.com/project/dreamerv3/)

</div>

## Follow-along experiment

<div class="grid cards" markdown>

-   :material-hammer-wrench:{ .lg .middle } __Walk through the World Models page__

    ---

    Follow the demos, draw "encoder → latent dynamics → controller", then write one page comparing it with your Manifold Dynamics work.

    [:octicons-arrow-right-24: Open page](https://worldmodels.github.io/)

</div>

## Self-check

??? question "How does a video predictor differ from a world model?"
    A world model is conditioned on actions and answers "what if I do a"; an unconditioned video predictor only says "what usually happens next" and cannot support decisions.

??? question "Why model dynamics in a latent space rather than pixels?"
    Pixels are high-dimensional and full of task-irrelevant detail (texture, lighting), making prediction hard and slow; a latent space is compact and keeps what matters for dynamics.

??? question "Why can planning be poor even when one-step error is low?"
    Errors accumulate over rollouts, and the policy actively seeks regions where the model is wrong but predicts high return (model exploitation).

??? question "What do Manifold Dynamics and Dreamer share, and how do they differ?"
    Both learn a latent state and a transition operator. Manifold Dynamics detects anomalies as deviation from the learned flow and is not action-conditioned; Dreamer is action-conditioned and trains a policy.

## Related topics

<div class="grid cards" markdown>

-   __Search, Planning & Decision Learning (AI)__

    ---

    Foundations of planning and reinforcement learning.

    [:octicons-arrow-right-24: Open](../../AI/DECISION-LEARNING/decision-learning.md)

-   __Representation Learning (AI)__

    ---

    JEPA and self-supervised representations.

    [:octicons-arrow-right-24: Open](../../AI/REPRESENTATION/representation.md)

</div>
