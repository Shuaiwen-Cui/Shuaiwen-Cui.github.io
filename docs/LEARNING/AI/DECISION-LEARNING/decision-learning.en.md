# Search, Planning & Decision Learning

!!! abstract "What this page solves"
    AI does not only predict; it chooses actions under states, constraints and feedback. This page separates search and planning, imitation and reinforcement learning — the common base of robot policies, world models and agents. **P1: know it**; RL is an optional deep dive.

## Knowledge map

| Level | Topic | What to understand |
|---|---|---|
| P1 | Classical decision and symbolic methods | Search, heuristics, constraints, knowledge representation |
| P1 | Reinforcement learning | MDP, value, policy, reward, exploration; DQN, policy gradient, PPO |
| P1 | Imitation and offline learning | Behaviour cloning, distribution shift; offline RL vs supervised learning |
| P2 | Models, reasoning and interaction | Model-based RL, world models, hierarchical skills, multi-agent |

## Core mechanics

**Markov decision process (MDP)**: state s, action a, transition P(s'|s,a), reward r, discount γ. The goal is a policy π(a|s) maximising expected cumulative discounted reward.

**Value functions**: V(s) is the expected return from state s; Q(s,a) is the expected return after taking action a and then following the policy. The Bellman equation writes them as "current reward + discounted value of the next state" — the core of almost every RL algorithm.

**Value-based methods**: Q-learning updates Q iteratively with the Bellman equation; DQN approximates Q with a network plus experience replay (breaks sample correlation) and a target network (stabilises the target). Suits discrete actions.

**Policy gradient and actor-critic**: optimise policy parameters directly, handling continuous actions; a critic's value estimate serves as a baseline to cut variance. PPO limits each policy update, is stable and practical, and is the default for robot RL and RLHF.

**Exploration vs exploitation**: always taking the best-known action misses better ones; ε-greedy, entropy regularisation and curiosity rewards are common exploration tools. With sparse rewards, exploration is the hardest part.

**Imitation vs reinforcement**: imitation learns from demonstrations, needs no reward but is limited by demo coverage; RL learns from reward, can exceed the demos, but is sample-inefficient and rewards are hard to design. Offline RL uses only a fixed dataset with no interaction; the difficulty is over-optimistic values for out-of-data actions.

**Model-based RL**: learn an environment model (world model) first and plan or generate data in it; sample-efficient, but the policy can exploit model errors.

## Curated resources

<div class="grid cards" markdown>

-   :material-compass-outline:{ .lg .middle } __Concept course · David Silver RL__

    ---

    The classic RL course: MDPs, value functions, policy iteration and temporal difference.

    [:octicons-arrow-right-24: Open course page](https://www.davidsilver.uk/teaching/)

-   :material-gamepad-variant:{ .lg .middle } __Official lab · PyTorch DQN CartPole__

    ---

    A complete small example of experience replay, ε-greedy exploration and evaluation (AI-4).

    [:octicons-arrow-right-24: Open tutorial](https://docs.pytorch.org/tutorials/intermediate/reinforcement_q_learning.html)

</div>

## Follow-along · AI-4 · CartPole with repeated evaluation

1. Run the official example; mark environment, observation, action, reward and termination.
2. Separate the exploration policy used in training from the evaluation policy; save checkpoints.
3. At least three training seeds; fix evaluation conditions and report the return distribution over many episodes.
4. Change only one hyperparameter or perturbation; record training and inference cost and failures.

**Done when**: starts on a CPU; one success video is not evidence of stability. For robotics, do imitation learning E3 first; RL is a separate branch.

## Self-check

??? question "Why does DQN need experience replay and a target network?"
    Consecutive samples are highly correlated, making direct training unstable; replay shuffles them. The target network fixes the TD target for a while so the network is not chasing its own tail.

??? question "Why report RL results over several random seeds?"
    RL is very sensitive to initialisation and exploration; one algorithm's spread across seeds can exceed the gap between algorithms, so a single seed is not trustworthy.

??? question "Why can't offline RL just use ordinary Q-learning?"
    For actions absent from the dataset, Q estimates may be overestimated; the policy picks those actions and never gets to correct them through interaction.

## Related topics

<div class="grid cards" markdown>

-   __Learned policies (Embodied)__

    ---

    Imitation and RL on robots.

    [:octicons-arrow-right-24: Open](../../EMBODIED-INTELLIGENCE/POLICY-LEARNING/policy-learning.md)

-   __World models (Foundation Models)__

    ---

    Model-based RL and planning in imagination.

    [:octicons-arrow-right-24: Open](../../FOUNDATION-MODELS/WORLD-MODELS/world-models.md)

-   __Applications & agents (Foundation Models)__

    ---

    Language models as decision-makers calling tools.

    [:octicons-arrow-right-24: Open](../../FOUNDATION-MODELS/APPLICATIONS/applications.md)

</div>
