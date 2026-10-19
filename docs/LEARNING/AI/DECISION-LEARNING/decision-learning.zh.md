# 搜索、规划与决策学习

!!! abstract "这一页解决什么"
    AI 不只做预测，还要在状态、约束和反馈下选择行动。这一页区分搜索规划、模仿和强化学习，是机器人策略、世界模型和智能体的共同基础。**P1 认识**；强化学习作为一个可选的深入方向。

## 知识地图

| 层次 | 主题 | 需要理解什么 |
|---|---|---|
| P1 | 经典决策与符号方法 | 搜索、启发式、约束、知识表示 |
| P1 | 强化学习 | MDP、价值、策略、奖励、探索；DQN、策略梯度、PPO |
| P1 | 模仿与离线学习 | 行为克隆、分布偏移；离线 RL 与监督学习的区别 |
| P2 | 模型、推理与交互 | 模型式 RL、世界模型、分层技能、多智能体 |

## 核心机制

**马尔可夫决策过程（MDP）**：状态 s、动作 a、转移 P(s'|s,a)、奖励 r、折扣 γ。目标是找一个策略 π(a|s) 最大化期望累计折扣奖励。

**价值函数**：V(s) 是从状态 s 出发的期望回报，Q(s,a) 是先做动作 a 再按策略走的期望回报。贝尔曼方程把它们写成"当前奖励 + 折扣后的下一状态价值"，是几乎所有 RL 算法的核心。

**基于价值的方法**：Q-learning 用贝尔曼方程迭代更新 Q；DQN 用神经网络近似 Q，加经验回放（打破样本相关性）和目标网络（稳定训练目标）。适合离散动作。

**策略梯度与 Actor-Critic**：直接优化策略参数，能处理连续动作；用 Critic 估计的价值作基线降低方差。PPO 限制每次策略更新的幅度，稳定好用，是机器人 RL 和 RLHF 的默认选择。

**探索与利用**：只利用已知最好的动作会错过更好的选择；ε-贪婪、熵正则、好奇心奖励是常见探索手段。奖励稀疏时探索是最大难点。

**模仿 vs. 强化**：模仿学习从演示学，不需要奖励但受演示覆盖限制；强化从奖励学，能超越示范但样本效率低、奖励难设计。离线 RL 只用固定数据集，不能与环境交互，难点是对数据外动作的价值过度乐观。

**模型式 RL**：先学环境模型（世界模型），在模型中规划或生成数据，样本效率高，但模型误差会被策略利用。

## 精选资源

<div class="grid cards" markdown>

-   :material-compass-outline:{ .lg .middle } __概念主课 · David Silver RL__

    ---

    经典强化学习课程，重点是 MDP、价值函数、策略迭代与时序差分。

    [:octicons-arrow-right-24: 打开课程主页](https://www.davidsilver.uk/teaching/)

-   :material-gamepad-variant:{ .lg .middle } __官方实验 · PyTorch DQN CartPole__

    ---

    经验回放、ε-贪婪探索与评估的完整小例子（AI-4）。

    [:octicons-arrow-right-24: 打开教程](https://docs.pytorch.org/tutorials/intermediate/reinforcement_q_learning.html)

</div>

## 跟随实验 · AI-4 · CartPole 与重复评估

1. 跑官方原例，标出环境、观测、动作、奖励与终止。
2. 区分训练时的探索策略与评估策略，保存 checkpoint。
3. 至少三个训练种子；固定评估条件，报告多次 episode 的回报分布。
4. 只改变一个超参或扰动，记录训练与推理成本及失败。

**验收**：CPU 可开始；一个成功视频不能作为稳定性证据。机器人方向先做模仿学习 E3，RL 作为另一分支。

## 自测

??? question "DQN 为什么需要经验回放和目标网络？"
    连续采到的样本高度相关，直接训练不稳定，经验回放打乱样本；目标网络让 TD 目标在一段时间内固定，避免"追着自己的尾巴跑"。

??? question "为什么 RL 结果要报告多个随机种子？"
    RL 对初始化和探索过程很敏感，同一算法不同种子的结果差异可能比不同算法之间还大；单个种子的结果不可信。

??? question "离线 RL 为什么不能直接用普通 Q-learning？"
    对数据集中没出现过的动作，Q 网络的估计可能被高估，策略会去选这些动作，而又没有机会通过交互纠正。

## 相关专题

<div class="grid cards" markdown>

-   __学习型策略（具身）__

    ---

    模仿与强化学习在机器人上的使用。

    [:octicons-arrow-right-24: 进入](../../EMBODIED-INTELLIGENCE/POLICY-LEARNING/policy-learning.md)

-   __世界模型（基础模型）__

    ---

    模型式 RL 与在想象中规划。

    [:octicons-arrow-right-24: 进入](../../FOUNDATION-MODELS/WORLD-MODELS/world-models.md)

-   __应用与智能体（基础模型）__

    ---

    语言模型作为决策者调用工具。

    [:octicons-arrow-right-24: 进入](../../FOUNDATION-MODELS/APPLICATIONS/applications.md)

</div>
