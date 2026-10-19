# 具身智能

从机器人系统全貌出发，理解本体、观测、决策与执行的关系，再通过仿真学习经典方法和数据驱动策略。

<div class="grid cards" markdown>

-   :material-rocket-launch:{ .lg .middle } __Get Started · 两周速成__

    ---

    14 天机器人感知闭环：坐标与相机、ROS 2、6DoF 位姿项目接入 ROS 2，再体验一次模仿学习。

    [开始速成](GET-STARTED/get-started.md)

-   :material-compass-outline:{ .lg .middle } __Overview · 全貌导读__

    ---

    机器人、具身智能、模仿学习、强化学习和 VLA 分别是什么？

    [从这里开始](OVERVIEW/overview.md)

-   :material-map-marker-path:{ .lg .middle } __Roadmap · 学习路线__

    ---

    从概念体验到仿真闭环，再选择感知、控制、策略或系统方向。

    [查看路线](ROADMAP/roadmap.md)

-   :material-hammer-wrench:{ .lg .middle } __Guided Projects · 跟随式项目__

    ---

    小海龟、二维机器人算法、仿真模仿学习，逐步增加复杂度。

    [开始实践](PROJECT-BASED-LEARNING/project-based-learning.md)

</div>

## 知识地图

| 单元 | 核心问题 | 学习层次 |
|---|---|---|
| [空间与本体](EMBODIMENT/embodiment.md) | 身体能做什么，坐标、状态和动作怎样描述？ | 基础；复杂动力学选学 |
| [软件系统与仿真](SYSTEMS/systems.md) | 数据、算法和执行器如何连成可运行系统？ | 基础与实践 |
| [感知与状态估计](PERCEPTION/perception.md) | 机器人在哪里，环境中有什么？ | 先导论；按任务深入 |
| [规划与控制](PLANNING-CONTROL/planning-control.md) | 怎样选择路径并修正执行误差？ | 基础；MPC、全身控制选学 |
| [学习型策略](POLICY-LEARNING/policy-learning.md) | 怎样从演示或交互学习动作？ | 先模仿学习；RL 与 VLA 进阶 |
| [仿真与虚实迁移](SIM-TO-REAL/sim-to-real.md) | 仿真经验怎样迁移到硬件，并保持可靠？ | 仿真早体验；真机迁移进阶 |

!!! tip "学习顺序"
    可以先运行简单策略，再按任务补运动学和控制。先掌握坐标、观测动作接口与基本反馈；完整 SLAM、MPC 和全身控制不作为模仿学习入门的统一门槛。

## 与基础模型的交汇

语言模型调用技能、VLA 生成动作、世界模型预测环境是不同接口。它们都需要与机器人的状态、控制和评估系统连接。

[基础模型全貌](../FOUNDATION-MODELS/OVERVIEW/overview.md) · [多模态与行动接口](../FOUNDATION-MODELS/MULTIMODAL/multimodal.md)
