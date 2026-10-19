# 软件系统与仿真

!!! abstract "这一页解决什么"
    把观测、算法和执行连成一个能运行、能调试、能复现的系统。工业界的机器人岗位几乎都要求 ROS 2；仿真决定你能不能在没有真机时迭代。**必修**；前置：Python、Linux 基本操作。

## 学到什么程度

| 层次 | 学什么 | 完成标准 |
|---|---|---|
| 扫盲 | 节点、消息、坐标、仿真器的作用 | 画出一个机器人任务的节点图 |
| 基础 | ROS 2 发布订阅、tf2、rosbag、RViz | 完成 E1；自己的坐标系在 RViz 中正确显示 |
| 进阶 | 仿真器接口、时延与实时性、系统集成 | 能从日志定位"没收到数据、坐标错、执行接口不匹配"三类故障 |

## 核心知识

**节点与通信方式**：ROS 2 把系统拆成节点。话题（topic）是持续的数据流，如相机图像；服务（service）是一问一答，如"重置"；动作（action）是带进度反馈的长任务，如"移动到某点"。选错通信方式是常见的设计问题。

**QoS**：规定可靠性（reliable / best-effort）、历史深度和持久性。传感器数据通常用 best-effort，丢一帧比堵塞好；控制命令和地图用 reliable。发布方和订阅方 QoS 不兼容时会静默收不到数据。

**tf2 与坐标树**：维护所有坐标系之间随时间变化的变换。查询"某时刻相机坐标系下的点在地图坐标系中在哪"时，tf2 会沿树相乘并按时间插值。时间戳不对是融合出错的头号原因。

**URDF**：用 XML 描述连杆、关节、惯量和外观，是 RViz 显示、运动学求解和仿真的共同输入。

**rosbag 记录与回放**：把传感器数据录下来，算法可以在完全相同的输入上反复调试，把"传感问题"和"算法问题"分开。

**仿真器**：MuJoCo（接触动力学准确，适合操作和学习）、Isaac Sim / Isaac Lab（GPU 并行、渲染真实，适合大规模训练和合成数据）、Gazebo（与 ROS 集成紧密，适合移动机器人）。选一个就够，按任务选。

## 路线之争 · ROS 2 vs. 自研中间件

ROS 2 生态大、招聘通用，但有学习成本和性能开销；很多量产公司（自动驾驶、消费机器人）自研轻量中间件以控制时延和依赖。学 ROS 2 的概念（节点、消息、tf、记录回放）在两种环境里都通用。

## 精选资源

<div class="grid cards" markdown>

-   :material-robot-industrial:{ .lg .middle } __主线 · ROS 2 官方教程__

    ---

    Beginner：CLI 工具 → Python 发布订阅；Intermediate：tf2。按你安装的发行版选文档。

    [:octicons-arrow-right-24: 打开教程](https://docs.ros.org/en/jazzy/Tutorials.html)

-   :material-axis-arrow:{ .lg .middle } __坐标 · tf2 教程__

    ---

    广播与监听坐标变换、时间旅行查询；感知与融合离不开它。

    [:octicons-arrow-right-24: 打开 tf2 教程](https://docs.ros.org/en/jazzy/Tutorials/Intermediate/Tf2/Tf2-Main.html)

-   :material-cube-scan:{ .lg .middle } __仿真 · MuJoCo 文档__

    ---

    E3 用到的时间步、控制与接触接口，按需查阅。

    [:octicons-arrow-right-24: 打开文档](https://mujoco.readthedocs.io/)

</div>

## 跟随实验

<div class="grid cards" markdown>

-   :material-hammer-wrench:{ .lg .middle } __E1 · 小海龟到点控制__

    ---

    订阅位置、发布速度，闭合比例反馈；多个起点重复并保存日志。

    [:octicons-arrow-right-24: 查看步骤](../PROJECT-BASED-LEARNING/project-based-learning.md)

-   :material-hammer-wrench:{ .lg .middle } __P1 接入 ROS 2__

    ---

    把 6DoF 位姿估计包成节点，订阅 RGB-D、发布物体 tf，在 RViz 显示。

    [:octicons-arrow-right-24: 查看项目](../PROJECT-BASED-LEARNING/pose-6dof.md)

</div>

## 自测

??? question "相机图像应该用话题、服务还是动作？为什么？"
    话题。图像是持续的数据流，不需要逐条确认；用服务会让发送方等待，用动作则没有意义。

??? question "订阅方一直收不到数据，但话题列表里能看到，最可能是什么？"
    QoS 不兼容（例如发布方 best-effort、订阅方要求 reliable），或者话题类型不一致、命名空间不对。

??? question "融合两个传感器时，tf 查询报外推错误，说明什么？"
    查询的时间戳超出了已有变换的时间范围，通常是时间戳不同步或某个变换发布频率太低。

??? question "为什么调试要先录 rosbag？"
    固定输入后，每次改算法都在同一数据上比较，结果可复现；也能判断问题出在传感器还是算法。

## 相关专题

<div class="grid cards" markdown>

-   __空间与本体__

    ---

    坐标系和 URDF 的数学基础。

    [:octicons-arrow-right-24: 进入](../EMBODIMENT/embodiment.md)

-   __仿真与虚实迁移__

    ---

    仿真里跑通之后，怎样到真机。

    [:octicons-arrow-right-24: 进入](../SIM-TO-REAL/sim-to-real.md)

</div>
