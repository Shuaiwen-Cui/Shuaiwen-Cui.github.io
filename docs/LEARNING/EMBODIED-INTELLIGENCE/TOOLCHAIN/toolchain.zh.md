# 机器人工具链与生态

!!! abstract "这一页解决什么"
    机器人岗位的 JD 里满是专有名词：ROS 2、Isaac Lab、MoveIt、Nav2、FAST-LIO、Pinocchio、EtherCAT……这一页把它们放回机器人闭环的各层，说清每个名词是什么、解决什么、学到什么程度。前面各专题讲原理，这一页讲<strong>工业界用什么工具把原理落地</strong>。

## 名词地图

| 层次 | 常见名词 | 是什么 | 建议程度 |
|---|---|---|---|
| 中间件 | ROS 2（Humble / Jazzy）、DDS / QoS、launch、rosbag2 | 节点通信与系统组织；DDS 是底层传输，QoS 决定可靠性与时效 | **必会**：节点、话题、服务、动作、参数、launch |
| 描述与坐标 | URDF / Xacro、tf2、RViz | 机器人模型与坐标树，可视化调试 | **必会** |
| 控制接口 | ros2_control、micro-ROS | 硬件抽象与控制器管理；MCU 接入 ROS 2 | 了解；micro-ROS 与你的 MCU 背景直接相关 |
| 运动规划 | MoveIt 2、OMPL | 机械臂逆运动学、碰撞检测与轨迹规划 | 跑通一个机械臂规划示例 |
| 导航 | Nav2、costmap、AMCL | 移动机器人定位、代价地图、全局 / 局部规划 | 跑通一次仿真导航 |
| 仿真 | Gazebo、Isaac Sim、Isaac Lab、MuJoCo / MJX、PyBullet、Genesis | Gazebo 配合 ROS；Isaac Sim 做逼真渲染与传感器；Isaac Lab 与 MuJoCo 做大规模并行强化学习 | 选一个用熟：做 RL 选 Isaac Lab 或 MuJoCo，做 ROS 系统选 Gazebo |
| 激光 / 视觉 SLAM | FAST-LIO2、LIO-SAM、ORB-SLAM3、VINS-Fusion | 激光-惯性、视觉-惯性里程计与建图 | 能讲清前端、后端与回环；跑通一个开源系统 |
| 优化与估计 | GTSAM / 因子图、Ceres、EKF | 后端优化与滤波 | 懂原理，会用一个库 |
| 点云与几何 | PCL、Open3D、ICP、手眼标定 | 点云处理、配准与相机-机器人标定 | 你已有 6DoF 项目基础 |
| 动力学与控制 | Pinocchio、WBC、MPC（OCS2 / Crocoddyl）、ZMP、阻抗控制 | 刚体动力学计算；人形全身控制；双足平衡判据；接触柔顺 | 懂概念；运控方向再深入 |
| 足式运动学习 | legged_gym、Isaac Lab locomotion、域随机化 | 仿真中用强化学习训练步态再上真机 | 了解流程，能讲 sim-to-real 的关键 |
| 学习型操作 | LeRobot、ACT、Diffusion Policy、OpenVLA / π0 | 模仿学习与 VLA 的开源实现 | 跑通一次 LeRobot 训练与评估 |
| 数据采集 | 遥操作（主从臂 / VR）、ALOHA、UMI | 收集高质量演示数据 | 了解；具身数据岗重点 |
| 硬件与通信 | EtherCAT、CAN / CAN FD、Jetson Orin / Thor、Isaac ROS、TensorRT | 实时现场总线；关节通信；端侧算力与 GPU 加速的 ROS 2 包 | 嵌入式背景加分，能讲清实时性 |

## 核心知识

**ROS 2 的通信**：话题（发布 / 订阅，持续数据）、服务（请求 / 应答，短操作）、动作（长任务，带反馈和取消）。底层是 DDS，QoS 决定「尽力而为还是可靠送达、只保留最新还是保留历史」，传感器流常用 best effort，控制指令要 reliable。

**tf2 与 URDF**：URDF 描述连杆与关节，`robot_state_publisher` 根据关节角发布坐标树；tf2 维护各坐标系随时间的变换，任何「坐标错了」的问题先查 tf 树。

**仿真器的分工**：Gazebo 与 ROS 2 集成最好，适合验证系统；Isaac Sim 渲染逼真，适合合成数据和视觉；Isaac Lab、MuJoCo（MJX 跑在 GPU / TPU）支持成千上万个环境并行，是训练 RL 策略的主力。

**SLAM 的结构**：前端做数据关联（特征匹配或点云配准），后端做状态估计（滤波或因子图优化），回环检测纠正累计漂移。FAST-LIO2 是紧耦合的激光-惯性滤波方法，LIO-SAM 用因子图，ORB-SLAM3 是特征点视觉 SLAM。

**人形与足式控制**：传统路线是模型预测控制 + 全身控制（WBC，用二次规划在关节力矩、接触力和任务目标之间分配）；学习路线是在 Isaac Lab 里用强化学习 + 域随机化训练策略，再迁移到真机。两者常结合使用。

**实时与通信**：EtherCAT 用一帧数据串行经过所有从站，周期可达毫秒以下，是关节驱动的主流总线；控制回路要区分硬实时和软实时，Linux 上常配合 PREEMPT_RT。

## 精选资源

<div class="grid cards" markdown>

-   :material-robot-industrial:{ .lg .middle } __中间件 · ROS 2 官方教程__

    ---

    从 CLI 工具、节点话题到 launch、tf2、URDF，按初级 → 中级顺序走。

    [:octicons-arrow-right-24: 打开教程](https://docs.ros.org/en/jazzy/Tutorials.html)

-   :material-gamepad-variant:{ .lg .middle } __仿真与强化学习 · Isaac Lab__

    ---

    NVIDIA 的机器人学习框架，含足式运动、操作等现成任务与并行训练。

    [:octicons-arrow-right-24: 打开文档](https://isaac-sim.github.io/IsaacLab/)

-   :material-cube-outline:{ .lg .middle } __物理仿真 · MuJoCo__

    ---

    接触动力学精确、轻量，学术界与 DeepMind 的主力仿真器。

    [:octicons-arrow-right-24: 打开文档](https://mujoco.readthedocs.io/)

</div>

## 跟随实验 · 一条最小的 ROS 2 系统链路

1. 用 URDF 描述一个两连杆机械臂，在 RViz 中显示，并用 `joint_state_publisher_gui` 拖动关节。
2. 写一个节点发布关节轨迹，录 rosbag2 再回放。
3. 把 6DoF 位姿项目的输出作为 tf 发布，在 RViz 里看到物体坐标系。
4. 选做：在 Isaac Lab 跑通一个足式运动示例任务，记录训练曲线。

**验收**：节点图（`rqt_graph`）、tf 树、回放视频；能解释每个话题用了什么 QoS。

## 自测

??? question "话题、服务、动作分别用在什么场景？"
    话题传持续数据流（传感器、状态）；服务做短的请求应答（查询、触发）；动作做长任务，需要中途反馈和取消（导航到某点、执行抓取）。

??? question "相机数据在 RViz 里位置不对，先查什么？"
    先查 tf 树：相机坐标系是否发布、父子关系和外参是否正确、时间戳是否对齐；再查消息里的 frame_id。

??? question "为什么 RL 训练常用 Isaac Lab / MuJoCo 而不是 Gazebo？"
    强化学习需要海量交互，Isaac Lab 和 MuJoCo（MJX）能在 GPU 上并行成千上万个环境；Gazebo 偏单环境系统验证，吞吐不够。

??? question "FAST-LIO2 和 LIO-SAM 的主要区别？"
    FAST-LIO2 用迭代扩展卡尔曼滤波紧耦合激光与 IMU，直接配准原始点，速度快；LIO-SAM 用因子图平滑，便于加入 GPS 和回环约束。

??? question "全身控制（WBC）解决什么问题？"
    人形机器人有多个同时要满足的任务（平衡、手的位置、关节限位、接触力约束），WBC 用二次规划按优先级求出满足约束的关节力矩或加速度。

## 相关专题

<div class="grid cards" markdown>

-   __软件系统与仿真__

    ---

    ROS 2 与仿真器的系统观。

    [:octicons-arrow-right-24: 进入](../SYSTEMS/systems.md)

-   __感知与状态估计__

    ---

    滤波、SLAM 与 6DoF 位姿的原理。

    [:octicons-arrow-right-24: 进入](../PERCEPTION/perception.md)

-   __规划与控制__

    ---

    MPC、全身控制与阻抗控制的原理。

    [:octicons-arrow-right-24: 进入](../PLANNING-CONTROL/planning-control.md)

</div>
