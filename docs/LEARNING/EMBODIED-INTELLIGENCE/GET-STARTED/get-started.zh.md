# 具身智能 · 两周速成

!!! abstract "这一页和系统学习的分工"
    系统学习（全貌、路线、E0–E5）从控制和策略学习出发，按周慢慢建立完整的机器人知识；这一页只做**两周内能补上的机器人感知最小闭环**：坐标与相机、ROS 2 基本功、一个 6DoF 位姿项目，并把它接进 ROS 2，最后对模仿学习有一次亲手体验。每天约 2–3 小时。

## 目标能力（从常见岗位要求倒推）

机器人感知、具身算法、机器人软件类岗位最常出现的要求：

| 能力 | 常见说法 | 你已有的基础 | 两周内要补到 |
|---|---|---|---|
| 空间表示 | 坐标变换、SE(3)、标定 | 几何、线性代数；硕士阶段三维几何与碰撞检测 | 熟练写出变换链，讲清内外参 |
| 3D 感知 | 检测、6DoF 位姿、点云、三维重建 | 多传感器融合、状态估计、信号处理 | 一个完整的 RGB-D 位姿项目，有基准指标 |
| 机器人软件 | ROS 2、tf、rosbag、RViz | 嵌入式实时系统、多节点通信 | 能写节点、发布 tf、录制回放 |
| 策略学习 | 模仿学习、强化学习、VLA | 闭环控制、贝叶斯优化、数字孪生到现场 | 跑通一次仿真模仿学习并评估成功率 |
| 部署 | 边缘推理、实时性 | **端侧部署是强项** | 在 Jetson 上测感知管线时延 |

## 14 天计划

| 天 | 学什么 | 动手 | 当天产出 |
|---|---|---|---|
| 1 | 坐标系与 SE(3) | Modern Robotics 第 3 章 3.1–3.3；用 numpy 写变换的复合与求逆 | 一页变换链笔记（相机、物体、机器人基座） |
| 2 | 相机模型与点云 | Open3D RGBD 教程：深度图反投影 | 手推并验证像素到三维点 |
| 3 | ROS 2 基本功 | ROS 2 教程：节点、话题、Python 发布订阅 | 一对自己写的发布 / 订阅节点 |
| 4 | tf2 与 rosbag | tf2 教程；录制与回放 | 在 RViz 里看到自己发布的坐标系 |
| 5–6 | 点云配准与位姿评估 | P1 · 6DoF 位姿估计：跑通真值掩码评估，读懂配准代码 | 复现基准结果；讲清 ICP 和 ADD(-S) |
| 7–8 | 学习型分割 | P1：训练分割模型，换掉真值掩码 | 完整系统的指标 |
| 9–10 | 接进 ROS 2 | 把位姿估计包成节点：订阅 RGB-D，发布物体 tf，在 RViz 显示模型 | 一段 RViz 录屏 |
| 11 | 边缘部署 | P1：Jetson 上 TensorRT 测速 | 整条管线的时延分解 |
| 12 | 模仿学习体验 | LeRobot 仿真模仿学习：跑通原例训练与评估（可用云端 GPU） | 固定评估次数下的成功率 |
| 13 | 整理作品 | README + 实拍视频（iPhone LiDAR） | GitHub 仓库首页 |
| 14 | 面试复盘 | 逐条回答下方问题 | 答题笔记 |

## 这两周用到的资源

<div class="grid cards" markdown>

-   :material-open-in-new:{ .lg .middle } __Open3D RGBD 教程__

    ---

    第 2 天 · 相机模型与点云

    [:octicons-arrow-right-24: 打开](https://www.open3d.org/docs/release/tutorial/geometry/rgbd_image.html)

-   :material-open-in-new:{ .lg .middle } __ROS 2 教程__

    ---

    第 3 天 · ROS 2 基本功

    [:octicons-arrow-right-24: 打开](https://docs.ros.org/en/jazzy/Tutorials.html)

-   :material-open-in-new:{ .lg .middle } __Python 发布订阅__

    ---

    第 3 天 · ROS 2 基本功

    [:octicons-arrow-right-24: 打开](https://docs.ros.org/en/jazzy/Tutorials/Beginner-Client-Libraries/Writing-A-Simple-Py-Publisher-And-Subscriber.html)

-   :material-open-in-new:{ .lg .middle } __tf2 教程__

    ---

    第 4 天 · tf2 与 rosbag

    [:octicons-arrow-right-24: 打开](https://docs.ros.org/en/jazzy/Tutorials/Intermediate/Tf2/Tf2-Main.html)

-   :material-open-in-new:{ .lg .middle } __rosbag 录制与回放__

    ---

    第 4 天 · tf2 与 rosbag

    [:octicons-arrow-right-24: 打开](https://docs.ros.org/en/jazzy/Tutorials/Beginner-CLI-Tools/Recording-And-Playing-Back-Data/Recording-And-Playing-Back-Data.html)

-   :material-book-open-variant:{ .lg .middle } __P1 · 6DoF 位姿估计__

    ---

    第 5–6 天 · 点云配准与位姿评估

    [:octicons-arrow-right-24: 打开](../PROJECT-BASED-LEARNING/pose-6dof.md)

-   :material-open-in-new:{ .lg .middle } __LeRobot 仿真模仿学习__

    ---

    第 12 天 · 模仿学习体验

    [:octicons-arrow-right-24: 打开](https://huggingface.co/docs/lerobot/main/il_sim)

</div>


## 作品：RGB-D 物体位姿估计，从基准到 ROS 2 与边缘设备

P1 项目加一个 ROS 2 节点：公开基准上的指标、Jetson 上的时延、RViz 中的实时显示、iPhone 实拍视频。它把感知、机器人软件和端侧部署三类要求连在一起。

## 面试必会（答题要点）

| 问题 | 要点 |
|---|---|
| 旋转有哪些表示？各自的问题？ | 旋转矩阵（9 个数、需正交约束）、欧拉角（万向锁、不连续）、四元数（无奇异、需归一化、q 与 −q 等价）、轴角；网络回归常用 6D 连续表示 |
| 相机内参、外参分别是什么？ | 内参：焦距、主点，把相机坐标投到像素；外参：相机相对世界或机器人的位姿；手眼标定求相机与机械臂之间的外参 |
| ICP 的前提和局限？ | 需要较好的初值，否则陷入局部最优；对离群点和部分重叠敏感；点到面收敛更快；对称物体有多解 |
| 已知 2D–3D 对应怎么求位姿？ | PnP（如 EPnP）+ RANSAC 剔除外点，再用非线性优化细化 |
| ADD 和 ADD-S 的区别？ | ADD 比对应点的平均距离；对称物体用 ADD-S，取最近点距离，避免把对称等价的位姿判错 |
| 多传感器融合要注意什么？ | 外参标定、时间同步、各自噪声模型；卡尔曼滤波按协方差加权预测与观测 |
| 行为克隆为什么会在闭环中失败？ | 误差累积导致进入训练数据没覆盖的状态（分布偏移）；DAgger、扩大数据覆盖、动作分块（ACT）可以缓解 |
| sim-to-real 的常见方法？ | 域随机化、系统辨识缩小差距、在真实数据上微调；先硬件在环再上真机 |
| VLA 是什么？ | 视觉语言模型加动作输出头，用机器人数据训练；关键在动作表示、数据规模和闭环评估 |
| ROS 2 里 tf 树和 QoS 是什么？ | tf 树维护各坐标系之间随时间变化的变换；QoS 规定可靠性、历史深度等，传感器数据常用 best-effort |

## 做完之后

回到 具身路线 的第二轮：感知方向继续深入估计与建图，策略方向接 E3–E4，部署方向接 E5。


<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __全貌（系统学习总览）__

    ---

    系统学习入口

    [:octicons-arrow-right-24: 打开](../OVERVIEW/overview.md)

-   :material-book-open-variant:{ .lg .middle } __具身路线__

    ---

    系统学习入口

    [:octicons-arrow-right-24: 打开](../ROADMAP/roadmap.md)

</div>
