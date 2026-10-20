# Robot Toolchain & Ecosystem

!!! abstract "What this page solves"
    Robotics job ads are full of names: ROS 2, Isaac Lab, MoveIt, Nav2, FAST-LIO, Pinocchio, EtherCAT… This page puts each back into its layer of the robot loop and says what it is, what it solves and how far to learn it. The topic pages explain principles; this page covers **the tools industry uses to put them into practice**.

## Term map

| Layer | Common terms | What it is | Suggested depth |
|---|---|---|---|
| Middleware | ROS 2 (Humble / Jazzy), DDS / QoS, launch, rosbag2 | Node communication and system organisation; DDS is the transport, QoS sets reliability and freshness | **Must**: nodes, topics, services, actions, parameters, launch |
| Description & frames | URDF / Xacro, tf2, RViz | Robot model and transform tree, visual debugging | **Must** |
| Control interface | ros2_control, micro-ROS | Hardware abstraction and controller management; MCUs inside ROS 2 | Know; micro-ROS ties directly to your MCU background |
| Motion planning | MoveIt 2, OMPL | Arm IK, collision checking and trajectory planning | Run one arm-planning example |
| Navigation | Nav2, costmap, AMCL | Mobile-robot localisation, cost maps, global / local planning | Run one simulated navigation |
| Simulation | Gazebo, Isaac Sim, Isaac Lab, MuJoCo / MJX, PyBullet, Genesis | Gazebo pairs with ROS; Isaac Sim for realistic rendering and sensors; Isaac Lab and MuJoCo for massively parallel RL | Master one: Isaac Lab or MuJoCo for RL, Gazebo for ROS systems |
| LiDAR / visual SLAM | FAST-LIO2, LIO-SAM, ORB-SLAM3, VINS-Fusion | LiDAR-inertial and visual-inertial odometry and mapping | Explain front end, back end and loop closure; run one open system |
| Optimisation & estimation | GTSAM / factor graphs, Ceres, EKF | Back-end optimisation and filtering | Understand the principles, use one library |
| Point clouds & geometry | PCL, Open3D, ICP, hand-eye calibration | Point-cloud processing, registration, camera-robot calibration | Your 6DoF project already covers the basics |
| Dynamics & control | Pinocchio, WBC, MPC (OCS2 / Crocoddyl), ZMP, impedance control | Rigid-body dynamics; humanoid whole-body control; biped balance criterion; compliant contact | Concepts; deeper only for control roles |
| Legged locomotion learning | legged_gym, Isaac Lab locomotion, domain randomisation | Train gaits with RL in simulation, then transfer to hardware | Know the pipeline; explain the keys to sim-to-real |
| Learned manipulation | LeRobot, ACT, Diffusion Policy, OpenVLA / π0 | Open implementations of imitation learning and VLA | Run one LeRobot training and evaluation |
| Data collection | Teleoperation (leader-follower / VR), ALOHA, UMI | Collecting high-quality demonstrations | Know; central for embodied-data roles |
| Hardware & comms | EtherCAT, CAN / CAN FD, Jetson Orin / Thor, Isaac ROS, TensorRT | Real-time fieldbus; joint communication; edge compute and GPU-accelerated ROS 2 packages | Embedded background helps; explain real-time behaviour |

## Core knowledge

**ROS 2 communication**: topics (publish / subscribe, continuous data), services (request / response, short operations), actions (long tasks with feedback and cancellation). Underneath is DDS; QoS decides best-effort vs reliable and keep-latest vs keep-history. Sensor streams usually use best effort, control commands reliable.

**tf2 and URDF**: URDF describes links and joints; `robot_state_publisher` turns joint angles into a transform tree; tf2 keeps frame transforms over time. Any "wrong coordinates" problem starts with the tf tree.

**Who does what in simulation**: Gazebo integrates best with ROS 2 for system validation; Isaac Sim renders realistically for synthetic data and vision; Isaac Lab and MuJoCo (MJX on GPU / TPU) run thousands of parallel environments and are the workhorses for training RL policies.

**SLAM structure**: the front end associates data (feature matching or point-cloud registration), the back end estimates state (filtering or factor-graph optimisation), and loop closure corrects accumulated drift. FAST-LIO2 is tightly coupled LiDAR-inertial filtering, LIO-SAM uses a factor graph, ORB-SLAM3 is feature-based visual SLAM.

**Humanoid and legged control**: the classical route is MPC plus whole-body control (WBC, a quadratic program allocating joint torques, contact forces and task goals); the learning route trains policies with RL and domain randomisation in Isaac Lab, then transfers to hardware. The two are often combined.

**Real time and comms**: EtherCAT passes one frame through every slave in sequence, with sub-millisecond cycles, and is the mainstream joint-drive bus; control loops must separate hard from soft real time, often with PREEMPT_RT on Linux.

## Curated resources

<div class="grid cards" markdown>

-   :material-robot-industrial:{ .lg .middle } __Middleware · ROS 2 tutorials__

    ---

    From CLI tools, nodes and topics to launch, tf2 and URDF, beginner to intermediate.

    [:octicons-arrow-right-24: Open tutorials](https://docs.ros.org/en/jazzy/Tutorials.html)

-   :material-gamepad-variant:{ .lg .middle } __Simulation & RL · Isaac Lab__

    ---

    NVIDIA's robot-learning framework with ready locomotion and manipulation tasks and parallel training.

    [:octicons-arrow-right-24: Open docs](https://isaac-sim.github.io/IsaacLab/)

-   :material-cube-outline:{ .lg .middle } __Physics simulation · MuJoCo__

    ---

    Accurate, lightweight contact dynamics; the workhorse of academia and DeepMind.

    [:octicons-arrow-right-24: Open docs](https://mujoco.readthedocs.io/)

</div>

## Follow-along · a minimal ROS 2 system chain

1. Describe a two-link arm in URDF, show it in RViz and move joints with `joint_state_publisher_gui`.
2. Write a node that publishes a joint trajectory; record and replay it with rosbag2.
3. Publish the 6DoF pose project's output as a tf and see the object frame in RViz.
4. Optional: run one legged-locomotion example task in Isaac Lab and record the training curve.

**Done when**: node graph (`rqt_graph`), tf tree and replay video; you can explain the QoS of each topic.

## Self-check

??? question "When do you use a topic, a service or an action?"
    Topics for continuous streams (sensors, state); services for short request-response (queries, triggers); actions for long tasks needing feedback and cancellation (navigate to a goal, execute a grasp).

??? question "Camera data shows up in the wrong place in RViz. What do you check first?"
    The tf tree: is the camera frame published, are parent-child relations and extrinsics right, are timestamps aligned; then the message frame_id.

??? question "Why is RL usually trained in Isaac Lab / MuJoCo rather than Gazebo?"
    RL needs huge amounts of interaction; Isaac Lab and MuJoCo (MJX) run thousands of environments in parallel on GPU, while Gazebo targets single-environment system validation and lacks the throughput.

??? question "Main difference between FAST-LIO2 and LIO-SAM?"
    FAST-LIO2 tightly couples LiDAR and IMU with an iterated EKF and registers raw points directly, so it is fast; LIO-SAM smooths with a factor graph, making GPS and loop-closure constraints easy to add.

??? question "What problem does whole-body control (WBC) solve?"
    A humanoid has several tasks to satisfy at once (balance, hand position, joint limits, contact-force constraints); WBC uses a quadratic program to find joint torques or accelerations that meet them by priority.

## Related topics

<div class="grid cards" markdown>

-   __Software Systems & Simulation__

    ---

    A systems view of ROS 2 and simulators.

    [:octicons-arrow-right-24: Open](../SYSTEMS/systems.md)

-   __Perception & State Estimation__

    ---

    Principles of filtering, SLAM and 6DoF pose.

    [:octicons-arrow-right-24: Open](../PERCEPTION/perception.md)

-   __Planning & Control__

    ---

    Principles of MPC, whole-body and impedance control.

    [:octicons-arrow-right-24: Open](../PLANNING-CONTROL/planning-control.md)

</div>
