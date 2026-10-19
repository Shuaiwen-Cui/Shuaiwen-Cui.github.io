# Robot Software & Simulation

!!! abstract "What this page solves"
    Connect observation, algorithms and actuation into a system you can run, debug and reproduce. Almost every industry robotics role asks for ROS 2; simulation decides whether you can iterate without a real robot. **Required**; prerequisites: Python, basic Linux.

## How far to go

| Level | Learn | Done when |
|---|---|---|
| Orientation | What nodes, messages, frames and simulators are for | You can draw the node graph of one robot task |
| Foundation | ROS 2 pub/sub, tf2, rosbag, RViz | E1 done; your frames display correctly in RViz |
| Advanced | Simulator interfaces, latency and real time, integration | You can tell "no data", "wrong frame" and "interface mismatch" apart from logs |

## Core knowledge

**Nodes and communication**: ROS 2 splits a system into nodes. Topics are continuous streams (camera images); services are request/response ("reset"); actions are long tasks with feedback ("move to a point"). Choosing the wrong one is a common design error.

**QoS**: sets reliability (reliable / best-effort), history depth and durability. Sensor streams usually use best-effort — dropping a frame beats blocking; commands and maps use reliable. Incompatible publisher and subscriber QoS fails silently.

**tf2 and the frame tree**: keeps every time-varying transform between frames. Asking "where is this camera point in the map frame at time t" multiplies along the tree and interpolates in time. Wrong timestamps are the top cause of broken fusion.

**URDF**: XML describing links, joints, inertia and visuals; the shared input to RViz, kinematics solvers and simulators.

**rosbag record and replay**: record sensor data once, then debug algorithms on identical inputs, separating sensing problems from algorithm problems.

**Simulators**: MuJoCo (accurate contact dynamics, good for manipulation and learning), Isaac Sim / Isaac Lab (GPU-parallel, realistic rendering, large-scale training and synthetic data), Gazebo (tight ROS integration, mobile robots). Pick one, by task.

## Debate · ROS 2 vs. in-house middleware

ROS 2 has a large ecosystem and is what job ads name, but carries learning and performance overhead; many production teams (autonomous driving, consumer robots) build lean middleware to control latency and dependencies. The concepts — nodes, messages, tf, record/replay — transfer either way.

## Curated resources

<div class="grid cards" markdown>

-   :material-robot-industrial:{ .lg .middle } __Main track · ROS 2 tutorials__

    ---

    Beginner: CLI tools → Python publisher/subscriber; Intermediate: tf2. Use the docs for your installed distribution.

    [:octicons-arrow-right-24: Open tutorials](https://docs.ros.org/en/jazzy/Tutorials.html)

-   :material-axis-arrow:{ .lg .middle } __Frames · tf2 tutorials__

    ---

    Broadcasting and listening to transforms, time-travel queries; essential for perception and fusion.

    [:octicons-arrow-right-24: Open tf2 tutorials](https://docs.ros.org/en/jazzy/Tutorials/Intermediate/Tf2/Tf2-Main.html)

-   :material-cube-scan:{ .lg .middle } __Simulation · MuJoCo docs__

    ---

    Time steps, control and contact interfaces used in E3; consult as needed.

    [:octicons-arrow-right-24: Open docs](https://mujoco.readthedocs.io/)

</div>

## Follow-along experiments

<div class="grid cards" markdown>

-   :material-hammer-wrench:{ .lg .middle } __E1 · Turtlesim go-to-goal__

    ---

    Subscribe to position, publish velocity, close a proportional loop; repeat from several starts and keep logs.

    [:octicons-arrow-right-24: See steps](../PROJECT-BASED-LEARNING/project-based-learning.md)

-   :material-hammer-wrench:{ .lg .middle } __P1 into ROS 2__

    ---

    Wrap 6DoF pose estimation as a node: subscribe RGB-D, publish the object tf, show it in RViz.

    [:octicons-arrow-right-24: See project](../PROJECT-BASED-LEARNING/pose-6dof.md)

</div>

## Self-check

??? question "Should camera images use a topic, a service or an action? Why?"
    A topic: images are a continuous stream that needs no per-message acknowledgement; a service would block the sender and an action makes no sense.

??? question "A subscriber never receives data although the topic is listed. Most likely cause?"
    Incompatible QoS (e.g. best-effort publisher, reliable subscriber), a mismatched message type, or a wrong namespace.

??? question "A tf lookup raises an extrapolation error while fusing two sensors. What does it mean?"
    The requested timestamp is outside the time range of the available transforms — usually unsynchronised timestamps or a transform published too rarely.

??? question "Why record a rosbag before debugging?"
    With fixed input every algorithm change is compared on the same data, results are reproducible, and you can tell sensor problems from algorithm problems.

## Related topics

<div class="grid cards" markdown>

-   __Embodiment & Space__

    ---

    The maths behind frames and URDF.

    [:octicons-arrow-right-24: Open](../EMBODIMENT/embodiment.md)

-   __Sim-to-Real__

    ---

    From working in simulation to working on the robot.

    [:octicons-arrow-right-24: Open](../SIM-TO-REAL/sim-to-real.md)

</div>
