# Embodied Intelligence · Two-Week Get Started

!!! abstract "How this page differs from the systematic track"
    The systematic track (overview, roadmap, E0–E5) starts from control and policy learning and builds full robotics knowledge week by week. This page is **the smallest robot-perception loop you can close in two weeks**: frames and cameras, ROS 2 basics, a 6DoF pose project wired into ROS 2, and one hands-on run of imitation learning. About 2–3 hours a day.

## Target skills (worked back from common role requirements)

| Skill | Usually phrased as | What you already have | Close the gap to |
|---|---|---|---|
| Spatial representation | Transforms, SE(3), calibration | Geometry, linear algebra; 3D geometry and collision detection from the master's | Write transform chains fluently; explain intrinsics and extrinsics |
| 3D perception | Detection, 6DoF pose, point clouds, reconstruction | Multi-sensor fusion, state estimation, signal processing | One complete RGB-D pose project with benchmark numbers |
| Robot software | ROS 2, tf, rosbag, RViz | Real-time embedded systems, multi-node communication | Write nodes, publish tf, record and replay |
| Policy learning | Imitation learning, RL, VLA | Closed-loop control, Bayesian optimisation, digital twin to field | Run one simulated imitation-learning task and report success rate |
| Deployment | Edge inference, real time | **Edge deployment is a strength** | Measure the perception pipeline's latency on a Jetson |

## 14-day plan

| Day | Learn | Do | Output |
|---|---|---|---|
| 1 | Frames and SE(3) | Modern Robotics ch. 3.1–3.3; compose and invert transforms in numpy | One page of transform-chain notes (camera, object, robot base) |
| 2 | Camera model and point clouds | Open3D RGBD tutorial: back-project a depth image | Pixel → 3D point derived and checked |
| 3 | ROS 2 basics | ROS 2 tutorials: nodes, topics, Python pub/sub | Your own publisher / subscriber pair |
| 4 | tf2 and rosbag | tf2 tutorials; record and replay | Your frames visible in RViz |
| 5–6 | Registration and pose metrics | P1 · 6DoF pose estimation: GT-mask evaluation, read the registration code | Benchmark reproduced; explain ICP and ADD(-S) |
| 7–8 | Learned segmentation | P1: train the segmentation model, replace GT masks | Full-system numbers |
| 9–10 | Into ROS 2 | Wrap pose estimation as a node: subscribe RGB-D, publish the object tf, show the model in RViz | An RViz screen recording |
| 11 | Edge deployment | P1: TensorRT timing on the Jetson | Latency breakdown of the whole pipeline |
| 12 | Imitation learning | LeRobot imitation learning in sim: run the original training and evaluation (cloud GPU is fine) | Success rate over a fixed number of trials |
| 13 | Package the project | README + real capture video (iPhone LiDAR) | Repository front page |
| 14 | Interview review | Answer every question below | Answer notes |

## Resources for these two weeks

<div class="grid cards" markdown>

-   :material-open-in-new:{ .lg .middle } __Open3D RGBD tutorial__

    ---

    Day 2 · Camera model and point clouds

    [:octicons-arrow-right-24: Open](https://www.open3d.org/docs/release/tutorial/geometry/rgbd_image.html)

-   :material-open-in-new:{ .lg .middle } __ROS 2 tutorials__

    ---

    Day 3 · ROS 2 basics

    [:octicons-arrow-right-24: Open](https://docs.ros.org/en/jazzy/Tutorials.html)

-   :material-open-in-new:{ .lg .middle } __Python pub/sub__

    ---

    Day 3 · ROS 2 basics

    [:octicons-arrow-right-24: Open](https://docs.ros.org/en/jazzy/Tutorials/Beginner-Client-Libraries/Writing-A-Simple-Py-Publisher-And-Subscriber.html)

-   :material-open-in-new:{ .lg .middle } __tf2 tutorials__

    ---

    Day 4 · tf2 and rosbag

    [:octicons-arrow-right-24: Open](https://docs.ros.org/en/jazzy/Tutorials/Intermediate/Tf2/Tf2-Main.html)

-   :material-open-in-new:{ .lg .middle } __rosbag record & replay__

    ---

    Day 4 · tf2 and rosbag

    [:octicons-arrow-right-24: Open](https://docs.ros.org/en/jazzy/Tutorials/Beginner-CLI-Tools/Recording-And-Playing-Back-Data/Recording-And-Playing-Back-Data.html)

-   :material-book-open-variant:{ .lg .middle } __P1 · 6DoF pose estimation__

    ---

    Day 5–6 · Registration and pose metrics

    [:octicons-arrow-right-24: Open](../PROJECT-BASED-LEARNING/pose-6dof.md)

-   :material-open-in-new:{ .lg .middle } __LeRobot imitation learning in sim__

    ---

    Day 12 · Imitation learning

    [:octicons-arrow-right-24: Open](https://huggingface.co/docs/lerobot/main/il_sim)

</div>


## Project: RGB-D object pose, from benchmark to ROS 2 and the edge

The P1 project plus a ROS 2 node: benchmark accuracy, Jetson latency, live display in RViz and an iPhone capture. It ties perception, robot software and edge deployment together.

## Must-know interview questions (answer key)

| Question | Key points |
|---|---|
| Rotation representations and their pitfalls? | Matrix (9 numbers, orthogonality constraint), Euler angles (gimbal lock, discontinuities), quaternions (no singularity, must normalise, q ≡ −q), axis-angle; networks often regress a continuous 6D representation |
| Intrinsics vs extrinsics? | Intrinsics: focal length and principal point, camera frame → pixels. Extrinsics: camera pose relative to the world or robot; hand–eye calibration finds camera-to-arm extrinsics |
| Assumptions and limits of ICP? | Needs a good initial guess or it falls into local minima; sensitive to outliers and partial overlap; point-to-plane converges faster; symmetric objects have multiple solutions |
| Pose from 2D–3D correspondences? | PnP (e.g. EPnP) + RANSAC for outliers, then non-linear refinement |
| ADD vs ADD-S? | ADD averages distances between corresponding points; symmetric objects use ADD-S (closest-point distance) so symmetry-equivalent poses are not marked wrong |
| What matters in sensor fusion? | Extrinsic calibration, time synchronisation, each sensor's noise model; a Kalman filter weights prediction and measurement by covariance |
| Why does behaviour cloning fail in closed loop? | Compounding errors drive the robot into states the demonstrations never covered (distribution shift); DAgger, wider data coverage, action chunking (ACT) help |
| Common sim-to-real methods? | Domain randomisation, system identification to shrink the gap, fine-tuning on real data; hardware-in-the-loop before the real robot |
| What is a VLA? | A vision-language model with an action head, trained on robot data; the hard parts are action representation, data scale and closed-loop evaluation |
| tf tree and QoS in ROS 2? | The tf tree keeps time-varying transforms between frames; QoS sets reliability, history depth etc., sensor streams usually use best-effort |

## Afterwards

Return to round two of the embodied roadmap: perception → estimation and mapping; policies → E3–E4; deployment → E5.


<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Overview (systematic track)__

    ---

    Systematic learning

    [:octicons-arrow-right-24: Open](../OVERVIEW/overview.md)

-   :material-book-open-variant:{ .lg .middle } __Embodied roadmap__

    ---

    Systematic learning

    [:octicons-arrow-right-24: Open](../ROADMAP/roadmap.md)

</div>
