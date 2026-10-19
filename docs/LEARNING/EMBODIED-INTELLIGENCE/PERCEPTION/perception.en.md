# Perception & State Estimation

!!! abstract "What this page solves"
    What the robot sees, where it is, and where the target is. Perception turns raw sensor data into quantities you can plan and control with; state estimation fuses noisy, incomplete observations into trustworthy estimates. **Required at entry level, depth by task**: manipulation goes "detection → 6DoF pose → grasping", mobility goes "odometry → fusion → SLAM".

## How far to go

| Level | Learn | Done when |
|---|---|---|
| Orientation | Images, depth, point clouds, detection, pose, localisation and mapping | You can say what observations mobile and grasping tasks each need |
| Foundation | Camera model, calibration, back-projection, filtering | Derive pixel → 3D point; explain one Kalman step |
| Advanced (manipulation) | Detection / segmentation, registration, 6DoF pose, grasping | P1 done with benchmark numbers |
| Advanced (mobility) | Visual / lidar odometry, multi-sensor fusion, SLAM | Run an open-source odometry and evaluate trajectory error |

## Core knowledge

**Camera model and calibration**: the pinhole model projects camera coordinates to pixels with intrinsics K (focal lengths fx, fy, principal point cx, cy); distortion parameters undo lens bending. Extrinsics are the camera pose relative to the robot or world; hand–eye calibration finds the camera-to-arm extrinsics. Calibration error becomes pose error directly.

**Depth and point clouds**: RGB-D cameras, stereo or lidar give depth. Depth z plus pixel (u, v) back-projects through the intrinsics to a 3D point. Common point-cloud operations: voxel downsampling, normal estimation, outlier removal, plane segmentation.

**Time synchronisation**: sensors sample at different moments and must be aligned before fusion. A few tens of milliseconds of offset is centimetres of error on a moving robot.

**Filtering and the Kalman family**: each Kalman step predicts with a motion model and corrects with the observation, weighting them by covariance; the EKF linearises non-linear models, and the error-state Kalman filter (ESKF) is the usual form for IMU fusion. Particle filters represent arbitrary distributions with samples and suit multi-modal problems such as global localisation.

**Odometry and SLAM**: odometry estimates relative motion between frames and drifts; SLAM builds a map while localising and removes drift with loop closure and graph optimisation. Visual-inertial (VIO) and lidar-inertial (LIO) odometry are the mainstream today.

**Object pose and affordances**: a 6DoF pose is the object's rotation and translation relative to the camera — the precondition for grasping. Geometric methods (registration, PnP) are interpretable and need no training; learned methods (e.g. FoundationPose) are stronger under occlusion and symmetry. Affordances answer "where can it be grasped and how is it used".

**3D scene representations**: occupancy grids (navigation), point clouds and meshes (geometry), NeRF and 3D Gaussian splatting (photorealistic rendering), semantic maps (object-level understanding). The downstream task decides.

## Debate · explicit geometry vs. learned implicit representations

Geometry: calibration, registration and optimisation are interpretable and verifiable, and you know why they fail. Learning: end-to-end networks are stronger under occlusion, weak texture and category generalisation, and foundation models handle new objects zero-shot. Industry is mostly hybrid: learning gives the initial guess or semantics, geometry refines and verifies.

## Curated resources

<div class="grid cards" markdown>

-   :material-video:{ .lg .middle } __Main course · Cyrill Stachniss__

    ---

    Probabilistic estimation, Kalman filters, ICP and SLAM, with slides and videos; choose lectures by task.

    [:octicons-arrow-right-24: Open channel](https://www.youtube.com/@CyrillStachniss)

-   :material-cube-outline:{ .lg .middle } __Practice · Open3D tutorials__

    ---

    RGB-D to point clouds, ICP and global registration; used directly in P1.

    [:octicons-arrow-right-24: Open tutorials](https://www.open3d.org/docs/release/tutorial/pipelines/icp_registration.html)

-   :material-trophy-outline:{ .lg .middle } __Benchmark · BOP__

    ---

    Unified data format, metrics and leaderboards for 6DoF pose estimation.

    [:octicons-arrow-right-24: Open BOP](https://bop.felk.cvut.cz/)

</div>

## Follow-along experiment

<div class="grid cards" markdown>

-   :material-hammer-wrench:{ .lg .middle } __P1 · RGB-D 6DoF pose estimation__

    ---

    Segmentation → back-projection → registration → evaluation, then Jetson deployment and an iPhone capture.

    [:octicons-arrow-right-24: See project](../PROJECT-BASED-LEARNING/pose-6dof.md)

</div>

## Self-check

??? question "Given pixel (u, v) and depth z, how do you get the 3D point in the camera frame?"
    x = (u − cx)·z / fx, y = (v − cy)·z / fy, z unchanged.

??? question "What happens if the observation noise covariance in a Kalman filter is set very large?"
    The filter trusts the prediction and nearly ignores observations: smooth but slow to respond. Very small makes it follow the observations closely and noisily.

??? question "Why does odometry drift, and how does SLAM fix it?"
    Small relative-motion errors accumulate. SLAM recognises revisited places (loop closure) and spreads the accumulated error over the whole trajectory with graph optimisation.

??? question "What is special about pose estimation for symmetric objects (e.g. a can)?"
    Rotations about the symmetry axis are indistinguishable from appearance, so several poses are equivalent; evaluate with ADD-S and train with symmetry-aware losses, or the network learns an average.

## Related topics

<div class="grid cards" markdown>

-   __Embodiment & Space__

    ---

    The maths behind intrinsics, extrinsics and poses.

    [:octicons-arrow-right-24: Open](../EMBODIMENT/embodiment.md)

-   __Vision & 3D Space (AI)__

    ---

    The model side: detection, segmentation, NeRF / 3DGS.

    [:octicons-arrow-right-24: Open](../../AI/CV/cv.md)

</div>
