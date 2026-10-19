# Embodiment & Space

!!! abstract "What this page solves"
    Where the robot is and how it can move. Frames, rotations and forward/inverse kinematics are the shared language of perception, planning, control and policy learning — get one frame wrong and every downstream module is wrong. **Required**; prerequisite: linear algebra (matrix products, orthogonal matrices).

## How far to go

| Level | Learn | Done when |
|---|---|---|
| Orientation | Body, joints, position, orientation, frames | You can mark the base, end-effector and camera frames on a robot drawing |
| Foundation | Homogeneous transforms, rotation representations, FK / IK | E2 two-link example done; explain reachable vs unreachable targets |
| Advanced | Jacobians, singularities, dynamics | Derive and implement for a chosen manipulation or mobile task |

## Core knowledge

**Rigid poses and SE(3)**: an object's pose relative to a frame is a rotation R (3×3 orthogonal, determinant 1) plus a translation t, packed into a 4×4 homogeneous matrix T. Questions like "where is the object the camera sees, in the robot base frame?" are chains of T: T_base_obj = T_base_cam · T_cam_obj.

**Rotation representations**: matrices are explicit but have 9 numbers and must stay orthogonal; Euler angles have 3 numbers but suffer gimbal lock and wrap-around jumps; quaternions have no singularity and interpolate smoothly, with q and −q the same rotation; axis-angle says "how much about which axis". In practice: quaternions for storage and interpolation, matrices for computation, a continuous 6D representation for network regression.

**Forward and inverse kinematics**: FK maps joint angles to the end-effector pose — unique and easy. IK maps a desired pose back to joint angles — it may have no solution (out of reach), several (elbow up / down), or be numerically unstable near singularities.

**Jacobian**: maps joint velocities to end-effector velocities. Where J loses rank the robot is singular: some end-effector motion is impossible whatever the joints do. Velocity control, force control and numerical IK all depend on it.

**Dynamics**: how forces and torques produce accelerations (M(q)q̈ + C(q,q̇)q̇ + g(q) = τ). Not needed for kinematic simulation; needed for force control, fast motion and legged robots.

**Embodiments**: arms (manipulation), mobile bases (navigation), legged and humanoid (terrain and generality), dexterous hands (fine manipulation). The body decides degrees of freedom, available sensors and which algorithms apply.

## Debate · general humanoids vs. specialised bodies

The humanoid camp argues the world is built for humans, so one general body reuses all human data and tools; the specialised camp argues a wheeled base plus an arm beats a humanoid on cost, reliability and energy, and most factory and warehouse tasks need no legs. Both coexist: research and capital chase humanoids, deployments are mostly specialised.

## Curated resources

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Main course · Modern Robotics__

    ---

    Ch. 2–3 for frames and rigid motion, ch. 4–6 for kinematics; book, videos and code share one notation.

    [:octicons-arrow-right-24: Open course site](https://modernrobotics.northwestern.edu/)

-   :material-github:{ .lg .middle } __Companion code · ModernRobotics__

    ---

    Python / MATLAB implementations to check transforms, FK / IK and Jacobians one by one.

    [:octicons-arrow-right-24: Open repository](https://github.com/NxRLab/ModernRobotics)

</div>

## Follow-along experiment

<div class="grid cards" markdown>

-   :material-hammer-wrench:{ .lg .middle } __E2 · Two-link kinematics__

    ---

    Change link lengths and targets, compare joint and end-effector positions, deliberately test an unreachable target.

    [:octicons-arrow-right-24: See steps](../PROJECT-BASED-LEARNING/project-based-learning.md)

</div>

## Self-check

??? question "How are T_cam_obj and T_obj_cam related, and how do you get one from the other?"
    They are inverses. The inverse has rotation Rᵀ and translation −Rᵀt; no general 4×4 inversion needed.

??? question "Why are Euler angles poor for interpolation or network regression?"
    Gimbal lock loses a degree of freedom in some poses, and angles jump at ±180°, giving unnatural interpolation and discontinuous losses.

??? question "Why can IK have several solutions or none?"
    None: the target is outside the workspace. Several: different joint combinations reach the same pose (elbow up / down). Near singularities the Jacobian loses rank and numerical solvers diverge.

??? question "Why check frames before algorithms when debugging a robot?"
    Wrong axis directions, units (m / mm, rad / deg) or transform order look like "the algorithm does not converge", but no algorithm change fixes them. One end-to-end check with a known pose is fastest.

## Related topics

<div class="grid cards" markdown>

-   __Perception & State Estimation__

    ---

    Camera intrinsics, extrinsics and pose estimation build on these transforms.

    [:octicons-arrow-right-24: Open](../PERCEPTION/perception.md)

-   __Planning & Control__

    ---

    Plan and track trajectories in joint or task space.

    [:octicons-arrow-right-24: Open](../PLANNING-CONTROL/planning-control.md)

</div>
