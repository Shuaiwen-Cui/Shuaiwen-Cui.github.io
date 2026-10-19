# Planning & Control

!!! abstract "What this page solves"
    Planning answers "which route and which trajectory"; control answers "how to correct errors while executing". This is the most classical, interpretable part of robotics and the baseline learned policies are compared against. **Required at entry level**; your feedback-control and system-identification background transfers directly.

## How far to go

| Level | Learn | Done when |
|---|---|---|
| Orientation | Paths vs trajectories, planning vs control, open vs closed loop | Explain "choosing a route" vs "correcting execution error" |
| Foundation | A*, PID and trajectory tracking | E1 and one E2 experiment done, with failures reported |
| Advanced | Sampling planners, trajectory optimisation, MPC, force / impedance control, whole-body control | Go deep on one, chosen by body and task |

## Core knowledge

**Paths and trajectories**: a path is a geometric sequence of points; a trajectory adds time (velocity, acceleration). Planners usually find a collision-free path first, then time-parameterise it within velocity and acceleration limits.

**Search-based planning**: A* finds shortest paths on grids or graphs and is optimal if the heuristic never overestimates. Great for 2D navigation; the state space explodes in high dimensions.

**Sampling-based planning**: RRT and PRM sample the continuous space and connect samples into a tree or graph — suited to high-dimensional arm configuration spaces; RRT* is asymptotically optimal. Results need smoothing.

**Trajectory optimisation**: treat the trajectory as decision variables, minimise a cost (time, energy, smoothness) under constraints (collision, joint limits). CHOMP and TrajOpt are typical; sensitive to initialisation.

**Feedback control**: PID acts on the error, its integral and its rate — the most common baseline. Trajectory tracking = feed-forward (from the desired trajectory) + feedback (correcting deviation).

**Model predictive control (MPC)**: each cycle, predict ahead with a model, solve a constrained optimisation online, execute only the first step, re-solve next cycle. Handles constraints explicitly at a computational cost; used widely in legged robots and autonomous driving.

**Force and impedance control**: contact tasks (peg-in-hole, wiping) cannot be position-controlled alone. Impedance control makes the end-effector behave like a spring-damper, compliant on contact instead of pushing hard.

**Whole-body control**: humanoids and legged robots must satisfy balance, contacts and several task goals at once, often as a prioritised quadratic programme.

## Debate · model-driven vs. learning-driven

Model-driven (planning + MPC) is interpretable, enforces safety constraints and needs little data, but requires accurate models and struggles with complex contact and perception uncertainty. Learning-driven (imitation, RL) captures behaviours that are hard to model, but needs lots of data and fails opaquely. The trend is layering: learning decides what to do, model-based controllers do it safely.

## Curated resources

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Main course · Modern Robotics__

    ---

    Ch. 9 trajectories, ch. 10 motion planning, ch. 11 control (incl. force control); pick chapters by task.

    [:octicons-arrow-right-24: Open course site](https://modernrobotics.northwestern.edu/)

-   :material-language-python:{ .lg .middle } __Example library · PythonRobotics__

    ---

    Small visual examples of A*, RRT, MPC and more; each runs as is and is easy to modify.

    [:octicons-arrow-right-24: Open docs](https://atsushisakai.github.io/PythonRobotics/index.html)

</div>

## Follow-along experiments

<div class="grid cards" markdown>

-   :material-hammer-wrench:{ .lg .middle } __E1 · Go-to-goal feedback__

    ---

    Replace constant velocity with position feedback; record final error and timeouts.

    [:octicons-arrow-right-24: See steps](../PROJECT-BASED-LEARNING/project-based-learning.md)

-   :material-hammer-wrench:{ .lg .middle } __E2 · A* or two-link arm__

    ---

    Change obstacles, start/goal or target; create one failure and explain it.

    [:octicons-arrow-right-24: See steps](../PROJECT-BASED-LEARNING/project-based-learning.md)

</div>

## Self-check

??? question "Why must the A* heuristic not overestimate?"
    Overestimating can make the search abandon the truly optimal path early, losing the optimality guarantee; an admissible heuristic guarantees the shortest path.

??? question "Why do arm planners use RRT rather than A*?"
    Arm configuration spaces are 6–7 dimensional; gridding them explodes exponentially, while sampling never discretises the whole space.

??? question "Advantages and costs of MPC over PID?"
    Advantages: explicit constraints, prediction with a model, multivariable coordination. Costs: needs a model and an optimisation every cycle, heavy computation and tight real-time demands.

??? question "Why do contact tasks need impedance control?"
    Pure position control turns small position errors into large contact forces that can damage the part or the robot; impedance control makes force proportional to deviation, so contact is compliant.

## Related topics

<div class="grid cards" markdown>

-   __Learned Policies__

    ---

    Learning actions from data, compared and combined with these methods.

    [:octicons-arrow-right-24: Open](../POLICY-LEARNING/policy-learning.md)

-   __Embodiment & Space__

    ---

    Planning happens in joint or task space and relies on kinematics.

    [:octicons-arrow-right-24: Open](../EMBODIMENT/embodiment.md)

</div>
