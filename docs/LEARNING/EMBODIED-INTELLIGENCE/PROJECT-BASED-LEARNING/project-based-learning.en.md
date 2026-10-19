# Embodied Intelligence · Guided Projects

!!! abstract "Explore systems before learning policies"
    E0–E2 build intuition; E3 connects data, training and execution. E4–E5 are electives. Follow upstream installation/code for your selected version; use the learning steps and acceptance criteria below.

| Project | Prerequisite | Environment | Depth |
|---|---|---|---|
| E0 · Task dissection | Introductory concepts | Browser, videos or demos | Orientation |
| E1 · ROS 2 turtlesim | Python and terminal basics | Supported ROS 2 system, CPU | Fundamentals |
| E2 · Small algorithms | Python, frames and arrays | CPU | Fundamentals |
| E3 · Simulated imitation | F1-equivalent PyTorch and observation/action semantics | Graphics; GPU or tutorial training notebook | Intermediate entry |
| E4 · Language and action | E3 and a fixed body/task | Skill-calling or VLA-dependent | Elective |
| E5 · Deployment and transfer | Working, evaluated policy | Target hardware and interfaces | Elective |
| [P1 · 6DoF pose estimation (hands-on)](pose-6dof.md) | Camera model, rigid transforms, PyTorch basics | RGB-D dataset or a LiDAR iPhone; GPU training, Jetson deployment | Perception track |

## E0

### Dissect a robot task

**Follow:** [Modern Robotics](https://modernrobotics.northwestern.edu/) Preview/chapter introductions; consult the selected [LeRobot tutorial](https://huggingface.co/spaces/lerobot/robot-learning-tutorial) for learning concepts.

Choose point navigation or manipulation. Identify the body, sensors, goal, observations, state, actions and success conditions; draw the feedback loop.

**Acceptance:** a one-page system explanation distinguishing perception, planning, control and policies.

## E1

### Turtlesim: messages and feedback

**Primary:** [ROS 2 Tutorials](https://docs.ros.org/en/jazzy/Tutorials.html), Beginner CLI Tools (turtlesim, nodes, topics), then Client Libraries (Python publisher/subscriber). Match the distribution; [FishROS](https://fishros.com/d2lros2/) is an optional Chinese explanation aid.

1. Run the original example, use keyboard control and inspect messages.
2. Publish velocity and identify units and frequency.
3. Subscribe to pose; log goal distance and heading error.
4. Turn toward the goal, then use simple proportional feedback; slow near the goal.
5. Add tolerance, velocity limits, timeout and stopping on missing pose messages.
6. Repeat from several initial positions and save trajectories/logs.

**Personal change:** replace fixed velocity with position feedback.

**Acceptance:** at least five starting points, endpoint errors and timeouts, plus an explanation of open-loop commands vs. feedback. Turtlesim demonstrates software and planar control, not physical robot dynamics or manipulation.

## E2

### Planning or arm motion: choose one first

**Follow:** [PythonRobotics](https://github.com/AtsushiSakai/PythonRobotics) / [documentation](https://atsushisakai.github.io/PythonRobotics/index.html), a collection of visual algorithm examples.

| Route | Example | Experiment |
|---|---|---|
| Mobile robot | [A* grid search](https://atsushisakai.github.io/PythonRobotics/modules/5_path_planning/grid_base_search/grid_base_search.html) | Vary goals, obstacles and robot radius; measure feasibility, path length and runtime |
| Arm | [Two-link example under Arm Navigation](https://atsushisakai.github.io/PythonRobotics/modules/7_arm_navigation/arm_navigation.html) | Vary link lengths and targets; compare joints/endpoints and unreachable goals |

1. Run the original and identify inputs, outputs, units and frames.
2. Change one condition at a time; save animations and results.
3. Create and explain a failure.

**Acceptance:** three conditions and one failure analysis. Explore both routes if useful; choose one body/task for a substantial project.

## E3

### Learn actions from simulated demonstrations

**Follow:** [LeRobot imitation learning in simulation](https://huggingface.co/docs/lerobot/main/il_sim) / [GitHub](https://github.com/huggingface/lerobot), with [ACT concepts](https://huggingface.co/docs/lerobot/act).

!!! note "Check the environment first"
    The simulation tutorial is on main. Match source and dependencies and record the commit. Verify graphics and keyboard/gamepad interaction first. Use CPU for initial environment/data checks; arrange GPU or the linked notebook for training. Allow installation/debugging time at this later stage.

1. Teleoperate the original task using keyboard or gamepad.
2. Record and inspect alignment of images, state, actions and time.
3. Reserve evaluation data by complete episode, not randomly mixed neighboring frames.
4. Train ACT; save configurations, curves and checkpoints.
5. Define success, failure, timeout and termination before autonomous execution.
6. Run at least 20 fixed evaluations; report successes/attempts, initial conditions and failure types.
7. Change demonstration count or coverage and compare under the same evaluation conditions.

**Personal change:** design a data-coverage experiment and explain cases with low training loss but failed execution.

**Acceptance:** data documentation, versions, training logs, videos, success rate and failures. Report offline action errors separately from task success; simulation results do not establish hardware performance.

## E4

### Connect language: choose one route

| Route | Reference | Sequence | Acceptance |
|---|---|---|---|
| Skill calling | [Applications & agents](../../FOUNDATION-MODELS/APPLICATIONS/applications.md), existing E1/E3 interface | Bounded skills/arguments → rule baseline → model output → validation and feedback | Normal, ambiguous and unsupported instructions; separate interpretation and execution success |
| VLA adaptation | [SmolVLA](https://huggingface.co/docs/lerobot/smolvla) | Check robot/compute/data support → collect or convert → fine-tune → fixed evaluation | Document data and action interfaces; test changed task conditions |

SmolVLA involves robot configuration and data; arbitrary E3 data is not automatically compatible. Check cameras, states, actions, annotations and controllers. If resources are limited, complete skill calling or study the VLA architecture first.

Both routes require closed-loop validation; a VLA is not necessarily a separate high-level planner.

## E5

### Edge deployment and sim-to-real

**References:** keep the E3 environment in [Simulation & sim-to-real](../SIM-TO-REAL/sim-to-real.md) and consult the selected robot's interfaces. Add platforms only after defining the transfer problem.

Fix the model/task. Measure software-loop latency, then connect target hardware in the loop. Check frames, units, control frequency, action limits and stopping. Validate gradually on hardware under controlled conditions and record domain gaps/recovery.

**Acceptance:** latency, resources, task success and transfer comparisons. Distinguish software-in-the-loop, hardware-in-the-loop and physical robot experiments; simulation alone does not complete transfer.

## Common record

Keep tutorial, versions, run instructions, task definitions, data sources, results, failures and personal changes. Label future work as planned and report only measured outcomes.

[Roadmap](../ROADMAP/roadmap.md)
