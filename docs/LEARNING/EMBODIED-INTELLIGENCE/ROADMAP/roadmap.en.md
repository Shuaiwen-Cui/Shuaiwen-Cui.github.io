# Embodied Intelligence · Roadmap

!!! tip "Short on time?"
    Do the [two-week Get Started](../GET-STARTED/get-started.md) first: it enters through perception and covers frames, ROS 2, [P1 · 6DoF pose](../PROJECT-BASED-LEARNING/pose-6dof.md) and one imitation-learning run; if you pick the perception branch in round two, continue from P1.

Use Modern Robotics as the backbone, then LeRobot for robot learning. Build on your geometry, signals and engineering background; fill specific spatial/interface gaps rather than restart mathematics and programming. Estimates assume 8–10 hours weekly; ROS 2 documentation and examples serve projects, not extra full courses.

## First pass: orientation and a small feedback loop

| Stage | Estimate | Minimum knowledge | Tutorial and project | Completion |
|---|---|---|---|---|
| 0 · Overview and spatial basics | 1–2 weeks | Body, observations, state, frames and feedback | Modern Robotics introduction and selected chapters 2–3; [Overview](../OVERVIEW/overview.md), E0 | Explain the loop, frames and actions |
| 1 · Software exploration | 2–3 weeks | Basic Linux, nodes, messages, position and velocity | [ROS 2 official tutorials](https://docs.ros.org/en/jazzy/Tutorials.html); E1 | Subscribe to pose and publish feedback velocity commands |
| 2 · Algorithm intuition | 1–2 weeks | Frames, kinematics, paths and errors | [PythonRobotics](https://atsushisakai.github.io/PythonRobotics/index.html); E2 | Use the two-link arm with Modern Robotics chapters 4–6; A* is an alternative |
| 3 · Data and policies | 3–5 weeks | PyTorch, demonstrations, observation/action interfaces and success rates | [LeRobot imitation in simulation](https://huggingface.co/docs/lerobot/main/il_sim); E3 | Go from data to a trained policy and repeated closed-loop evaluation |

Complete foundation [F1](../../FOUNDATION-MODELS/PROJECT-BASED-LEARNING/project-based-learning.md#f1) before policy training, or verify that you can read a PyTorch training loop. Initial simulations can use ground-truth state without a full visual SLAM system.

## Second pass: deepen one task

| Branch | Focus | Project extension |
|---|---|---|
| Perception and estimation | Calibration, detection, pose and fusion, chosen for the task | Replace ground truth with sensor estimates and measure errors |
| Planning and control | Tracking, constraints and stability; MPC/force control as needed | Compare controllers under disturbances and delay |
| Software and systems | ROS 2, replay, interfaces, latency and recovery | Make the task reproducible and debuggable |
| Learned policies | Imitation, data coverage, action representations, RL/VLA | E4: adapt and evaluate a language-conditioned policy |
| Deployment and sim-to-real | Identification, domain gaps, runtimes and resources | E5: hardware-in-the-loop, then gradual robot validation |

Choose one small example first, then one body and task for a substantial project. Legged/humanoid robots, dexterous hands, whole-body control, full SLAM and world models are electives.

## Three foundation-model connections

| Route | Output interface | Evaluate |
|---|---|---|
| LLM skill selection | Goals or skill names and parameters | Instruction understanding, interface constraints, feedback and recovery |
| VLA policy | Actions or action sequences | Adaptation data, timing, control interfaces and success |
| World-model-assisted decisions | Predicted states, observations or latent features | Prediction errors, planning benefit and distribution shifts |

Choose one route at a time. Adding an action head alone does not establish a usable robot policy.

## Project records

- Versions: tutorial date, OS, dependencies, simulator, code commit and model/data revisions.
- Task: observation/action semantics, control frequency, success, timeout and termination.
- Evaluation: fixed initial states or seeds, successes/attempts, failure types and latency.
- Contribution: what changed, why, observed results and limitations.

After E1–E3, use target roles to select a branch for a complete project.

[E0–E5 project steps](../PROJECT-BASED-LEARNING/project-based-learning.md) · [Software and simulation](../SYSTEMS/systems.md)
