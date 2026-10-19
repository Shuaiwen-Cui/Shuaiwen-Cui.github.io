# Embodied Intelligence

Start with the robot system: body, observations, decisions and execution. Explore classical methods and learned policies through simulation.

<div class="grid cards" markdown>

-   :material-rocket-launch:{ .lg .middle } __Get Started · Two Weeks__

    ---

    A 14-day robot-perception loop: frames and cameras, ROS 2, a 6DoF pose project wired into ROS 2, and one run of imitation learning.

    [Start here](GET-STARTED/get-started.md)

-   :material-compass-outline:{ .lg .middle } __Overview__

    ---

    Understand robotics, embodied intelligence, imitation learning, reinforcement learning and VLA.

    [Start here](OVERVIEW/overview.md)

-   :material-map-marker-path:{ .lg .middle } __Roadmap__

    ---

    Move from orientation to a simulated control loop, then choose perception, control, policies or systems.

    [Follow the roadmap](ROADMAP/roadmap.md)

-   :material-hammer-wrench:{ .lg .middle } __Guided Projects__

    ---

    Turtlesim, small robot algorithms and simulated imitation learning, with increasing complexity.

    [Start building](PROJECT-BASED-LEARNING/project-based-learning.md)

</div>

## Knowledge map

| Unit | Main question | Depth |
|---|---|---|
| [Space & embodiment](EMBODIMENT/embodiment.md) | What can the body do, and how are frames, states and actions described? | Fundamentals; complex dynamics optional |
| [Software systems & simulation](SYSTEMS/systems.md) | How do data, algorithms and actuators connect? | Fundamentals and practice |
| [Perception & state estimation](PERCEPTION/perception.md) | Where is the robot and what surrounds it? | Orientation, then task-specific depth |
| [Planning & control](PLANNING-CONTROL/planning-control.md) | How are paths chosen and execution errors corrected? | Fundamentals; MPC and whole-body control optional |
| [Learned policies](POLICY-LEARNING/policy-learning.md) | How can demonstrations or interaction teach actions? | Imitation first; RL and VLA advanced |
| [Simulation & sim-to-real](SIM-TO-REAL/sim-to-real.md) | How can simulation experience transfer reliably to hardware? | Simulation early; hardware transfer advanced |

!!! tip "Learning order"
    Run a simple policy and fill in mechanics and control as needed. Frames, observation/action interfaces and basic feedback come first; full SLAM, MPC and whole-body control are not universal prerequisites for introductory imitation learning.

## Connection to foundation models

Language-based skill calling, VLA action generation and world-model prediction are different interfaces. Each must connect to robot state, control and evaluation.

[Foundation overview](../FOUNDATION-MODELS/OVERVIEW/overview.md) · [Multimodality & action interfaces](../FOUNDATION-MODELS/MULTIMODAL/multimodal.md)
