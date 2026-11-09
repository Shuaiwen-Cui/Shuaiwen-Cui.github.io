# __PH.D. STUDIES__

![IOT-SHM](IOT-SHM.jpg){width=100%}

## __Research Overview__

!!! info "RESEARCH TOPICS -- PUSHING STRUCTURAL HEALTH MONITORING TO THE EDGE DEVICES"
    My Ph.D. research topic is __Distributed Edge Intelligence Framework for IoT-based Structural Health Monitoring__, with a focus on developing a distributed edge-intelligence enabling framework for resource-constrained, low-cost edge devices, to promote a paradigm shift from cloud to edge and from knowledge-driven to AI-driven (especially agent-driven) approaches, and to apply it to typical IoT-based structural health monitoring applications.

!!! info "KEY TASKS"
    - IoT System Development: Edge Device (MCU Level, Low-cost & Resource-constrained) and Cloud Platform Development
    - Distributed Edge Intelligence Enabling Framework Development: TinySHM (Vector & Matrix Operations / Digital Signal Processing / Machine Learning & Artificial Intelligence Libraries)
    - Structural Health Monitoring Applications: Measurement / System Identification / Damage Detection / Damage Localization / Damage Evaluation

<div class="rs-road" markdown>

<div class="rs-step" markdown>

![Platform](ROAD-PLATFORM.png)

<div class="rs-step-k">Foundation</div>
<div class="rs-step-t">LiftNode + LiftHub</div>
<div class="rs-step-d">Where edge intelligence runs</div>

</div>

<div class="rs-step" markdown>

![Framework](ROAD-FRAMEWORK.png)

<div class="rs-step-k">Foundation</div>
<div class="rs-step-t">TinySHM</div>
<div class="rs-step-d">How edge computation is organized</div>

</div>

<div class="rs-step" markdown>

![Level 1](ROAD-L1.png)

<div class="rs-step-k">Level 1</div>
<div class="rs-step-t">Single-Node Intelligence</div>
<div class="rs-step-d">Decisions on one node</div>

</div>

<div class="rs-step" markdown>

![Level 2](ROAD-L2.png)

<div class="rs-step-k">Level 2</div>
<div class="rs-step-t">Multi-Node Cooperation</div>
<div class="rs-step-d">Rule-based coordination</div>

</div>

<div class="rs-step" markdown>

![Level 3](ROAD-L3.png)

<div class="rs-step-k">Level 3</div>
<div class="rs-step-t">Multi-Node Collaboration</div>
<div class="rs-step-d">Coordination learned from data</div>

</div>

</div>

<div class="rs-two" markdown>

<figure class="rs-fig" markdown="span">
  [![Paradigm shift](SUMMARY.jpg)](SUMMARY.jpg)
  <figcaption>From centralized to distributed computing</figcaption>
</figure>

<figure class="rs-fig" markdown="span">
  [![Review graphical abstract](REVIEW-GA.jpg)](REVIEW-GA.jpg)
  <figcaption>Graphical abstract of the edge-to-cloud review</figcaption>
</figure>

</div>

<div class="grid cards" markdown>

-   :material-file:{ .lg .middle } __Review Paper - Ubiquitous Computing and Intelligence in IoT-based Structural Health Monitoring (From Edge to Cloud)__

    ---

    **Cui, S.**, Fu, Y.*, Fu, H., & Shen, W. (2026). Edge-to-Cloud Computing and Intelligence for IoT-based Structural Health Monitoring: A Comprehensive Review. Advanced Engineering Informatics. 

    [:octicons-arrow-right-24: <a href="https://doi.org/10.1016/j.aei.2025.104300" target="_blank"> DOI </a>](#)

</div>

## __I IoT-based SHM System Development__

The platform connects sensor nodes, a cloud service and the operator. It is a research prototype, built to the level needed to implement and measure edge-intelligence methods.

<figure class="rs-fig" markdown="span">
  [![Platform architecture](PLATFORM-ARCH.png)](PLATFORM-ARCH.png)
  <figcaption>LiftNode and LiftHub platform architecture</figcaption>
</figure>

### __1.1 Edge Device__

LiftNode pairs an always-on sentinel accelerometer with a 20-bit acquisition accelerometer, so the node can wait cheaply and record well. We developed MCU nodes on both STM32 and ESP32; the current focus is ESP32.

<div class="rs-stats six">
  <div class="rs-stat"><b>ESP32-S3</b><span>2 cores, 240 MHz, FPU, SIMD</span></div>
  <div class="rs-stat"><b>8 MB</b><span>PSRAM, with 16 MB flash</span></div>
  <div class="rs-stat"><b>20-bit</b><span>ADXL355, 22.5 µg/√Hz</span></div>
  <div class="rs-stat"><b>0.89 µA</b><span>ADXL367 sentinel (datasheet)</span></div>
  <div class="rs-stat"><b>10.5 Ah</b><span>3 × 18650, MPPT solar input</span></div>
  <div class="rs-stat"><b>≈ USD 140</b><span>bill of materials per node</span></div>
</div>

<figure class="rs-fig" markdown="span">
  [![LiftNode composition](LIFTNODE-PARTS.jpg)](LIFTNODE-PARTS.jpg)
  <figcaption>LiftNode composition: assembly, functional architecture and assembly sequence</figcaption>
</figure>

<figure class="rs-fig" markdown="span">
  [![LiftNode design features](LIFTNODE-FEATURES.jpg)](LIFTNODE-FEATURES.jpg)
  <figcaption>Design features. Currents in panel (c) are datasheet values, not whole-board measurements.</figcaption>
</figure>

<div class="rs-two" markdown>

<figure class="rs-fig" markdown="span">
  [![Node design](Node-Overview.png)](Node-Overview.png)
  <figcaption>Node design</figcaption>
</figure>

<figure class="rs-fig" markdown="span">
  [![Wireless sensor network](WSN.jpg)](WSN.jpg)
  <figcaption>Assembled nodes</figcaption>
</figure>

<figure class="rs-fig" markdown="span">
  [![PCB design](PCB.png)](PCB.png)
  <figcaption>PCB design</figcaption>
</figure>

<figure class="rs-fig" markdown="span">
  [![PCB](PCB.jpg)](PCB.jpg)
  <figcaption>Fabricated boards</figcaption>
</figure>

</div>

<div class="rs-two">

<div class="rs-video">
<p>Sensor Streaming Demo</p>
<iframe src="https://www.youtube-nocookie.com/embed/f1rnn7kyeOM" frameborder="0" allowfullscreen></iframe>
</div>

<div class="rs-video">
<p>Sensor Synchronization Demo</p>
<iframe src="https://www.youtube-nocookie.com/embed/IdZ6cGDO1go" frameborder="0" allowfullscreen></iframe>
</div>

</div>

<div class="grid cards" markdown>

-   :simple-github:{ .lg .middle } __LiftNode_ESP32🎯🏆__

    ---

    MCU IoT Node with High Performance Edge Computing, Based on ESP32

    <!-- ROLLBACK (github link hidden): [:octicons-arrow-right-24: <a href="https://github.com/Shuaiwen-Cui/LiftNode_ESP32.git" target="_blank"> Code </a>](#) -->

    [:octicons-arrow-right-24: <a href="http://www.cuishuaiwen.com:8100/" target="_blank"> Online Doc (Hosted on Personal Server) </a>](#)

    <!-- ROLLBACK (github link hidden): [:octicons-arrow-right-24: <a href="https://shuaiwen-cui.github.io/LiftNode_ESP32/" target="_blank"> Online Doc (Hosted on Github Pages) </a>](#) -->

</div>

### __1.2 Cloud Platform__

LiftHub forwards MQTT messages, stores data and provides digital-twin visualization. It also exposes the nodes to conversational-agent control. It is an enabling layer for research demonstrations; the research focus stays on the edge.

<div class="grid cards" markdown>

-   :simple-github:{ .lg .middle } __LiftHub 🎯🏆__

    ---

    Acts as the server-side to interact with LiftNode, providing data storage and analysis functions.

    <!-- ROLLBACK (github link hidden): [:octicons-arrow-right-24: <a href="https://github.com/Shuaiwen-Cui/LiftHub.git" target="_blank"> Code </a>](#) -->

    [:octicons-arrow-right-24: <a href="http://www.cuishuaiwen.com:8200/" target="_blank"> Online Doc (Hosted on Personal Server) </a>](#)

    <!-- ROLLBACK (github link hidden): [:octicons-arrow-right-24: <a href="https://shuaiwen-cui.github.io/LiftHub/" target="_blank"> Online Doc (Hosted on Github Pages) </a>](#) -->

</div>

## __II Distributed Edge Intelligence Framework for Resource-Constrained Devices__

TinySHM is MCU-native middleware that sits between the board drivers and SHM applications. It is not a general TinyML library: learning is one capability alongside acquisition, identification and damage screening. It runs without Linux or Python.

<div class="rs-stats">
  <div class="rs-stat"><b>9</b><span>modules in 4 groups</span></div>
  <div class="rs-stat"><b>≈ 410</b><span>public API functions</span></div>
  <div class="rs-stat"><b>≈ 20.5k</b><span>lines of implementation</span></div>
  <div class="rs-stat"><b>165</b><span>unit-test functions</span></div>
  <div class="rs-stat"><b>6 + 2</b><span>modal methods + pairwise variants</span></div>
</div>

<figure class="rs-fig" markdown="span">
  [![TinySHM stack and modules](TINYSHM-STACK.png)](TINYSHM-STACK.png)
  <figcaption>LiftNode firmware stack and the nine TinySHM modules</figcaption>
</figure>

<div class="grid cards" markdown>

-   :material-puzzle-outline:{ .lg .middle } __Composable__

    ---

    Each module offers a portable C/C++ baseline, with hardware acceleration added where available.

-   :material-transit-connection-variant:{ .lg .middle } __End to End__

    ---

    Modules span the pipeline from acquisition through processing, identification and damage screening.

-   :material-chip:{ .lg .middle } __MCU First__

    ---

    Float32 numerics, Jacobi SVD and pre-derived backpropagation follow directly from the hardware limits.

-   :material-access-point-network:{ .lg .middle } __IoT Native__

    ---

    Functions are exposed as commands over ESP-NOW and MQTT, so a node can be reconfigured at runtime.

</div>

### __2.1 Nine Modules__

<div class="rs-layer">
  <div class="rs-layer-name">Task Orchestration</div>
  <div class="rs-mods">
    <div class="rs-mod"><b>TinyBench</b><span>Named data slots, typed tool registry, remote query</span></div>
    <div class="rs-mod"><b>TinyOrch</b><span>Named flows of typed steps, run by command</span></div>
  </div>
</div>
<div class="rs-layer">
  <div class="rs-layer-name">SHM Specific</div>
  <div class="rs-mods">
    <div class="rs-mod"><b>TinyMeasurement</b><span>Streaming, burst and real-time acquisition; clock sync</span></div>
    <div class="rs-mod"><b>TinySysID</b><span>PP, FDD, EFDD, ITD, ERA, SSI; pairwise; voting</span></div>
    <div class="rs-mod"><b>TinyDamage</b><span>Damage index, COMAC localization, severity</span></div>
  </div>
</div>
<div class="rs-layer">
  <div class="rs-layer-name">General Computation</div>
  <div class="rs-mods">
    <div class="rs-mod"><b>TinyMath</b><span>Linear algebra, decompositions, eigensolvers</span></div>
    <div class="rs-mod"><b>TinyDSP</b><span>FFT, filter design, wavelets, ICA</span></div>
    <div class="rs-mod"><b>TinyAI</b><span>Layers, training, inference, quantization</span></div>
  </div>
</div>
<div class="rs-layer">
  <div class="rs-layer-name">Platform Adaptation</div>
  <div class="rs-mods">
    <div class="rs-mod"><b>TinyToolbox</b><span>ESP-IDF adaptation, node identity, microsecond time, boot self-test</span></div>
  </div>
</div>

<div class="rs-two" markdown>

<figure class="rs-fig" markdown="span">
  [![TinyMath](MOD-MATH.png)](MOD-MATH.png)
  <figcaption><b>TinyMath.</b> LU, Cholesky and QR; Jacobi and randomized SVD; shifted QR; Arnoldi and Lanczos.</figcaption>
</figure>

<figure class="rs-fig" markdown="span">
  [![TinyDSP](MOD-DSP.png)](MOD-DSP.png)
  <figcaption><b>TinyDSP.</b> FFT, FIR and IIR design, correlation, resampling, wavelets and ICA.</figcaption>
</figure>

<figure class="rs-fig" markdown="span">
  [![TinyMeasurement](MOD-MEASURE.png)](MOD-MEASURE.png)
  <figcaption><b>TinyMeasurement.</b> Streaming, burst and real-time acquisition, with gateway-beacon clock synchronization.</figcaption>
</figure>

<figure class="rs-fig" markdown="span">
  [![TinyDamage](MOD-DAMAGE.png)](MOD-DAMAGE.png)
  <figcaption><b>TinyDamage.</b> Screening against a healthy baseline: detect, locate, rate.</figcaption>
</figure>

<figure class="rs-fig" markdown="span">
  [![TinyBench and TinyOrch](MOD-ORCH.png)](MOD-ORCH.png)
  <figcaption><b>TinyBench + TinyOrch.</b> Pipelines as data: flows are defined and run by command, without reflashing.</figcaption>
</figure>

<figure class="rs-fig" markdown="span">
  [![Lines of code per module](TINYSHM-SCALE.png)](TINYSHM-SCALE.png)
  <figcaption><b>Scale of the code.</b> Approximate counts from the source tree; vendor libraries excluded.</figcaption>
</figure>

</div>

### __2.2 On-Device Training: TinyAI__

Inference-only TinyML leaves model maintenance outside the node. TinyAI runs forward propagation, loss evaluation, backpropagation and parameter updates on the MCU. Each operator carries an analytically derived backward method, so no automatic-differentiation graph is needed at runtime.

<div class="rs-split" markdown>

<figure class="rs-fig" markdown="span">
  [![TinyAI](MOD-AI.png)](MOD-AI.png)
</figure>

<div markdown>

<div class="rs-stats">
  <div class="rs-stat"><b>0.30 ms</b><span>forward pass, 16-64-8 MLP</span></div>
  <div class="rs-stat"><b>18.8 ms</b><span>training epoch, 150 samples</span></div>
</div>

Layers cover dense, 1-D and 2-D convolution, pooling, normalization and multi-head self-attention. An inference-only build compiles out the backward paths.

</div>

</div>

<figure class="rs-fig" markdown="span">
  [![TinyAI training](TINYAI-TRAIN.png)](TINYAI-TRAIN.png)
  <figcaption>TinyAI training step, illustrative SRAM/PSRAM placement and optimizer buffers</figcaption>
</figure>

### __2.3 On-Node System Identification: TinySysID__

TinySysID moves output-only operational modal analysis onto the MCU. Six method families and eight variants share one interface, so their estimates can be compared directly.

<div class="rs-split" markdown>

<figure class="rs-fig" markdown="span">
  [![TinySysID](MOD-SYSID.png)](MOD-SYSID.png)
</figure>

<div markdown>

<div class="rs-stats">
  <div class="rs-stat"><b>4–47 KiB</b><span>on-chip peak heap across variants</span></div>
  <div class="rs-stat"><b>&lt; 0.4–7.0 s</b><span>per call, 5 channels, 60 s at 50 Hz</span></div>
</div>

Pairwise processing lowers peak memory by solving channel pairs in sequence. It trades away joint spatial information, and runtime does not always improve.

</div>

</div>

<figure class="rs-fig" markdown="span">
  [![Eigen pipeline](EIGEN-PIPELINE.png)](EIGEN-PIPELINE.png)
  <figcaption>Modal parameters are eigen-quantities: TinyMath provides a solver for each matrix type</figcaption>
</figure>

<figure class="rs-fig" markdown="span">
  [![TinySysID pairwise](TINYSYSID-PAIRS.png)](TINYSYSID-PAIRS.png)
  <figcaption>Pairwise division, resource cost and pair-and-order voting (configuration-specific library benchmark)</figcaption>
</figure>

### __2.4 Compared with Embedded Software__

<div class="rs-table">
<table>
  <thead>
    <tr><th>Software</th><th>SHM Tasks</th><th>Scope</th><th>Modal ID</th><th>On-Node AI</th><th>Workflow Control</th></tr>
  </thead>
  <tbody>
    <tr><td>CMSIS-DSP</td><td>Not supported</td><td>Numerical kernels</td><td>Not supported</td><td>Classical inference</td><td>N.A.</td></tr>
    <tr><td>ESP-DSP</td><td>Not supported</td><td>DSP kernels</td><td>Not supported</td><td>N.A.</td><td>N.A.</td></tr>
    <tr><td>TFLite Micro</td><td>Not supported</td><td>Neural inference</td><td>Not supported</td><td>Inference only</td><td>N.A.</td></tr>
    <tr><td>ISHMP Toolsuite</td><td>ID, detection</td><td>SHM services</td><td>ERA</td><td>N.A.</td><td>App. commands</td></tr>
    <tr class="rs-hl"><td>TinySHM</td><td>Acquisition, ID, assessment</td><td>SHM pipeline</td><td>Six methods</td><td>Inference and training</td><td>Command-defined flows</td></tr>
  </tbody>
</table>
</div>

<div class="rs-scope"><strong>Scope.</strong> The table compares built-in, documented functions, not accuracy or maturity. TinySHM is evaluated on ESP32-S3 through workload replays. Large-network, long-term field and whole-board energy validation are not established. The SSI routines use a zero-lag block construction that differs from standard SSI-COV.</div>

<div class="grid cards" markdown>

-   :simple-github:{ .lg .middle } __TinySHM🎯🏆__

    ---

    Currently Supported Platforms:

    - ESP32

    <!-- ROLLBACK (github link hidden): [:octicons-arrow-right-24: <a href="https://github.com/Shuaiwen-Cui/TinySHM.git" target="_blank"> Code </a>](#) -->

    [:octicons-arrow-right-24: <a href="http://www.cuishuaiwen.com:8300/" target="_blank"> Online Doc (Hosted on Personal Server) </a>](#)

    <!-- ROLLBACK (github link hidden): [:octicons-arrow-right-24: <a href="https://shuaiwen-cui.github.io/TinySHM/" target="_blank"> Online Doc (Hosted on Github Pages) </a>](#) -->

</div>

## __III Emerging ML/AI Powered Typical Structural Health Monitoring Applications__

The applications progress from decisions on a single node, to rule-based cooperation across nodes, to collaboration learned from data.

<div class="rs-table">
<table>
  <thead>
    <tr><th>SHM Task</th><th>Single-Node</th><th>Multi-Node Cooperation (Rule-Based)</th><th>Multi-Node Collaboration (Learned)</th></tr>
  </thead>
  <tbody>
    <tr><td><strong>Measurement</strong></td><td>Adaptive Trigger Sensing; On-Node Learning</td><td>-</td><td>-</td></tr>
    <tr><td><strong>System Identification</strong></td><td>On-Board Modal Identification</td><td>-</td><td>LCSI</td></tr>
    <tr><td><strong>Damage Assessment</strong></td><td>-</td><td>Rapid Condition Assessment</td><td>Manifold Dynamics</td></tr>
  </tbody>
</table>
</div>

### __3.1 Single-Node Intelligence__

<div class="rs-level" markdown>

![Level 1](ROAD-L1.png)

The evidence for the decision lies in one node's own record, so the decision stays on that node.

</div>

<div class="rs-entry" markdown>

<div class="rs-head">
  <span class="rs-title">Smart Adaptive Trigger Sensing Powered by Edge Intelligence and Digital Twin</span>
  <span class="rs-pill pub">Published</span>
</div>

<div class="rs-tags">
  <span class="rs-tag on">Measurement</span><span class="rs-tag">System Identification</span><span class="rs-tag">Damage Detection</span><span class="rs-tag">Damage Localization</span><span class="rs-tag">Damage Evaluation</span>
</div>
<div class="rs-tags">
  <span class="rs-tag kw">Trigger Sensing</span><span class="rs-tag kw">Closed-Loop Feedback Control</span><span class="rs-tag kw">Edge Intelligence</span><span class="rs-tag kw">Bayesian Optimization</span>
</div>

<p class="rs-lead">A simple threshold-and-duration trigger stays cheap during surveillance. Compact surrogates and Bayesian optimization recalibrate its two parameters when needed.</p>

<figure class="rs-fig" markdown="span">
  [![SATM graphical abstract](SATM.jpg)](SATM.jpg)
</figure>

<div class="rs-split" markdown>

<figure class="rs-fig" markdown="span">
  [![SATM testbed](SATM-SETUP.jpg)](SATM-SETUP.jpg)
  <figcaption>Laboratory testbed and event types</figcaption>
</figure>

<div markdown>

<div class="rs-stats">
  <div class="rs-stat"><b>20</b><span>appended Stage II evaluations</span></div>
  <div class="rs-stat"><b>1,530</b><span>grid points avoided</span></div>
  <div class="rs-stat"><b>0.9396</b><span>surrogate F3 (grid: 0.9351)</span></div>
</div>

<figure class="rs-fig" markdown="span">
  [![SATM loop](SATM-LOOP.png)](SATM-LOOP.png)
  <figcaption>Feedback loop of the adaptive trigger</figcaption>
</figure>

</div>

</div>

<div class="rs-scope"><strong>Scope.</strong> Validated in the laboratory with one node and one acceleration axis on a five-storey frame. Long-term field validation, net energy saving, a complete MCU feedback-loop trace and the benefit of field recalibration are not established.</div>

<div class="grid cards" markdown>

-   :material-file:{ .lg .middle } __Journal Paper - Mechanical System and Signal Processing__

    ---

    **Cui, S.**, Fu, Y.*, Fu, H., Yu, X., & Shen, W. (2025). Smart Adaptive Trigger Sensing Powered by Edge Intelligence and Digital Twin for Energy-Efficient Wireless Structural Health Monitoring.  Mechanical Systems and Signal Processing, Volume 241, 2025, 113537.

    [:octicons-arrow-right-24: <a href="https://doi.org/10.1016/j.ymssp.2025.113537" target="_blank"> DOI </a>](#)

-   :material-file:{ .lg .middle } __Conference Paper - SHMII-13__

    ---

    **Cui, S.**, Yu, X., & Fu, Y.* (2025). Smart adaptive triggering strategy for edge intelligence enabled energy-efficient sensing. In *Proceedings of the 13th International Conference on Structural Health Monitoring of Intelligent Infrastructure (SHMII-13)*, pp. 609–616. Graz, Austria: Verlag der TU Graz. (🏆 **Best Conference Paper Award** 1st/202)

    [:octicons-arrow-right-24: <a href="https://doi.org/10.3217/978-3-99161-057-1-094" target="_blank"> DOI </a>](#)

-   :material-file:{ .lg .middle } __Singapore Provisional Patent Application - 10202502426R__

    ---

    Adaptive Triggering Mechanism for Time-Series Data Sensing on Edge Devices, Singapore provisional patent application number 10202502426R, 2025

</div>

</div>

<div class="rs-entry" markdown>

<div class="rs-head">
  <span class="rs-title">On-Board Modal Identification and Learning on a Microcontroller Node</span>
  <span class="rs-pill review">Under Review</span>
</div>

<div class="rs-tags">
  <span class="rs-tag on">Measurement</span><span class="rs-tag on">System Identification</span><span class="rs-tag">Damage Detection</span><span class="rs-tag">Damage Localization</span><span class="rs-tag">Damage Evaluation</span>
</div>
<div class="rs-tags">
  <span class="rs-tag kw">Operational Modal Analysis</span><span class="rs-tag kw">On-Device Training</span><span class="rs-tag kw">TinySHM</span><span class="rs-tag kw">LiftNode</span>
</div>

<p class="rs-lead">One ESP32-S3 acts as both sensor and computer: it reproduces workstation modal identification and trains its own event classifier.</p>

<figure class="rs-fig" markdown="span">
  [![Decisions beside data graphical abstract](DECISIONS-GA.jpg)](DECISIONS-GA.jpg)
  <figcaption>In this graphic, 100% and 10.45× refer to the held-out Z24 set with ±10% bands in steady operation. Case 2 photographs: Yu et al. (2026), <a href="https://doi.org/10.1016/j.measurement.2025.119308" target="_blank">Measurement</a>.</figcaption>
</figure>

<div class="rs-stats">
  <div class="rs-stat"><b>957 / 957</b><span>comparable groups agree with PC</span></div>
  <div class="rs-stat"><b>14.5 s</b><span>on-board training, 18,000 windows</span></div>
  <div class="rs-stat"><b>98.86%</b><span>F1 (PC-trained: 98.78%)</span></div>
  <div class="rs-stat"><b>2.33 ms</b><span>features + inference per 1-s window</span></div>
</div>

<div class="rs-split" markdown>

<figure class="rs-fig" markdown="span">
  [![On-board modal identification](ONBOARD-MODAL.jpg)](ONBOARD-MODAL.jpg)
</figure>

<div markdown>

**Modal identification.** Eight estimators run on the board. All 957 target groups with a PC reference (960 planned) meet a frequency difference of at most 2% and a complex MAC of at least 0.95.

**On-node learning.** A 66-parameter classifier trains on recorded metro-tunnel vibration in 40 epochs, with a 46.7 KiB peak TinyAI heap.

</div>

</div>

<figure class="rs-fig" markdown="span">
  [![Tunnel inference](ONBOARD-TUNNEL.png)](ONBOARD-TUNNEL.png)
  <figcaption>On-board inference on recorded tunnel vibration</figcaption>
</figure>

<div class="rs-scope"><strong>Scope.</strong> Both cases replay recorded data. They verify that the board reproduces the workstation, not physical modal accuracy or field damage detection. Tunnel labels come from automatic acceleration and audio thresholds, not annotated train events.</div>

<div class="grid cards" markdown>

-   :material-file:{ .lg .middle } __Manuscript Under Review__

    ---

    **Cui, S.**, Fu, Y.*, Wang, Y., & Fu, H. Decisions beside data: LiftNode and TinySHM for AIoT-native learning and system identification in infrastructure monitoring. Manuscript under review.

</div>

</div>

### __3.2 Multi-Node Cooperation (Rule-Based Coordination)__

<div class="rs-level" markdown>

![Level 2](ROAD-L2.png)

Nodes and gateway divide the work through predefined rules and interfaces.

</div>

<div class="rs-entry" markdown>

<div class="rs-head">
  <span class="rs-title">Adaptive Edge Intelligence for Rapid Structural Condition Assessment Using a Wireless Smart Sensor Network</span>
  <span class="rs-pill pub">Published</span>
</div>

<div class="rs-tags">
  <span class="rs-tag on">Measurement</span><span class="rs-tag">System Identification</span><span class="rs-tag on">Damage Detection</span><span class="rs-tag">Damage Localization</span><span class="rs-tag">Damage Evaluation</span>
</div>
<div class="rs-tags">
  <span class="rs-tag kw">Data-driven Anomaly Detection</span><span class="rs-tag kw">Gaussian Process Regression (GPR)</span><span class="rs-tag kw">Stochastic Process Control (SPC)</span>
</div>

<p class="rs-lead">Nodes compute displacement features, and the gateway combines them through a predefined command workflow.</p>

<figure class="rs-fig" markdown="span">
  [![Rapid condition assessment](PAPER1.png)](PAPER1.png)
</figure>

<p><strong>Contribution.</strong> The reference-free displacement method comes from prior work. Hoang's doctoral research applied it on Xnode, developed the GPR-SPC formulation and collected the railroad-bridge records. My work adapts and optimizes these methods for microcontroller execution and organizes them as a node-gateway workflow.</p>

<div class="rs-stats">
  <div class="rs-stat"><b>8.74 s</b><span>to process a 182.01-s record</span></div>
  <div class="rs-stat"><b>≈ 10 nm</b><span>maximum difference from PC</span></div>
  <div class="rs-stat"><b>1.70 ms</b><span>GPR-SPC numerical kernel</span></div>
  <div class="rs-stat"><b>3,960 B</b><span>dynamic workspace</span></div>
</div>

<div class="rs-two" markdown>

<figure class="rs-fig" markdown="span">
  [![Node-gateway workflow](COOP-WORKFLOW.png)](COOP-WORKFLOW.png)
  <figcaption>Designed node-gateway workflow</figcaption>
</figure>

<figure class="rs-fig" markdown="span">
  [![Displacement replay](DISPLACEMENT-REPLAY.png)](DISPLACEMENT-REPLAY.png)
  <figcaption>LiftNode displacement replay against the PC reference</figcaption>
</figure>

</div>

<div class="rs-scope"><strong>Scope.</strong> These are separate replays of stored records. The bridge campaigns are historical Xnode deployments by Hoang et al. A new field deployment and end-to-end wireless timing are not established.</div>

<div class="grid cards" markdown>

-   :material-file:{ .lg .middle } __Journal Paper - Engineering Structures__

    ---

    **Cui, S.**, Hoang, T., Mechitov, K., Fu, Y.*, & Spencer Jr, B. F. (2025). Adaptive edge intelligence for rapid structural condition assessment using a wireless smart sensor network. Engineering Structures, 326, 119520.

    [:octicons-arrow-right-24: <a href="https://doi.org/10.1016/j.engstruct.2024.119520" target="_blank"> DOI </a>](#)

</div>

</div>

### __3.3 Multi-Node Collaboration (Learned Coordination)__

<div class="rs-level" markdown>

![Level 3](ROAD-L3.png)

The coordination itself is learned from data: first what each node should send, then what state the nodes should share.

</div>

<div class="rs-entry" markdown>

<div class="rs-head">
  <span class="rs-title">Learned-Support Collaborative System Identification (LCSI)</span>
  <span class="rs-pill review">Under Review</span>
</div>

<div class="rs-tags">
  <span class="rs-tag">Measurement</span><span class="rs-tag on">System Identification</span><span class="rs-tag">Damage Detection</span><span class="rs-tag">Damage Localization</span><span class="rs-tag">Damage Evaluation</span>
</div>
<div class="rs-tags">
  <span class="rs-tag kw">Gateway-Learned Spectral Support</span><span class="rs-tag kw">Band Coding</span><span class="rs-tag kw">Equal-Byte Evaluation</span><span class="rs-tag kw">Z24 Bridge</span>
</div>

<p class="rs-lead">The gateway learns which frequency bands the modes occupy and sends them to the nodes. Each node transmits only those DCT coefficients, and paired checks are designed to keep or renew the bands.</p>

<figure class="rs-fig" markdown="span">
  [![LCSI path](LCSI-PATH.png)](LCSI-PATH.png)
</figure>

<div class="rs-stats">
  <div class="rs-stat"><b>10.45×</b><span>smaller per channel (757 of 8,192)</span></div>
  <div class="rs-stat"><b>355 / 355</b><span>modal pairs kept on held-out set</span></div>
  <div class="rs-stat"><b>57–65%</b><span>kept by node-chosen codecs, equal bytes</span></div>
  <div class="rs-stat"><b>−90.4%</b><span>leaf uplink, steady (−88.3% with checks)</span></div>
</div>

<figure class="rs-fig" markdown="span">
  [![LCSI operating cycle](LCSI-CYCLE.png)](LCSI-CYCLE.png)
  <figcaption>Operating cycle: commissioning, stable operation and periodic checks</figcaption>
</figure>

<div class="rs-two" markdown>

<figure class="rs-fig" markdown="span">
  [![Equal-byte comparison](LCSI-EQUAL-BYTES.png)](LCSI-EQUAL-BYTES.png)
  <figcaption>Equal-byte comparison with node-chosen codecs</figcaption>
</figure>

<figure class="rs-fig" markdown="span">
  [![Traffic](LCSI-TRAFFIC.png)](LCSI-TRAFFIC.png)
  <figcaption>Application bytes on each link</figcaption>
</figure>

</div>

<div class="rs-scope"><strong>Scope.</strong> Pre-damage Z24 records are replayed on one ESP32-S3 in sequence, not on a live multi-node network. Savings are application bytes; energy and latency were not measured. Commissioning costs 9.6% more than raw upload and breaks even at window 14. Support renewal under larger structural change has not been tested.</div>

<div class="grid cards" markdown>

-   :material-file:{ .lg .middle } __Manuscript Under Review__

    ---

    Same manuscript as the on-board modal identification entry in Section 3.1.

</div>

</div>

<div class="rs-entry" markdown>

<div class="rs-head">
  <span class="rs-title">Manifold Dynamics: A Healthy-Only Diagnostic State for Distributed Monitoring</span>
  <span class="rs-pill thesis">Ph.D. Thesis · Unpublished</span>
</div>

<div class="rs-tags">
  <span class="rs-tag">Measurement</span><span class="rs-tag on">System Identification</span><span class="rs-tag on">Damage Detection</span><span class="rs-tag on">Damage Localization</span><span class="rs-tag on">Damage Evaluation</span>
</div>
<div class="rs-tags">
  <span class="rs-tag kw">Healthy-Only Learning</span><span class="rs-tag kw">Latent Flow</span><span class="rs-tag kw">Distance, Direction, Relation, Trend</span><span class="rs-tag kw">Leaf-Gateway Representation</span>
</div>

<p class="rs-lead">Channel-wise codes are fused at the gateway into one diagnostic state, trained without damage labels. Identification, prediction, detection, localization and degradation trend become readouts of that state. Information flows one way, and the models stay frozen after training.</p>

<div class="rs-two" markdown>

<figure class="rs-fig" markdown="span">
  [![Manifold concept](MANIFOLD-CONCEPT.jpg)](MANIFOLD-CONCEPT.jpg)
  <figcaption>Healthy states move along a learned manifold</figcaption>
</figure>

<figure class="rs-fig" markdown="span">
  [![Manifold framework](MANIFOLD-FRAMEWORK.png)](MANIFOLD-FRAMEWORK.png)
  <figcaption>Proposed logical flow, not a measured live multi-board implementation</figcaption>
</figure>

</div>

<div class="rs-stats">
  <div class="rs-stat"><b>4 / 4</b><span>column-removal scenarios localized</span></div>
  <div class="rs-stat"><b>0.731</b><span>window top-1 (healthy false alarm 0.112)</span></div>
  <div class="rs-stat"><b>0.810</b><span>Z24 Stage 02 AUROC (static: 0.683)</span></div>
  <div class="rs-stat"><b>34.552 ms</b><span>per state pair, ESP32-S3 replay</span></div>
</div>

<figure class="rs-fig" markdown="span">
  [![State walks](MANIFOLD-STATES.png)](MANIFOLD-STATES.png)
  <figcaption>State walks on the three validation systems (visualization only)</figcaption>
</figure>

<div class="rs-two" markdown>

<figure class="rs-fig" markdown="span">
  [![NTU localization](MANIFOLD-NTU.png)](MANIFOLD-NTU.png)
  <figcaption>Four-storey frame: localization against label-free baselines</figcaption>
</figure>

<figure class="rs-fig" markdown="span">
  [![Z24 readouts](MANIFOLD-Z24.png)](MANIFOLD-Z24.png)
  <figcaption>Z24: readouts at the earliest recorded departure</figcaption>
</figure>

</div>

<p>On FEMTO bearings, the causal degradation-trend readout reaches a mean monotonicity-trendability score of 0.453. This is not remaining-useful-life prediction.</p>

<div class="rs-scope"><strong>Scope.</strong> Unpublished work reported in my Ph.D. thesis. "Distributed" means a proposed leaf-gateway split, evaluated on centrally acquired and benchmark data rather than a live multi-board network. Z24 Stage 02 includes the settlement-system installation and possible environmental change, so it is not confirmed damage detection. The board result uses preloaded features.</div>

<div class="grid cards" markdown>

-   :material-file:{ .lg .middle } __Ph.D. Thesis - Chapter 6 (Unpublished)__

    ---

    **Cui, S.** Distributed Edge Intelligence for IoT-Based Wireless Structural Health Monitoring. Ph.D. Thesis, Nanyang Technological University.

</div>

</div>

## __IV Teaching and Pedagogical Research__

### ArduinoNode

A simple WSN node developed for SHM teaching, based on Arduino UNO R4 WIFI, successfully used in the APESS2025 Summer School held at The Hong Kong Polytechnic University.

![](apess2025.jpg)

<div class="grid cards" markdown>

-   :simple-github:{ .lg .middle } __APESS2025_ARDUINONODE__

    ---

    SHM Teaching WSN Node Based on Arduino UNO R4 WIFI

    [:octicons-arrow-right-24: <a href="https://github.com/Shuaiwen-Cui/APESS2025_ArduinoNode.git" target="_blank"> Code </a>](#)

    [:octicons-arrow-right-24: <a href="https://shuaiwen-cui.github.io/APESS2025_ArduinoNode/" target="_blank"> Online Doc </a>](#)

-   :material-file:{ .lg .middle } __Journal Paper - IEEE Transactions on Education__

    ---

    **Cui, S.**, Fu, Y.*, Xia, Y., Zhang Q., & Li, S. (2026). A Class-Lab-Field Pedagogical Framework for Structural Health Monitoring using Ultra-Low-Cost Wireless IoT Prototypes. IEEE Transactions on Education. 

    [:octicons-arrow-right-24: <a href="https://doi.org/10.1109/TE.2026.3717462" target="_blank"> DOI </a>](#)

</div>
