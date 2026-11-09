# __攻博期间研究__

![IOT-SHM](IOT-SHM.jpg){width=100%}

## __研究概览__

!!! info "研究主题 -- 将结构健康监测推向边缘设备"
    我的博士研究课题是 __面向物联网无线结构健康监测的分布式边缘智能__，重点是开发一种分布式边缘智能使能框架，面向资源受限的低成本边缘设备，以促进从云到边和从知识驱动到人工智能驱动（特别是智能体驱动）的范式转变，并将其应用于物联网结构健康监测的典型应用。

!!! info "主要任务"
    - 物联网系统开发：边缘设备（MCU级别，低成本和资源受限）和云平台开发
    - 分布式边缘智能使能框架开发：TinySHM (向量、矩阵运算/数字信号处理/机器学习与人工智能算法库)
    - 结构健康监测应用：测量/系统识别/损伤检测/损伤定位/损伤量化

![REVIEW-GA](REVIEW-GA.jpg){width=100%}

<div class="grid cards" markdown>

-   :material-file:{ .lg .middle } __综述论文-物联网结构健康监测中的泛在计算与智能（从端到云）__

    ---

    **Cui, S.**, Fu, Y.*, Fu, H., & Shen, W. (2026). Edge-to-Cloud Computing and Intelligence for IoT-based Structural Health Monitoring: A Comprehensive Review. Advanced Engineering Informatics. 

    [:octicons-arrow-right-24: <a href="https://doi.org/10.1016/j.aei.2025.104300" target="_blank"> DOI </a>](#)

    [:octicons-arrow-right-24: <a href="https://mp.weixin.qq.com/s/ySQzpHCOhTIcSIN8gqK1bA" target="_blank"> 公众号推送 </a>](#)

</div>

## __I SHM物联网系统开发__

![](platform.png){width=100%}

### __1.1 边缘设备__

!!! note "边缘设备开发"
    为实现边缘智能计算，我们开发了两种类型的 MCU 节点，基于 STM32 和 ESP32。这些节点具有高性能的边缘计算能力，可用于物联网、智能家居、智慧城市等应用场景。目前的开发重点是 ESP32。

<div class="grid cards" markdown>

-   :simple-github:{ .lg .middle } __LiftNode_ESP32🎯🏆__

    ---

    基于 ESP32 的 MCU IoT 节点，具有高性能边缘计算


    <!-- ROLLBACK (github link hidden): [:octicons-arrow-right-24: <a href="https://github.com/Shuaiwen-Cui/LiftNode_ESP32.git" target="_blank"> 代码 </a>](#) -->

    [:octicons-arrow-right-24: <a href="http://www.cuishuaiwen.com:8100/" target="_blank"> 在线文档（个人服务器 托管） </a>](#)

    <!-- ROLLBACK (github link hidden): [:octicons-arrow-right-24: <a href="https://shuaiwen-cui.github.io/LiftNode_ESP32/" target="_blank"> 在线文档（Github 托管） </a>](#) -->

</div>

<!-- <iframe width="800" height="450" src="https://www.youtube-nocookie.com/embed/O2b3-Bjhhws" frameborder="0" allowfullscreen></iframe> -->

<!-- 主控板

![AIoTNode-Base](Main-Board.png){width=100%}

扩展板

![AIoTNode-Extension](Extension-Board.png){width=100%} -->

节点设计

![Node-Design](Node-Overview.png){width=100%}

![WSN](WSN.jpg){width=100%}

PCB 设计

![PCB-Design](PCB.png){width=100%}

![PCB](PCB.jpg){width=100%} 

传感器流式采样Demo

<iframe width="800" height="450" src="//player.bilibili.com/player.html?bvid=BV1Vm7Y6wEH2&autoplay=0" frameborder="0" allowfullscreen></iframe>

传感器同步Demo

<iframe width="800" height="450" src="//player.bilibili.com/player.html?bvid=BV1q75i6qE8r&autoplay=0" frameborder="0" allowfullscreen></iframe>

### __1.2 云平台__

<div class="grid cards" markdown>

-   :simple-github:{ .lg .middle } __LiftHub 🎯🏆__

    ---

    作为服务器端与LiftNode交互，提供数据存储和分析等功能。

    <!-- ROLLBACK (github link hidden): [:octicons-arrow-right-24: <a href="https://github.com/Shuaiwen-Cui/LiftHub.git" target="_blank"> 代码 </a>](#) -->

    [:octicons-arrow-right-24: <a href="http://www.cuishuaiwen.com:8200/" target="_blank"> 在线文档（个人服务器 托管） </a>](#)

    <!-- ROLLBACK (github link hidden): [:octicons-arrow-right-24: <a href="https://shuaiwen-cui.github.io/LiftHub/" target="_blank"> 在线文档（Github 托管） </a>](#) -->

</div>


## __II 资源受限设备分布式边缘智能使能框架__

![](SUMMARY.jpg){width=100%}

<div class="grid cards" markdown>

-   :simple-github:{ .lg .middle } __TinySHM🎯🏆__

    ---

    当前支持平台：

    - ESP32

    <!-- ROLLBACK (github link hidden): [:octicons-arrow-right-24: <a href="https://github.com/Shuaiwen-Cui/TinySHM.git" target="_blank"> 代码 </a>](#) -->

    [:octicons-arrow-right-24: <a href="http://www.cuishuaiwen.com:8300/" target="_blank"> 在线文档（个人服务器 托管） </a>](#)

    <!-- ROLLBACK (github link hidden): [:octicons-arrow-right-24: <a href="https://shuaiwen-cui.github.io/TinySHM/" target="_blank"> 在线文档（Github 托管） </a>](#) -->

</div>

## __III 新兴的机器学习/人工智能驱动的典型结构健康监测应用__

!!! info "三个智能层级"
    应用从单节点决策，推进到多节点间规则驱动的协同，再到从数据中学习的协作。协同依靠预先设计的规则和接口来协调节点，协作则从数据中学习协调方式。

<table>
  <thead>
    <tr><th>SHM任务</th><th>单节点</th><th>多节点协同（规则驱动）</th><th>多节点协作（数据学习驱动）</th></tr>
  </thead>
  <tbody>
    <tr><td><strong>测量</strong></td><td>自适应触发传感；节点端学习</td><td>-</td><td>-</td></tr>
    <tr><td><strong>系统识别</strong></td><td>板端模态识别</td><td>-</td><td>LCSI</td></tr>
    <tr><td><strong>损伤评估</strong></td><td>-</td><td>结构快速状况评估</td><td>流形动力学</td></tr>
  </tbody>
</table>

### __3.1 单节点独立应用__

**- 基于边缘智能与数字孪生的智能自适应触发传感用于能源高效的无线结构健康监测**

🏷️ SHM分类：**测量**/系统识别/损伤检测/损伤定位/损伤量化

- 触发传感
- 闭环反馈控制
- 边缘智能
- 贝叶斯优化

![](SATM.jpg)

**主要结果**

- 在统一的 F3 代理目标下，第二阶段只追加 20 次评估，无需搜索 1,530 点网格，得分 0.9396，网格搜索为 0.9351。

!!! note "适用范围"
    在实验室中用单节点、单轴加速度和五层框架模型验证。长期现场验证、净节能效果、完整的 MCU 反馈闭环记录以及现场重新标定的收益均尚未确立。

![SATM-SETUP](SATM-SETUP.jpg){width=100%}

<div class="grid cards" markdown>

-   :material-file:{ .lg .middle } __会议论文 - 13th International Conference on Structural Health Monitoring of Intelligent Infrastructure (SHMII-13)__

    ---

    **Cui, S.**, Yu, X., & Fu, Y.* (2025). Smart adaptive triggering strategy for edge intelligence enabled energy-efficient sensing. In *Proceedings of the 13th International Conference on Structural Health Monitoring of Intelligent Infrastructure (SHMII-13)*, pp. 609–616. Graz, Austria: Verlag der TU Graz. (🏆 **最佳会议论文** 1st/202)

    [:octicons-arrow-right-24: <a href="https://doi.org/10.3217/978-3-99161-057-1-094" target="_blank"> DOI </a>](#)

</div>

<div class="grid cards" markdown>

-   :material-file:{ .lg .middle } __期刊论文 - Mechanical System and Signal Processing__

    ---

    **Cui, S.**, Fu, Y.*, Fu, H., Yu, X., & Shen, W. (2025). Smart Adaptive Trigger Sensing Powered by Edge Intelligence and Digital Twin for Energy-Efficient Wireless Structural Health Monitoring.  Mechanical Systems and Signal Processing, Volume 241, 2025, 113537.

    [:octicons-arrow-right-24: <a href="https://doi.org/10.1016/j.ymssp.2025.113537" target="_blank"> DOI </a>](#)

    [:octicons-arrow-right-24: <a href="https://mp.weixin.qq.com/s/FnFBCZ_S0t8YfapvUsaKeQ" target="_blank"> 公众号推送 </a>](#)


</div>

<div class="grid cards" markdown>

-   :material-file:{ .lg .middle } __新加坡临时专利申请 - 10202502426R__

    ---

    Adaptive Triggering Mechanism for Time-Series Data Sensing on Edge Devices, 新加坡临时专利申请号 10202502426R, 2025

</div>

**- 微控制器节点上的板端模态识别与学习**

🏷️ SHM分类：**测量**/**系统识别**/损伤检测/损伤定位/损伤量化

- 运行模态分析
- 设备端训练
- TinySHM
- LiftNode

![DECISIONS-GA](DECISIONS-GA.jpg){width=100%}

*图形摘要中的 100% 和 10.45× 对应 Z24 留出数据、±10% 频带下的稳态运行。案例 2 照片来源：Yu et al. (2026), [Measurement](https://doi.org/10.1016/j.measurement.2025.119308)。*

**主要结果**

- 8 种输出型模态识别方法运行在一块 ESP32-S3 LiftNode 上。有 PC 参考的 957 组目标（计划 960 组）全部满足频率差不超过 2%、复数 MAC 不低于 0.95。
- 基于地铁隧道振动记录，一个 66 参数的分类器在板端训练 14.5 s（18,000 个窗口，40 轮，TinyAI 堆峰值 46.7 KiB）。F1 为 98.86%，PC 训练的网络为 98.78%。
- 特征提取加推理每 1 s 窗口耗时 2.33 ms。

!!! note "适用范围"
    两个案例均为已有记录的回放，验证的是板端能复现工作站的计算结果，不代表物理模态精度或现场损伤检测能力。隧道标签来自加速度和音频的自动阈值，并非人工标注的列车事件。

四层框架模型上的板端模态识别

![ONBOARD-MODAL](ONBOARD-MODAL.jpg){width=100%}

隧道振动记录上的板端推理

![ONBOARD-TUNNEL](ONBOARD-TUNNEL.png){width=100%}

<div class="grid cards" markdown>

-   :material-file:{ .lg .middle } __在审稿件__

    ---

    **Cui, S.**, Fu, Y.*, Wang, Y., & Fu, H. Decisions beside data: LiftNode and TinySHM for AIoT-native learning and system identification in infrastructure monitoring. Manuscript under review.

</div>

### __3.2 多节点协同（规则驱动）__

**- 基于无线智能传感网络的自适应边缘智能用于结构快速状况评估**

🏷️ SHM分类：**测量**/系统识别/**损伤检测**/损伤定位/损伤量化

- 数据驱动的异常检测
- 高斯过程回归
- 随机过程控制

![PAPER1](PAPER1.png){width=100%}

**协调方式：** 规则驱动。节点计算位移特征，网关按预先设计的命令流程汇总。

**本人贡献：** 无参考位移方法来自已有工作。Hoang 的博士研究将其应用于 Xnode，提出 GPR-SPC 评估方法，并采集了铁路桥数据。我的工作是把这些方法适配、优化到微控制器上执行，并组织成节点—网关工作流。

**主要结果**

- 当前 LiftNode 处理 182.01 s 的位移记录耗时 8.74 s，每个输出点与 PC 参考的差异都在约 10 nm 以内。
- 30 个训练点、10 个测试点的 GPR-SPC 数值核心耗时 1.70 ms，动态工作区 3,960 B。

!!! note "适用范围"
    以上是对已存记录分别进行的回放。桥梁现场工作是 Hoang 等人早期的 Xnode 部署。新的现场部署和端到端无线时延尚未确立。

LiftNode 位移回放与 PC 参考对比：五秒细节和逐点残差

![DISPLACEMENT-REPLAY](DISPLACEMENT-REPLAY.png){width=100%}

<div class="grid cards" markdown>

-   :material-file:{ .lg .middle } __期刊论文 - Engineering Structures__

    ---

    **Cui, S.**, Hoang, T., Mechitov, K., Fu, Y.*, & Spencer Jr, B. F. (2025). Adaptive edge intelligence for rapid structural condition assessment using a wireless smart sensor network. Engineering Structures, 326, 119520.

    [:octicons-arrow-right-24: <a href="https://doi.org/10.1016/j.engstruct.2024.119520" target="_blank"> DOI </a>](#)

    [:octicons-arrow-right-24: <a href="https://mp.weixin.qq.com/s/KHquagqxXvckCuE57ua8YA" target="_blank"> 公众号推送 </a>](#)

</div>

### __3.3 多节点协作（数据学习驱动）__

**- 学习支撑协同系统识别（LCSI）**

🏷️ SHM分类：测量/**系统识别**/损伤检测/损伤定位/损伤量化

- 网关学习频谱支撑
- 频带编码
- 等字节评估
- Z24 桥

![LCSI-PATH](LCSI-PATH.png){width=100%}

**协调方式：** 学习“每个节点该传什么”。网关从同步记录中学习频谱支撑并下发给节点，配对校验用于决定保留或更新。

**主要结果**

- 每个节点只发送学习频带内的 DCT 系数：8,192 个中的 757 个，单通道缩小 10.45 倍。
- 在 Z24 留出数据上，同样字节数下 LCSI 保留 355/355 个目标模态对，三种由节点自选的编码方式仅保留 57–65%。
- 节点上行流量在稳态下减少 90.4%，在包含一次配对校验的 48 窗口周期内减少 88.3%。

!!! note "适用范围"
    Z24 损伤前记录在一块 ESP32-S3 上依次回放，并非实时多节点网络。节省的是应用层字节，未测量能耗和时延。调试阶段比直接上传原始数据多 9.6%，到第 14 个窗口回本。结构发生较大变化时的频带更新尚未测试。

LCSI 与节点自选编码的等字节对比

![LCSI-EQUAL-BYTES](LCSI-EQUAL-BYTES.png){width=100%}

三条处理路线在各链路上的应用层字节数

![LCSI-TRAFFIC](LCSI-TRAFFIC.png){width=100%}

<div class="grid cards" markdown>

-   :material-file:{ .lg .middle } __在审稿件__

    ---

    与 3.1 节板端模态识别条目为同一稿件。

</div>

**- 流形动力学：面向分布式监测的仅健康数据诊断状态**

🏷️ SHM分类：测量/**系统识别**/**损伤检测**/**损伤定位**/**损伤量化**

- 仅健康数据学习
- 隐空间流
- 距离、方向、关系和趋势读出
- 叶节点—网关表示

![MANIFOLD-CONCEPT](MANIFOLD-CONCEPT.jpg){width=100%}

**协调方式：** 学习共享表示。各通道的编码在网关融合成一个诊断状态。信息从节点单向流向网关，模型训练后冻结。

**主要结果**

- 一个不依赖损伤标签训练的紧凑状态，同时支持识别、预测、检测、定位和退化趋势读出。
- 四层框架：4 种柱移除工况在 8 个随机种子下全部定位正确，窗口 top-1 为 0.731，健康误报率为 0.112。
- Z24 Stage 02（相对参考状态最早记录到的偏离）：流方向 AUROC 为 0.810 ± 0.010，最佳静态基线为 0.683。
- FEMTO 轴承：因果退化趋势读出的单调性与趋势性平均得分为 0.453，这不是剩余寿命预测。
- 一块 ESP32-S3 以每个状态对 34.552 ms 复现预载特征的五通道数值路径，与 PC 结果接近。

!!! note "适用范围"
    未发表工作，见博士论文。“分布式”指提出的叶节点—网关划分，在集中采集数据和公开基准数据上评估，并非实时多板网络。Z24 Stage 02 包含沉降系统安装和可能的环境变化，因此不是确认的损伤检测。

提出的表示—读出架构（逻辑流程，并非实测的实时多板实现）

![MANIFOLD-FRAMEWORK](MANIFOLD-FRAMEWORK.png){width=100%}

三个验证系统上的状态轨迹（仅用于可视化）

![MANIFOLD-STATES](MANIFOLD-STATES.png){width=100%}

Z24 最早阶段的读出结果

![MANIFOLD-Z24](MANIFOLD-Z24.png){width=100%}

<div class="grid cards" markdown>

-   :material-file:{ .lg .middle } __博士论文 - 第 6 章（未发表）__

    ---

    **Cui, S.** Distributed Edge Intelligence for IoT-Based Wireless Structural Health Monitoring. Ph.D. Thesis, Nanyang Technological University.

</div>


## __IV 教研工作__

### ArduinoNode

设计并开发了一个基于 Arduino UNO R4 WIFI 的简单无线智能传感器网络节点，成功应用于香港理工大学举办的 APESS2025 暑期学校的结构健康监测教学中。

![](apess2025.jpg)

<div class="grid cards" markdown>

-   :simple-github:{ .lg .middle } __APESS2025_ARDUINONODE__

    ---

    SHM教学用WSN节点--基于Arduino UNO R4 WIFI

    [:octicons-arrow-right-24: <a href="https://github.com/Shuaiwen-Cui/APESS2025_ArduinoNode.git" target="_blank"> 代码 </a>](#)

    [:octicons-arrow-right-24: <a href="https://shuaiwen-cui.github.io/APESS2025_ArduinoNode/" target="_blank"> 在线文档 </a>](#)


</div>

<div class="grid cards" markdown>

-   :material-file:{ .lg .middle } __期刊论文 - IEEE Transactions on Education__

    ---

    **Cui, S.**, Fu, Y.*, Xia, Y., Zhang Q., & Li, S. (2026). A Class-Lab-Field Pedagogical Framework for Structural Health Monitoring using Ultra-Low-Cost Wireless IoT Prototypes. IEEE Transactions on Education. 

    [:octicons-arrow-right-24: <a href="https://doi.org/10.1109/TE.2026.3717462" target="_blank"> DOI </a>](#)

    [:octicons-arrow-right-24: <a href="https://mp.weixin.qq.com/s/d7ucOQSB-XSaVzWxegbzBg?scene=1&click_id=524330546" target="_blank"> 公众号推送 </a>](#)

</div>