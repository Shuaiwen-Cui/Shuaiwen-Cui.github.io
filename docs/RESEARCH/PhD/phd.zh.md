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

![](PLATFORM-ARCH.png){width=100%}

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

TinySHM 是运行在 MCU 上的中间件，位于板级驱动和 SHM 应用之间，覆盖采集、信号处理、系统识别、损伤筛查和设备端学习，不依赖 Linux 或 Python。

![](TINYSHM-STACK.png){width=100%}

!!! note "九个模块"
    - 平台适配：TinyToolbox
    - 通用计算：TinyMath / TinyDSP / TinyAI
    - SHM 专用：TinyMeasurement / TinySysID / TinyDamage
    - 任务编排：TinyBench / TinyOrch

设备端训练（TinyAI）

![](MOD-AI.png){width=100%}

![](TINYAI-PROPAGATION.png){width=100%}

节点端系统识别（TinySysID）

![](MOD-SYSID.png){width=100%}

<div class="grid cards" markdown>

-   :simple-github:{ .lg .middle } __TinySHM🎯🏆__

    ---

    当前支持平台：

    - ESP32

    <!-- ROLLBACK (github link hidden): [:octicons-arrow-right-24: <a href="https://github.com/Shuaiwen-Cui/TinySHM.git" target="_blank"> 代码 </a>](#) -->

    [:octicons-arrow-right-24: <a href="http://www.cuishuaiwen.com:8300/" target="_blank"> 在线文档（个人服务器 托管） </a>](#)

    <!-- ROLLBACK (github link hidden): [:octicons-arrow-right-24: <a href="https://shuaiwen-cui.github.io/TinySHM/" target="_blank"> 在线文档（Github 托管） </a>](#) -->

</div>

## __III 面向结构健康监测的机器学习与人工智能算法__

算法从单节点决策，推进到多节点间规则驱动的协同，再到从数据中学习的协作。

### __3.1 单节点智能__

#### 基于边缘智能与数字孪生的智能自适应触发传感用于能源高效的无线结构健康监测

🏷️ SHM分类：**测量**/系统识别/损伤检测/损伤定位/损伤量化

- 触发传感
- 闭环反馈控制
- 边缘智能
- 贝叶斯优化

![](SATM.jpg)

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

#### 微控制器节点上的板端模态识别与学习

🏷️ SHM分类：**测量**/**系统识别**/损伤检测/损伤定位/损伤量化

- 运行模态分析
- 设备端训练
- TinySHM

![](ONBOARD-MODAL.jpg){width=100%}

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

- 网关与节点协作
- 通信高效的系统识别

![](LCSI-PATH.png){width=100%}

**- 流形动力学**

🏷️ SHM分类：测量/**系统识别**/**损伤检测**/**损伤定位**/**损伤量化**

- 仅健康数据学习
- 学习得到的诊断状态
- 叶节点—网关表示

![](MANIFOLD-CONCEPT.jpg){width=100%}

三个验证系统上的状态轨迹（仅用于可视化）

![](MANIFOLD-STATES.png){width=100%}


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

## __V 实地部署与工程应用__

### __印度尼西亚 Flores__

2026 年 11 月，两台 LiftNode 部署于 Flores 的 RSUD Aeramo 医院，用于支持余震的持续远程监测。

<div class="grid cards deployment-gallery" markdown>

-   [![RSUD Aeramo 医院现场的部署团队](deployment/liftnode-deployment-team.webp){ loading=lazy width=1600 height=1200 }](deployment/liftnode-deployment-team.webp)

    **部署团队合照**

-   [![楼梯上下两处的 LiftNode 节点](deployment/liftnode-staircase-two-nodes.webp){ loading=lazy width=1080 height=1920 }](deployment/liftnode-staircase-two-nodes.webp)

    **楼梯间双节点布置**

-   [![柱脚处的 LiftNode 外壳与天线](deployment/liftnode-column.webp){ .deployment-detail loading=lazy width=1280 height=960 }](deployment/liftnode-column.webp)

    **柱脚测点**

-   [![楼梯平台栏杆立柱旁的 LiftNode](deployment/liftnode-stair-landing.webp){ .deployment-detail loading=lazy width=1280 height=720 }](deployment/liftnode-stair-landing.webp)

    **楼梯平台测点**

</div>

<figure markdown="span">
  [![Flores 楼梯双测点候选模态振型动画](deployment/flores-staircase-mode-shapes.gif){ loading=lazy width=1340 height=420 style="width:100%;height:auto;" }](deployment/flores-staircase-mode-shapes.gif)
  <figcaption>楼梯双测点模态动画 · 候选振型</figcaption>
</figure>

<p markdown="span" style="text-align:center;">
  [:material-newspaper-variant-outline: NTU CEE 官方报道 · 2026-11-15（英文）](https://www.ntu.edu.sg/cee/news/ntu-cee-researchers-conduct-field-deployment-in-flores--indonesia){ .md-button }
</p>

#### __2026 年 11 月 18 日 M5.8 地震附近时刻的实测记录__

11 月 18 日，部署在医院的主节点 **000** 于 **12:42:07.502（Flores／新加坡时间，UTC+8）** 记录到一段三轴振动，幅值随时间明显衰减。[BMKG 官方通报](https://www.bmkg.go.id/gempabumi/20261118114412)报告同一区域发生 **M5.8、深度 10 km** 的地震，换算到相同时区的发震时刻为 **12:41:32–33**；[GFZ GEOFON 地震目录](https://geofon.gfz.de/eqinfo/)也列出了约 04:41:32 UTC 的这次区域地震。独立台网时间与现场响应支持两者的候选关联，后续仍需附近参考台站波形进一步核对。

| 独立地震台网信息 | LiftNode 现场观测 |
| :-- | :-- |
| BMKG 发震时刻为 12:41:32–33，UTC+8 | 节点于 12:42:07.502 触发，约晚 35 秒 |
| 初报震中为 8.17°S、120.54°E；更新通报为 8.21°S、120.55°E | 测点位于 Flores 的 RSUD Aeramo 医院，记录来自主节点 000 |
| 区域发生 M5.8 地震，深度 10 km | 三轴记录在保存窗口内呈明显的振幅衰减 |
| BMKG 与 GFZ GEOFON 独立报告该事件 | ADXL367 预缓存与 ADXL355 采集数据保存在同一 SD 文件中，随后远程取回 |

<figure markdown="span" class="field-event-figure">
  [![2026年11月18日现场三轴融合时域波形与仅基于ADXL355的PSD](deployment/flores-20261118-waveform-psd.png){ loading=lazy width=800 height=600 style="width:100%;max-width:800px;height:auto;background:#fff;" }](deployment/flores-20261118-waveform-psd.svg)
  <figcaption>主节点 000 的现场记录。左侧为包含 ADXL367 前段和 ADXL355 主体的融合时域，虚线与点线标记触发、停止时刻；右侧为仅基于 ADXL355 的 Welch PSD，采用 dB 标度。三轴主要谱峰均在约 2.69 Hz，该结果表示响应谱峰，尚不是经独立辨识确认的结构模态。</figcaption>
</figure>

**这次记录验证了什么。** 现场节点完成了自主触发、触发前缓存、板载 AI 推理、SD 本地保存，以及真实振动记录的后续远程取回。完整文件为 **478,248 字节**，包含 **3,168 个 ADXL367 样本与 3,762 个 ADXL355 样本**；文件 SHA-256、两路传感器 CRC 及预缓存 CRC 均通过验证。节点实验模型将初始窗口标记为 *earthquake*，输出分数 **0.921**；独立云端模型输出 *unknown*。两项结果均作为参考保留，不用于确认事件归属。

<details markdown="1">
<summary>数据处理与解释边界</summary>

- 首个保留样本为 **12:42:06.535**，ADXL355 数据覆盖 **12:42:07.505 至 12:42:37.503**。文件仅包含触发前约 **0.97 秒** 的数据。30 秒记录长度由设定的最短事件时长控制，不等于地震持续时长；发震至触发的间隔也不是已拾取的地震波传播时间。
- 时域采用平台中参考 [Fu 等人方法](https://doi.org/10.3390/s18124480)实现的融合，输出为 100 Hz、40 Hz 低通，并去除记录末段中位数基线。轴向按 X–X、Y–Y、Z–Z 对应；重叠段检查提示轴向不匹配，因此当前为**强制同轴合并**，安装方向尚未经独立校验。
- PSD **只使用 ADXL355**，按记录时间戳重采样，采用 Hann 窗与 50% 重叠。4 个启动样本保留在原始文件中，在融合主体绘图与 PSD 中排除。ADXL355 有两处时间间隔超过 12 ms，最大为 **16.916 ms**；文件完整性验证不代表 ADC 零丢样。
- 本例来自 **000**，不声称已获得同一时刻的 001 记录或完成参考台站波形匹配。
- BMKG 初报为 11:41:32 WIB，更新通报为 11:41:33 WIB；加一小时即为上文的 Flores／新加坡时间。台网参数后续可能更新。

</details>

[下载 PNG 图片](deployment/flores-20261118-waveform-psd.png) · [SVG 矢量图](deployment/flores-20261118-waveform-psd.svg) · [原始波形 TXT](deployment/000_20261118T124207502_T.txt)

<details markdown="1">
<summary>原始文件 SHA-256</summary>

`8470de6a621cec303de7170e725bfea6cd122aa68f3d58a76fddab3a883da402`

</details>
