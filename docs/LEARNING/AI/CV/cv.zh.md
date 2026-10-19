# 视觉、多模态与三维空间

!!! abstract "这一页解决什么"
    视觉任务从"是什么"（分类）扩展到"在哪里"（检测、分割）、"多远多大"（深度、三维）和"怎么用语言描述"（视觉语言）。这是机器人感知和很多工业 AI 岗位的核心。**P1 必修任务与指标，P2 按方向深入三维**。

## 知识地图

| 层次 | 主题 | 需要理解什么 |
|---|---|---|
| P1 | 分类、检测与分割 | 标签、框、像素掩码及对应指标；CNN / ViT、YOLO / DETR、提示式分割 |
| P1 | 视觉语言与视频 | 图文对齐、VLM、时序上下文 |
| P2 | 三维几何与表示 | 深度、点云、位姿、占据 / SDF、NeRF、3D 高斯泼溅 |
| P2 | 空间与交互 | 几何重建、语义、可供性与动作条件 |

## 核心机制

**分类到检测到分割**：分类输出一个标签；检测输出框加类别（YOLO 单阶段、快；Faster R-CNN 两阶段、准；DETR 用 Transformer 直接预测框集合、无需 NMS）；语义分割给每个像素分类，实例分割还区分同类的不同个体（Mask R-CNN、YOLO-seg）。

**评价指标**：IoU 是预测与真值区域的交并比；mAP 是按置信度排序得到 PR 曲线后积分、再对类别和 IoU 阈值平均。分割常用 mIoU。

**ViT**：把图像切成 16×16 的块当 token 送进 Transformer。数据少时不如 CNN（缺少局部性先验），数据多时更强；端侧部署 CNN 仍更友好。

**提示式分割（SAM）**：给一个点或框，分割出对应物体，不需要针对类别训练。适合标注和机器人场景中的新物体。

**深度与三维**：单目深度估计（如 Depth Anything）、双目、RGB-D；三维表示有点云、网格、体素、SDF；NeRF 用网络隐式表示场景、可渲染新视角；3D 高斯泼溅用大量带颜色的高斯椭球显式表示，渲染快得多。

**视觉语言**：CLIP 对齐图文空间，VLM 把视觉 token 接入语言模型（见基础模型的多模态专题）。

## 精选资源

<div class="grid cards" markdown>

-   :material-school:{ .lg .middle } __主课选读 · Stanford CS231n 2026__

    ---

    先掌握分类与迁移学习，再延伸到检测分割、自监督与 3D 视觉。

    [:octicons-arrow-right-24: 打开课表](https://cs231n.stanford.edu/schedule)

-   :material-target:{ .lg .middle } __实践 · Ultralytics 检测与分割__

    ---

    YOLO 的训练、验证、导出一条龙，最快跑通检测和实例分割。

    [:octicons-arrow-right-24: 打开文档](https://docs.ultralytics.com/tasks/segment/)

-   :material-cube-scan:{ .lg .middle } __三维选修 · Nerfstudio Splatfacto__

    ---

    体验 3D 高斯泼溅重建，测新视角质量与渲染时延。

    [:octicons-arrow-right-24: 打开文档](https://docs.nerf.studio/nerfology/methods/splat.html)

</div>

## 跟随实验

1. 先完成表示学习页的 AI-2；不要再建立一套重复的分类训练项目。
2. 选检测、分割或三维任务之一，画出输入输出与评估流程。
3. 选三维时按 Splatfacto 官方数据运行原例，留出新视角不参与训练。
4. 测新视角效果、渲染时间和失败区域，区分视觉逼真与几何准确。

**验收**：三维报告须说明相机、数据、算力和评估视角，不能仅凭好看的重建图宣称空间理解。机器人方向的完整感知项目见下方 P1。

## 自测

??? question "检测模型的 mAP@0.5 和 mAP@0.5:0.95 有什么区别？"
    前者只要求 IoU ≥ 0.5 就算对；后者在 0.5 到 0.95 的十个阈值上取平均，对定位精度要求严格得多。

??? question "为什么小数据集上 ViT 常不如 CNN？"
    CNN 自带局部性和平移不变性先验，少量数据就能学好；ViT 先验弱，需要大数据或强预训练才能学到这些规律。

??? question "NeRF 和 3D 高斯泼溅的主要区别？"
    NeRF 用 MLP 隐式表示，渲染要对每条光线采样多次网络，慢；3DGS 用显式高斯椭球加光栅化，实时渲染，但存储更大。

## 相关专题

<div class="grid cards" markdown>

-   __P1 · 6DoF 位姿估计（具身）__

    ---

    检测分割加深度，得到可抓取的物体位姿。

    [:octicons-arrow-right-24: 进入](../../EMBODIED-INTELLIGENCE/PROJECT-BASED-LEARNING/pose-6dof.md)

-   __多模态与行动接口（基础模型）__

    ---

    视觉语言模型与 VLA。

    [:octicons-arrow-right-24: 进入](../../FOUNDATION-MODELS/MULTIMODAL/multimodal.md)

-   __表示学习__

    ---

    AI-2 迁移学习实验。

    [:octicons-arrow-right-24: 进入](../REPRESENTATION/representation.md)

</div>
