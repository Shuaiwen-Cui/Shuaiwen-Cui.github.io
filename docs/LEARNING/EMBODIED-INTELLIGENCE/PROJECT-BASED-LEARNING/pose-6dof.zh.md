# 实战 · RGB-D 物体 6DoF 位姿估计

!!! abstract "一句话"
    给一帧 RGB-D 图像和物体的三维模型，算出物体在相机坐标系下的旋转和平移——机器人抓取之前必须先回答的问题。这个项目从公开基准走到边缘设备，最后用 iPhone 实拍收尾。

**代码：**[Pose6D-Edge](https://github.com/Shuaiwen-Cui/Pose6D-Edge)（本地 `E:\DEV\Pose6D-Edge`）

## 在学习地图中的位置

它是[感知与状态估计](../PERCEPTION/perception.md)里「操作：检测、6DoF 位姿与抓取」这条选修的第一个完整项目，也对应[视觉与三维空间](../../AI/CV/cv.md)的 P2 层次（深度、点云、位姿）。和 E0–E5 的区别：E 系列从控制与策略入手，本项目从**感知**入手，直接对齐工业机器人感知岗位（检测、6DoF 位姿、三维重建、抓取，以及模型量化部署）。

```
RGB ──► 实例分割（YOLO-seg，或用框提示 SAM）──► 物体掩码
                                                  │
深度 ─► 用相机内参反投影成点云 ◄──────────────────┘
         │
         ▼
  FPFH + RANSAC 全局配准 ─► 多尺度点到面 ICP（另试绕主轴翻转 180°）
         │
         ▼
  用渲染一致性（z-buffer 轮廓 + 深度）给候选打分 ─► 位姿 (R, t)
         │
         ▼
  视频：上一帧位姿做初值的 ICP 跟踪，分数低时用 SAM 重新初始化
```

## 先学什么（按项目进度排序）

| 顺序 | 知识点 | 为什么先学 | 一套资源 | 自测 |
|---|---|---|---|---|
| 1 | 针孔相机模型、内参、反投影 | 深度图变点云的唯一入口 | [Open3D RGBD 与点云教程](https://www.open3d.org/docs/release/tutorial/geometry/rgbd_image.html)；OpenCV 相机标定文档的模型部分 | 手推像素 (u, v, z) 到相机坐标 (x, y, z) |
| 2 | 刚体变换 SE(3)、旋转的几种表示 | 位姿就是一个 SE(3)；评估、翻转、跟踪全是变换的复合 | Modern Robotics 第 3 章 3.1–3.3（与具身主课一致） | 解释 `T_cam_obj` 与 `T_obj_cam` 的关系 |
| 3 | 点云配准：ICP（点到点 / 点到面）、FPFH + RANSAC | 本项目的几何核心；面试必问 ICP 的初值依赖 | [Open3D ICP](https://www.open3d.org/docs/release/tutorial/pipelines/icp_registration.html) 与 [全局配准](https://www.open3d.org/docs/release/tutorial/pipelines/global_registration.html)；Cyrill Stachniss 的 ICP 讲座 | 说清 ICP 为什么会陷入局部最优、对称物体为什么会翻转 |
| 4 | 位姿评估：ADD / ADD-S、BOP 指标 | 不会评估就不能声称结果 | [BOP 基准](https://bop.felk.cvut.cz/) 的挑战说明 | 解释对称物体为什么用 ADD-S |
| 5 | 实例分割微调与合成数据 | 从真值掩码换成学习到的掩码，才是完整系统 | Ultralytics Segment 文档；CS231n 检测与分割讲 | 说清为什么用 PBR 合成图训练、测试集不能进训练 |
| 6 | TensorRT FP16 / INT8 量化与测速 | 工业岗位常见的加分项，也是你的强项 | Ultralytics 导出文档；NVIDIA TensorRT INT8 校准说明 | 报告量化前后的时延与精度变化 |
| 7 | 学习型位姿方法的全貌 | 面试时要知道几何基线的边界在哪 | FoundationPose（CVPR 2024）、MegaPose、GDR-Net 论文的方法图 | 说清学习型方法在遮挡与对称上强在哪里 |

!!! tip "可迁移"
    系统辨识与状态估计的经验直接对应「配准 = 求解最优刚体变换」「跟踪 = 用上一时刻状态做先验」。TinyAuton 的量化与算子经验对应第 6 步。

## 十天计划

| 天 | 任务 | 验收 |
|---|---|---|
| 1 | 装环境，下载 LM-O，用真值掩码跑 `eval_pose.py`；读 `geometry.py`、`registration.py` | 复现 README 中真值掩码的结果 |
| 2–3 | 下载 LM 的 PBR 合成图，`prepare_yolo.py` 转格式，在本机显卡训练 YOLO-seg | 验证集掩码 mAP，训练曲线 |
| 4 | `eval_pose.py --masks yolo`，得到完整系统的指标 | 真值掩码与学习掩码两行对比 |
| 5–6 | Jetson 上导出 TensorRT FP16 / INT8，`bench_latency.py` 测速，再测量化后精度 | 时延与精度对比表 |
| 7–8 | iPhone 扫描物体建模，用 Stray Scanner 录 RGB-D，跑 `iphone_demo.py` | 一段带位姿叠加的视频与 GIF |
| 9–10 | 整理 README、推 GitHub、录一段讲解视频 | 仓库首页能看到结果表与 GIF |

## 已有结果（云端 CPU 实测）

| 设置 | ADD(-S) < 0.1d |
|---|---|
| LM-O 全部 1,445 个目标，真值掩码 | **0.643**（中位误差 7.6° / 8.6 mm） |
| 同 150 个目标：只按点拟合选假设 | 0.52 |
| 同 150 个目标：加 180° 翻转假设 + 渲染一致性打分 | **0.62** |
| 同 20 个目标：框提示 SAM 掩码（平均 IoU 0.75） | 0.55 |
| 跟踪一步（初值扰动 5° / 10 mm） | 0.85 |

失败主要来自两类：近似对称物体在单侧视角下的 180° 翻转，以及可见比例低于 30% 的重遮挡。这正是 FoundationPose 等学习型方法的优势所在；本项目是透明的几何基线加边缘部署链路，不声称达到最新水平。

## 记录要求

沿用[具身项目的统一记录要求](project-based-learning.md)：版本与依赖、数据来源、固定的评估集、成功率与失败样例、个人修改。离线指标、仿真或数据集结果、真机实拍分别报告，不互相替代。

[返回跟随式项目](project-based-learning.md) · [感知与状态估计](../PERCEPTION/perception.md)
