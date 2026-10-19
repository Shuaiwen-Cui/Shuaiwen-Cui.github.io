# Hands-on · RGB-D 6DoF Object Pose Estimation

!!! abstract "In one sentence"
    Given an RGB-D frame and a 3D model of an object, recover the object's rotation and translation in the camera frame — the question a robot must answer before it can grasp. This project goes from a public benchmark to an edge device and ends with a real iPhone capture.

**Code:** [Pose6D-Edge](https://github.com/Shuaiwen-Cui/Pose6D-Edge)

## Where it sits in the map

It is the first complete project of the "manipulation: detection, 6DoF pose and grasping" elective in [Perception & State Estimation](../PERCEPTION/perception.md), and the P2 level (depth, point clouds, pose) of [Vision & 3D Space](../../AI/CV/cv.md). Unlike E0–E5, which start from control and policies, this project starts from **perception** and lines up directly with industrial robot perception roles (detection, 6DoF pose, 3D reconstruction, grasping, plus model quantisation for deployment).

```
RGB ──► instance segmentation (YOLO-seg, or SAM from a box) ──► object mask
                                                                   │
Depth ─► back-projection with camera intrinsics ◄──────────────────┘
          │
          ▼
   FPFH + RANSAC global registration ─► multi-scale point-to-plane ICP (+ 180° flips)
          │
          ▼
   rank hypotheses by render consistency (z-buffer silhouette + depth) ─► pose (R, t)
          │
          ▼
   video: ICP tracking from the previous pose, re-initialised from SAM when the score drops
```

## What to learn first (in project order)

| # | Topic | Why now | One resource | Self-check |
|---|---|---|---|---|
| 1 | Pinhole camera, intrinsics, back-projection | The only way from a depth image to a point cloud | [Open3D RGBD tutorial](https://www.open3d.org/docs/release/tutorial/geometry/rgbd_image.html); the model section of the OpenCV calibration docs | Derive (u, v, z) → (x, y, z) by hand |
| 2 | Rigid transforms SE(3), rotation representations | A pose is an SE(3); evaluation, flips and tracking are all compositions | Modern Robotics ch. 3.1–3.3 | Explain `T_cam_obj` vs `T_obj_cam` |
| 3 | Registration: ICP (point-to-point / point-to-plane), FPFH + RANSAC | The geometric core; ICP's dependence on initialisation is a classic interview question | [Open3D ICP](https://www.open3d.org/docs/release/tutorial/pipelines/icp_registration.html) and [global registration](https://www.open3d.org/docs/release/tutorial/pipelines/global_registration.html); Cyrill Stachniss' ICP lecture | Why ICP falls into local minima; why symmetric objects flip |
| 4 | Pose metrics: ADD / ADD-S, BOP | No evaluation, no claim | [BOP benchmark](https://bop.felk.cvut.cz/) challenge description | Why symmetric objects use ADD-S |
| 5 | Fine-tuning instance segmentation, synthetic data | Replacing GT masks with learned ones makes it a full system | Ultralytics Segment docs; CS231n detection & segmentation lecture | Why train on PBR renders and never on the test split |
| 6 | TensorRT FP16 / INT8 and latency measurement | A common plus in industry roles and an existing strength | Ultralytics export docs; NVIDIA TensorRT INT8 calibration guide | Report latency and accuracy before and after quantisation |
| 7 | Landscape of learned pose estimators | Know where the geometric baseline stops | FoundationPose (CVPR 2024), MegaPose, GDR-Net method figures | Where learned methods win under occlusion and symmetry |

!!! tip "What transfers"
    System identification and state estimation map directly onto "registration = solving for the best rigid transform" and "tracking = using the previous state as a prior". TinyAuton's quantisation and kernel work maps onto step 6.

## Ten-day plan

| Day | Task | Done when |
|---|---|---|
| 1 | Set up, download LM-O, run `eval_pose.py` with GT masks; read `geometry.py` and `registration.py` | GT-mask result in the README reproduced |
| 2–3 | Download LM PBR renders, convert with `prepare_yolo.py`, train YOLO-seg on the desktop GPU | Validation mask mAP and training curves |
| 4 | `eval_pose.py --masks yolo` for the full system | GT-mask vs learned-mask rows side by side |
| 5–6 | Export TensorRT FP16 / INT8 on the Jetson, run `bench_latency.py`, re-evaluate accuracy | Latency and accuracy table |
| 7–8 | Scan an object with the iPhone, record RGB-D with Stray Scanner, run `iphone_demo.py` | A video and GIF with the pose overlay |
| 9–10 | Polish the README, push to GitHub, record a short walkthrough | Results table and GIF on the repo front page |

## Results so far (measured on a cloud CPU)

| Setting | ADD(-S) < 0.1d |
|---|---|
| All 1,445 LM-O targets, GT masks | **0.643** (median error 7.6° / 8.6 mm) |
| Same 150 targets, hypotheses ranked by point fit only | 0.52 |
| Same 150 targets, + 180° flip hypotheses and render-consistency ranking | **0.62** |
| Same 20 targets, box-prompted SAM masks (mean IoU 0.75) | 0.55 |
| One tracking step from a 5° / 10 mm perturbation | 0.85 |

Failures are mostly 180° flips of near-symmetric objects seen from one side, and heavily occluded instances (visible fraction below 30 %). That is exactly where learned methods such as FoundationPose are strong; this project is a transparent geometric baseline with an edge-deployment path, not a state-of-the-art claim.

## Record keeping

Follow the [common record rules of the embodied projects](project-based-learning.md): versions, data sources, a fixed evaluation set, success rates with failure examples, and your own changes. Report offline metrics, dataset results and real captures separately.

[Back to guided projects](project-based-learning.md) · [Perception & State Estimation](../PERCEPTION/perception.md)
