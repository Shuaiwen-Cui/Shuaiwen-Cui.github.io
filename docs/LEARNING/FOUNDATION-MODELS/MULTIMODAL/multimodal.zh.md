# 多模态与行动接口

!!! abstract "这一页解决什么"
    模型怎样"看图、听声音"，以及怎样从"会描述"走到"会动手"：视觉语言模型（VLM）如何接上机器人的动作输出（VLA）。这是基础模型与具身智能的交汇点。**导论必修，VLA 适配选修**。

## 学到什么程度

| 层次 | 学什么 | 完成标准 |
|---|---|---|
| 扫盲 | 图文对齐、VLM、VLA 的输入输出 | 说清 VLM 和 VLA 的区别 |
| 基础 | 运行一个 VLM，测试空间与计数能力 | 完成 F3，记录至少三类错误 |
| 进阶 | 动作表示、机器人数据适配、闭环评估 | 完成具身 E4 的一条路线 |

## 核心知识

**对比式图文对齐（CLIP）**：图像编码器和文本编码器把配对的图和文映射到同一空间，拉近匹配对、推远不匹配对。结果是零样本分类和图文检索，也是很多 VLM 的视觉编码器。

**视觉语言模型（VLM）的结构**：视觉编码器（ViT）把图像变成 patch 向量，投影层把它们映射到语言模型的嵌入空间，当作"图像 token"和文字一起输入 LLM。

**VLM 的局限**：擅长描述和常识问答，但精确计数、细粒度空间关系、度量（距离、尺寸）和三维推理常出错。用于机器人前必须单独验证空间输出。

**从 VLM 到 VLA**：在 VLM 上加动作输出，用机器人数据（图像 + 指令 + 动作）训练。动作可以离散成 token 和文字一起生成（如 RT-2、OpenVLA），也可以用连续生成头（扩散或流匹配，如 π0）输出平滑的动作块。

**动作表示**：关节角或末端位姿、绝对或增量、单步或动作块。表示方式决定数据能否跨机器人复用，以及控制是否平滑。

**跨本体与数据**：不同机器人的关节数、相机位置、控制频率都不同。Open X-Embodiment 等数据集试图跨本体训练，但适配到你的机器人仍需本体数据。

**两种接入机器人的方式**：一是 LLM 只选择技能和参数，由已有控制器执行（可解释、易验证）；二是 VLA 直接输出动作（泛化潜力大、难调试）。

## 路线之争 · 离散动作 token vs. 连续动作生成

离散派把动作当成词表里的 token，能直接复用语言模型的训练和推理机制；连续派认为离散化损失精度，高频精细操作需要连续、平滑的动作块。目前高性能 VLA 多转向连续生成头，但离散方案实现最简单。

## 精选资源

<div class="grid cards" markdown>

-   :hugging:{ .lg .middle } __VLM 实践 · Image-text-to-text__

    ---

    Hugging Face 官方的图文问答任务指南，F3 直接使用。

    [:octicons-arrow-right-24: 打开指南](https://huggingface.co/docs/transformers/main/tasks/image_text_to_text)

-   :material-robot:{ .lg .middle } __VLA 实践 · SmolVLA__

    ---

    LeRobot 生态的小型 VLA，从数据到微调与评估的完整流程。

    [:octicons-arrow-right-24: 打开教程](https://huggingface.co/docs/lerobot/smolvla)

</div>

## 跟随实验

<div class="grid cards" markdown>

-   :material-hammer-wrench:{ .lg .middle } __F3 · VLM 能力与限制__

    ---

    10 张固定图片，测物体、数量、空间关系，记录错误类型。

    [:octicons-arrow-right-24: 查看步骤](../PROJECT-BASED-LEARNING/project-based-learning.md)

-   :material-hammer-wrench:{ .lg .middle } __E4 · 语言与动作__

    ---

    技能调用或 VLA 适配二选一，闭环验证。

    [:octicons-arrow-right-24: 查看步骤](../../EMBODIED-INTELLIGENCE/PROJECT-BASED-LEARNING/project-based-learning.md)

</div>

## 自测

??? question "CLIP 和 VLM 的区别？"
    CLIP 只把图和文映射到同一空间用于匹配，不能生成文字；VLM 把图像特征送进语言模型，可以对图像进行开放式问答。

??? question "为什么 VLM 说"杯子在盘子左边"不能直接用于抓取？"
    VLM 的空间描述是定性的、可能出错，且不给出度量位姿；抓取需要毫米级的三维位置与姿态，要由几何感知给出或验证。

??? question "VLA 的动作为什么常预测一个动作块而不是单步？"
    大模型推理慢，单步预测跟不上控制频率；动作块让一次推理覆盖多个控制周期，动作也更连贯。

??? question "把 LLM 用作技能选择器有什么好处？"
    动作由经过验证的控制器执行，安全可控；语言模型只负责理解指令和组合技能，失败容易定位。

## 相关专题

<div class="grid cards" markdown>

-   __学习型策略（具身）__

    ---

    VLA 在策略学习谱系中的位置。

    [:octicons-arrow-right-24: 进入](../../EMBODIED-INTELLIGENCE/POLICY-LEARNING/policy-learning.md)

-   __世界模型__

    ---

    预测环境变化，而不是直接输出动作。

    [:octicons-arrow-right-24: 进入](../WORLD-MODELS/world-models.md)

</div>
