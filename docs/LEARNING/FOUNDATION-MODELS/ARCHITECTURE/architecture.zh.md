# 架构与表示

!!! abstract "这一页解决什么"
    大模型的内部是怎么组成的：文字怎样变成向量、注意力怎样让 token 互相"看"、一层层模块怎样堆叠成 GPT。看懂这一页，后面的训练、推理、多模态都有了共同的骨架。**必修**；前置：矩阵乘法、PyTorch 张量。

## 学到什么程度

| 层次 | 学什么 | 完成标准 |
|---|---|---|
| 扫盲 | token、嵌入、上下文、注意力、生成 | 画出输入到输出的过程，完成 F0 |
| 基础 | 残差、归一化、位置编码、训练循环 | 跟随 F4 实现小 GPT，解释每个张量的形状 |
| 进阶 | 长上下文、MoE、状态空间模型、线性注意力 | 针对一个问题设计比较实验 |

## 核心知识

**分词（tokenization）**：BPE 等算法把文本切成子词单元。词表大小决定嵌入矩阵大小和序列长度；中文、代码、数字的切分方式直接影响成本和能力（例如数位被切碎会让算术变难）。

**嵌入**：每个 token 查表得到一个向量；最后的输出层把向量映射回词表上的概率。很多模型让输入嵌入和输出层共享权重。

**自注意力**：每个 token 生成 Q、K、V；softmax(QKᵀ/√d)·V 让每个位置按内容相关性汇聚其他位置的信息。因果掩码保证生成时只看过去。多头注意力让不同头关注不同关系。

**Transformer 块**：注意力 + 前馈网络（FFN，参数量的大头），每个子层外包残差连接和归一化。现代 LLM 普遍用 Pre-Norm（RMSNorm）和 SwiGLU 激活，训练更稳。

**位置编码**：注意力本身不感知顺序。绝对位置编码给每个位置加一个向量；RoPE 把 Q、K 按位置旋转，使注意力只依赖相对距离，是当前主流。

**长上下文**：注意力是 O(n²)，KV cache 随长度线性增长。手段包括滑动窗口、稀疏注意力、FlashAttention（不改复杂度，减少显存读写）、RoPE 外推与插值。

**混合专家（MoE）**：把 FFN 换成多个专家，路由器每个 token 只激活其中几个。总参数大、每 token 计算少；难点是负载均衡和显存（所有专家都要放得下）。

**替代路线**：状态空间模型（Mamba）和线性注意力把历史压缩成固定大小的状态，推理每步 O(1)，长序列更省；代价是精确回忆长距离细节较弱，混合架构正在流行。

**多模态统一表示**：图像切成 patch、音频切成帧，经编码器变成与文字同维的 token，再一起送进 Transformer。

## 路线之争 · 自回归 vs. 扩散

自回归逐 token 生成，训练简单、与语言天然匹配，是 LLM 的主流；扩散模型并行去噪生成，在图像、视频和机器人动作上占优，也有人尝试扩散式语言模型以并行解码提速。分歧本质是"顺序生成的简单性"和"并行生成的效率与多峰表达"之间的取舍。

## 精选资源

<div class="grid cards" markdown>

-   :material-code-braces:{ .lg .middle } __主线 · Karpathy Zero to Hero__

    ---

    makemore 前两讲 → "Let's build GPT"；不懂梯度先补 micrograd。

    [:octicons-arrow-right-24: 打开课程](https://karpathy.ai/zero-to-hero.html)

-   :material-book-open-page-variant:{ .lg .middle } __文字替代 · LLMs-from-scratch__

    ---

    第 2–5 章逐章 Notebook，从分词到 GPT 预训练；和上一条二选一。

    [:octicons-arrow-right-24: 打开仓库](https://github.com/rasbt/LLMs-from-scratch)

-   :hugging:{ .lg .middle } __接口查阅 · HF LLM Course__

    ---

    第 2 章 tokenizer 与模型加载，遇到接口问题时查。

    [:octicons-arrow-right-24: 打开课程](https://huggingface.co/learn/llm-course/en/chapter1/1)

</div>

## 跟随实验

<div class="grid cards" markdown>

-   :material-hammer-wrench:{ .lg .middle } __F0 · token 与上下文__

    ---

    比较中英文分词，只改上下文观察回答变化。

    [:octicons-arrow-right-24: 查看步骤](../PROJECT-BASED-LEARNING/project-based-learning.md)

-   :material-hammer-wrench:{ .lg .middle } __F4 · 小型 GPT__

    ---

    跟写 GPT，改上下文长度或模型宽度中的一个因素。

    [:octicons-arrow-right-24: 查看步骤](../PROJECT-BASED-LEARNING/project-based-learning.md)

</div>

## 自测

??? question "注意力为什么要除以 √d？"
    Q·K 的方差随维度 d 增大，数值过大时 softmax 接近 one-hot、梯度消失；除以 √d 让方差保持在 1 左右。

??? question "因果掩码怎么实现？"
    在 softmax 之前把注意力分数矩阵的上三角（未来位置）设为负无穷，使其权重为 0。

??? question "一个 Transformer 块里参数主要在哪里？"
    在前馈网络：两个线性层，中间维度通常是隐藏维度的 4 倍左右（SwiGLU 约 8/3 倍、三个矩阵），参数约为注意力投影的两倍。

??? question "MoE 为什么"总参数大但推理便宜"？显存呢？"
    每个 token 只经过少数专家，计算量按激活参数算；但所有专家权重都要加载，显存按总参数算。

## 相关专题

<div class="grid cards" markdown>

-   __Transformer（AI）__

    ---

    张量骨架与数学细节。

    [:octicons-arrow-right-24: 进入](../../AI/TRANSFORMER/transformer.md)

-   __训练与适配__

    ---

    同一个架构怎样被训练和改造。

    [:octicons-arrow-right-24: 进入](../TRAINING/training.md)

</div>
