# 随机子空间识别（SSI）

本页对 [SHM roadmap](../shm.zh.md) 中的 **2.8 SSI** 做展开： **随机子空间识别（Stochastic Subspace Identification, SSI）** 是 **子空间方法** 中面向「仅输出、随机激励」的一支，也是结构健康监测中最重要、应用最广的仅输出（或输入–输出）模态识别方法之一。从输出时间序列出发，SSI 构造块 Hankel/Toeplitz 矩阵，通过 **投影** 与 **SVD** 得到 **可观性子空间**（及可控性子空间），恢复离散状态空间模型，并提取 **固有频率**、**阻尼比** 与 **振型**。两种主要实现形式为 **SSI-COV**（协方差驱动）：用输出协方差滞后构成块 Toeplitz 矩阵； **SSI-DATA** （数据驱动）：用原始数据块 Hankel 与斜投影。两者最终得到同类型的状态空间模型与模态参数；选择哪种主要取决于实现方式、数值性质与数据使用上的权衡。

---

## 动机与背景：为什么需要 SSI？

### 实际约束催生「仅输出」识别

- **环境激励不可测**：大型结构（桥梁、建筑、风机塔架）的激励来自风、交通、地脉动等，难以施加已知力或测量。传统模态分析依赖 **已知输入**（锤击、激振器）做频响或脉冲识别，在现场往往不可行。
- **运行状态需求**：关心的是结构在 **真实运行/服役** 下的动力特性（刚度、阻尼、振型），而非实验室条件下的结果。只有响应传感器（加速度、应变等）时，需要 **仅凭输出时间序列** 估计模态参数。
- **目标**：在 **不已知、不测量输入** 的前提下，从多通道输出 \(\mathbf{y}[k]\) 估计 **固有频率 \(f_r\)、阻尼比 \(\zeta_r\)、振型 \(\mathbf{\phi}_r\)**，即 **运行模态分析（OMA）**。

### 为什么选用状态空间模型？

- **统一框架**：离散 LTI 状态空间 \((\mathbf{A},\mathbf{C})\) 天然包含动力演化（\(\mathbf{A}\)）与观测映射（\(\mathbf{C}\)）；模态参数（极点、振型）由其特征值/特征向量唯一决定。
- **与随机激励相容**：若将未知环境激励与量测噪声建模为 **白噪声**，则状态方程 (1a)(1b) 恰好描述「白噪声驱动 + 仅观测输出」的设定；系统仍是 LTI，便于理论分析。
- **可辨识性**：在可观性成立时，**输出的二阶统计（多滞后协方差）** 由 \((\mathbf{A},\mathbf{C})\) 唯一决定（差一相似变换）。因此从协方差或从数据块的几何结构（投影）可以 **反推** 可观性子空间，再恢复 \(\mathbf{A},\mathbf{C}\)。这是 SSI 可行的根本原因。

### 小结

SSI 的动机可概括为：**在仅能观测响应、激励未知或未测的现实下，利用「输出协方差/数据块由可观性矩阵唯一决定」这一数学事实，从输出数据恢复状态空间与模态参数。** 下面给出支撑这一做法的核心数学原理。

---

## 核心数学原理：为什么仅输出就能做系统识别？

### 最核心的数学原理（一句话 + 两条根本原理）

如果只说一个核心，那就是：**可观性子空间（observability subspace）可以从输出数据中恢复出来。**

再展开一点，是两条根本原理：

- **原理 1：动态系统生成的数据具有低秩/低维子空间结构。**  
  由输出构成的块矩阵（协方差的块 Toeplitz、或数据的块 Hankel/投影）在理想情况下**秩为 \(n\)**（系统阶次），其行/列空间由可观性矩阵张成。因此用 SVD 即可从数据中「挖出」可观性子空间的一组基。
- **原理 2：状态空间模型具有移位不变性。**  
  可观性矩阵的块行满足「下一块行 = 当前块行 × \(\mathbf{A}\)」。一旦得到可观性子空间的一组基（即可观性矩阵的估计），即可通过该移位关系**唯一恢复** \(\mathbf{A}\) 和 \(\mathbf{C}\)。

这两条合起来，就足够支撑整个子空间识别：**先由数据得到可观性子空间（原理 1），再从移位关系得到状态空间参数（原理 2）**。下面分别给出更具体的表述与推导。

---

### 原理一：输出协方差由可观性与可控性唯一决定

在模型 (1) 中，设 \(\mathbf{w}[k],\mathbf{v}[k]\) 为零均值、宽平稳白噪声，且与初始状态不相关。则输出 \(\mathbf{y}[k]\) 为宽平稳过程，其 **滞后 \(\ell\) 的协方差** 为
\[
\mathbf{R}_\ell = \mathbb{E}[\mathbf{y}[k+\ell]\mathbf{y}[k]^T].
\]
可以证明（见下方「块 Toeplitz 分解」），存在与噪声协方差相关的 **可控性型矩阵** \(\mathcal{C}_i\)（块结构，与 \(\mathcal{O}_i\) 对称），使得 **块 Toeplitz 矩阵**（由 \(\mathbf{R}_1,\ldots,\mathbf{R}_{2i-1}\) 按 (4) 排成）满足
\[
\mathbf{T}_{1|i} = \mathcal{O}_i \, \mathcal{C}_i.
\]
即：**多滞后输出协方差排成的块 Toeplitz 矩阵，可分解为可观性矩阵 × 可控性型矩阵。** 因此：
- \(\mathbf{T}_{1|i}\) 的 **行空间**（左奇异向量张成的空间）与 \(\mathcal{O}_i\) 的 **列空间** 一致（差一个可逆变换）；
- 对 \(\mathbf{T}_{1|i}\) 做 **SVD**，取左奇异向量（或 \(\mathbf{U}\mathbf{\Sigma}^{1/2}\)）即得到 **可观性矩阵列空间的一组基**，从而得到 \(\mathcal{O}_i\) 的估计（差一个基的选取与相似变换）。这就是 **SSI-COV 中「对 T 做 SVD 就能得到可观性」** 的数学依据。

### 原理二：可观性矩阵的移位结构唯一决定 (A, C)

可观性矩阵 \(\mathcal{O}_i\) 的块行满足递推关系：
\[
\text{第 } \alpha+1 \text{ 块行} = \mathbf{C}\mathbf{A}^\alpha = (\text{第 } \alpha \text{ 块行})\,\mathbf{A}.
\]
因此若记：
- \(\mathcal{O}_{i,\mathrm{up}}\)：\(\mathcal{O}_i\) 去掉 **最后一块行**；
- \(\mathcal{O}_{i,\mathrm{down}}\)：\(\mathcal{O}_i\) 去掉 **第一块行**，

则有 **移位关系** \(\mathcal{O}_{i,\mathrm{down}} = \mathcal{O}_{i,\mathrm{up}} \, \mathbf{A}\)。故
\[
\mathbf{A} = \mathcal{O}_{i,\mathrm{up}}^\dagger \, \mathcal{O}_{i,\mathrm{down}}.
\]
而 **输出矩阵** \(\mathbf{C}\) 就是 \(\mathcal{O}_i\) 的 **第一块行**。因此：**一旦从数据得到可观性矩阵 \(\mathcal{O}_i\)（或其列空间的一组基），\(\mathbf{A}\) 和 \(\mathbf{C}\) 即可通过上述移位与第一块行唯一恢复。** 这是所有 SSI 变体（COV/DATA）在「恢复 A、C」一步共用的核心数学事实。

### 原理三：离散极点与连续模态参数的对应

离散状态矩阵 \(\mathbf{A}\) 的特征值 \(\mu_j\) 与连续时间极点 \(\lambda_j\) 在采样间隔 \(\Delta t\) 下满足 \(\mu_j = e^{\lambda_j \Delta t}\)（零阶保持或精确离散化）。因此 \(\lambda_j = \ln(\mu_j)/\Delta t\)。对共轭对 \(\lambda_r,\lambda_r^*\)：
- \(|\lambda_r|\) 为固有圆频率 \(\omega_r\)；
- \(-\mathrm{Re}(\lambda_r)/|\lambda_r|\) 为阻尼比 \(\zeta_r\)；
- \(\mathbf{C}\) 与 \(\mathbf{A}\) 对应特征向量的乘积给出传感器处振型。  
因此 **从 \((\mathbf{A},\mathbf{C})\) 到模态参数** 的步骤是确定性的，无需额外假设。

### 原理四：数据驱动形式（SSI-DATA）与投影的等价性

在 SSI-DATA 中，用 **块 Hankel 矩阵** 的「过去」与「未来」块做 **斜投影**（未来对过去的回归），可以证明：该投影的列空间与 **可观性子空间** 一致。直观上，「由过去可预测的那部分未来」恰好由系统的状态（可观性）所张成；对投影做 SVD 即得到可观性矩阵的基。因此 DATA 与 COV 在数学上 **同源**：一个从协方差（二阶统计）出发，一个从原始数据的几何（投影）出发，都得到同一个可观性子空间，进而用同一套移位关系得到 \(\mathbf{A},\mathbf{C}\)。

### 小结：SSI 可行的核心链条

1. **LTI + 白噪声** → 输出协方差序列有确定结构；
2. **块 Toeplitz = \(\mathcal{O}_i \mathcal{C}_i\)** → 从协方差或数据投影可恢复 \(\mathcal{O}_i\) 的列空间（SVD）；
3. **移位关系** → 从 \(\mathcal{O}_i\) 唯一恢复 \(\mathbf{A},\mathbf{C}\)；
4. **特征值/特征向量** → 从 \(\mathbf{A},\mathbf{C}\) 得到 \(f_r,\zeta_r,\mathbf{\phi}_r\)。

因此：**「可观性 + 移位结构 + 离散/连续极点对应」** 是 SSI 能进行仅输出系统识别的 **最核心的数学原理**。

---

## 子空间方法的数学本质

**「子空间方法」** 泛指一类系统识别算法：不从参数空间直接拟合，而是先从一个 **由数据构造的块结构矩阵**（Hankel、Toeplitz 或类似）中，用 **线性代数**（SVD、投影、QR 等）提取出 **可观性子空间**（或可控性子空间），再由此子空间 **唯一恢复** 状态空间参数 \((\mathbf{A},\mathbf{C})\)。SSI 是其中面向 **随机、仅输出** 的一支；与它对应的是 **确定性子空间方法**（输入已知）。二者共享同一套数学本质，区别仅在于 **数据从何而来、块矩阵如何构造、分解式中的「驱动侧」是什么**。

**直观地说**：子空间方法就是在 **找出张成「可观性子空间」的一组基**。我们不直接估计 \((\mathbf{A},\mathbf{C})\)，而是先利用数据矩阵的分解结构，用 SVD（或投影 + SVD）得到该子空间的一组基向量（例如 \(\mathbf{U}_n \mathbf{\Sigma}_n^{1/2}\) 的列，或投影的左奇异向量）。有了这组基，就等价于知道了可观性矩阵的列张成的空间；再利用可观性矩阵的 **移位结构**，从这组基唯一恢复 \(\mathbf{A}\) 和 \(\mathbf{C}\)。因此可以概括为两步：**先找子空间的基，再从基恢复状态空间参数**。

### 统一数学框架：块矩阵分解 + 子空间提取 + 移位恢复

1. **块矩阵的分解结构**  
   在 LTI 状态空间模型下，总可构造一个与数据相关的块矩阵 \(\mathbf{H}\)（行/列按时间或滞后分块），使其在理想情况下满足
   \[
   \mathbf{H} = \mathcal{O} \cdot \mathcal{Z}.
   \]
   - \(\mathcal{O}\)：**可观性矩阵**（或它的列张成的子空间），只依赖 \((\mathbf{A},\mathbf{C})\)，与输入/噪声的具体形式无关。
   - \(\mathcal{Z}\)：**「驱动」矩阵**，依赖输入或噪声如何驱动状态（确定性输入下的可控性块，或随机下的 \(\mathbf{G}\)/新息相关块）。
   因此 \(\mathbf{H}\) 的 **行空间**（或列空间，视分块方式）由 \(\mathcal{O}\) 的列张成；即 **可观性子空间** 嵌在数据矩阵的结构里。

2. **子空间提取**  
   对 \(\mathbf{H}\) 做 **SVD**（或对数据块做 **正交/斜投影** 再 SVD），取前 \(n\) 个左奇异向量（或投影的左奇异向量）张成的子空间，即得到 **可观性矩阵列空间的一组基**（差一可逆变换）。阶次 \(n\) 由奇异值截断或稳定图确定。

3. **移位恢复 (A, C)**  
   可观性矩阵 \(\mathcal{O}\) 的块行满足 \(\text{下一块行} = \text{当前块行} \times \mathbf{A}\)。故从估计出的 \(\mathcal{O}\) 中取「去掉首块行」与「去掉末块行」的两块，用 **移位关系** \(\mathcal{O}_{\mathrm{down}} = \mathcal{O}_{\mathrm{up}} \mathbf{A}\) 得到 \(\mathbf{A} = \mathcal{O}_{\mathrm{up}}^\dagger \mathcal{O}_{\mathrm{down}}\)；\(\mathbf{C}\) 为 \(\mathcal{O}\) 的 **第一块行**。这一步与输入/噪声是确定性还是随机 **无关**，是子空间方法共有的第二步。

**本质概括**：子空间方法 = **用数据矩阵的几何（子空间）代替直接参数估计**；数学上依赖 (1) 块矩阵 = 可观性 × 驱动矩阵，(2) SVD/投影提取可观性子空间，(3) 可观性的移位结构唯一决定 \((\mathbf{A},\mathbf{C})\)。

---

## 确定性子空间方法 vs 随机子空间方法

| 维度 | **确定性子空间方法** | **随机子空间方法**（如 SSI） |
|------|----------------------|------------------------------|
| **输入** | **已知、可测**（锤击、激振器、已知力信号） | **未知或未测**，建模为白噪声（或宽平稳噪声） |
| **数据** | 输入–输出对 \((\mathbf{u}[k],\mathbf{y}[k])\)，或 **脉冲响应** / 阶跃响应采样 | **仅输出** \(\mathbf{y}[k]\)，或输出的多滞后 **协方差** \(\mathbf{R}_\ell\) |
| **块矩阵** | 多由 **Markov 参数**（脉冲响应 \(h[\ell] = \mathbf{C}\mathbf{A}^{\ell-1}\mathbf{B}\)）排成 **块 Hankel**；或由输入–输出数据块构成 Hankel | **块 Toeplitz**（协方差 \(\mathbf{R}_1,\ldots,\mathbf{R}_{2i-1}\)）或 **块 Hankel**（输出数据）+ **未来对过去的投影** |
| **分解中的「驱动侧」** | **可控性矩阵** \(\mathcal{C}_i = [\mathbf{B},\,\mathbf{A}\mathbf{B},\,\ldots,\,\mathbf{A}^{i-1}\mathbf{B}]\)，由 **已知输入矩阵 \(\mathbf{B}\)** 与 \(\mathbf{A}\) 决定 | **可控性型矩阵**（如 \([\mathbf{G},\mathbf{A}\mathbf{G},\ldots]\)），由 **状态–输出互协方差 \(\mathbf{G}\)**、噪声协方差等决定；无 \(\mathbf{B}\)，仅有输出统计 |
| **子空间含义** | 同一可观性 \(\mathcal{O}_i\)；可控性由输入通道张成 | 同一可观性 \(\mathcal{O}_i\)；「驱动侧」由噪声/新息统计张成，不涉及显式输入 |
| **典型算法** | **ERA**（Eigensystem Realization Algorithm）、**Ho–Kalman**、**N4SID**（确定性部分）、**PO-MOESP** 等 | **SSI-COV**、**SSI-DATA**、**N4SID**（随机部分）、**CVA** 等 |
| **应用场景** | 实验室/可控激励：锤击、激振器、已知输入的模态试验 | 现场/环境激励：仅响应传感器、运行模态分析（OMA） |
| **假设** | 已知输入；线性；可得到脉冲响应或足够丰富的输入–输出数据 | 平稳性；激励可视为白噪声（或宽频）；仅输出即可辨识 \((\mathbf{A},\mathbf{C})\)（可观性 + 噪声统计结构） |

**联系**：两类方法都基于 **同一可观性矩阵与移位结构**，从「数据块矩阵 = 可观性 × 驱动矩阵」中通过 SVD/投影得到可观性子空间，再恢复 \(\mathbf{A},\mathbf{C}\) 与模态参数。**区别**：确定性方法用 **已知输入** 构造块矩阵，驱动侧是 **可控性（\(\mathbf{B}\)）**；随机方法用 **仅输出或输出协方差** 构造块矩阵，驱动侧是 **噪声相关的 \(\mathbf{G}\)/新息**，不出现 \(\mathbf{B}\)。SSI 属于 **随机** 这一支，专为「无已知输入、仅输出」的 OMA 设计。

---

## 概念

**核心思路** — 在 **环境激励**（未知或未测量）下，结构响应可视为由白噪声驱动的 **线性时不变（LTI）离散状态空间系统** 的输出。该系统的 **可观性矩阵** 张成的子空间可从仅输出数据中恢复：原因是在白噪声激励与 \((\mathbf{A},\mathbf{C})\) 可观的前提下，输出的二阶统计（多滞后协方差）或输出数据块的结构（块 Toeplitz、块 Hankel 的投影）由可观性矩阵唯一决定，因此用仅输出数据估计协方差或构造投影、再做 SVD，即可得到可观性子空间的一组基。**SSI** 的做法是：

1. **SSI-COV：** 估计多滞后输出的 **协方差矩阵**，排成 **块 Toeplitz** 矩阵，对其做 SVD 得到可观性矩阵的列空间；再由可观性矩阵恢复状态矩阵 \(\mathbf{A}\) 与输出矩阵 \(\mathbf{C}\)。
2. **SSI-DATA：** 用输出数据构造 **块 Hankel 矩阵**，划分为「过去」与「未来」块，计算未来行空间对过去行空间的 **斜投影**，对（加权形式的）该投影做 SVD 得到可观性矩阵；与 COV 相同方式恢复 \(\mathbf{A}\)、\(\mathbf{C}\)。

由离散状态矩阵 \(\mathbf{A}\) 的特征值得到 **离散极点**，再换算为 **连续时间极点**（固有频率 \(\omega_r\)、阻尼比 \(\zeta_r\)）。由 \(\mathbf{C}\) 与 \(\mathbf{A}\) 的特征向量得到传感器处的 **振型** 。无需已知输入力（仅输出运行模态分析）。

**状态空间模型（离散时间）** — 底层 LTI 模型为

\[\mathbf{x}[k+1] = \mathbf{A} \mathbf{x}[k] + \mathbf{w}[k], \tag{1a}\]

\[\mathbf{y}[k] = \mathbf{C} \mathbf{x}[k] + \mathbf{v}[k]. \tag{1b}\]

状态按 (1a) 演化，我们只能观测到输出 (1b)。SSI 的目标即从输出数据恢复 \(\mathbf{A}\)、\(\mathbf{C}\)，进而得到模态参数。

**式 (1) 符号含义：**

- \(\mathbf{x}[k] \in \mathbb{R}^n\)：状态向量；\(n\) 为 **模型阶次**（状态维数）。
- \(\mathbf{A} \in \mathbb{R}^{n \times n}\)：状态矩阵；其特征值为离散极点，决定模态频率与阻尼。
- \(\mathbf{C} \in \mathbb{R}^{p \times n}\)：输出矩阵；将状态映射到测量输出，用于求振型。
- \(\mathbf{y}[k] \in \mathbb{R}^p\)：输出向量（如 \(p\) 个传感器加速度）。
- \(\mathbf{w}[k], \mathbf{v}[k]\)：过程噪声与量测噪声（常假设为白噪声、零均值）；输入隐含在 \(\mathbf{w}\) 中。

**可观性** — **可观性矩阵**（\(\mathbf{C}\) 与 \(\mathbf{A}\) 幂次之块堆叠）为

\[\mathcal{O}_i = \begin{bmatrix} \mathbf{C} \\ \mathbf{C}\mathbf{A} \\ \vdots \\ \mathbf{C}\mathbf{A}^{i-1} \end{bmatrix}. \tag{2}\]

SSI 恢复的是 \(\mathcal{O}_i\) **列空间** 的一组基（差一个相似变换）。一旦得到 \(\mathcal{O}_i\)，\(\mathbf{A}\)、\(\mathbf{C}\) 即可推出：其第一块行为 \(\mathbf{C}\)；\(\mathbf{A}\) 由移位结构恢复（如 \(\mathcal{O}_{i,\mathrm{up}}^\dagger \mathcal{O}_{i,\mathrm{down}}\)，其中「up」去掉最后一块行，「down」去掉第一块行）。因此 SSI 的核心步骤是从数据得到 \(\mathcal{O}_i\)，两种实现形式的区别即在这一步的做法不同。

**SSI-COV 与 SSI-DATA** — **SSI-COV** 基于 **协方差估计**：用输出协方差 \(\mathbf{R}_\ell = \mathbb{E}[\mathbf{y}_{k+\ell}\mathbf{y}_k^T]\)（或样本估计）构成块 Toeplitz 矩阵，对该矩阵（或其因子）做一次 SVD 即得可观性。**SSI-DATA** 基于 **原始数据**：用 \(\mathbf{y}[k]\) 构造块 Hankel 矩阵，计算「未来」输出对「过去」输出的正交/斜投影，对该投影做 SVD 得到可观性。COV 先将数据压缩为协方差（数据量小、规模固定）；DATA 使用完整数据块，并可加权重（如 CVA）以改善数值与统计性质。两者终点一致：可观性 → \(\mathbf{A},\mathbf{C}\) → 模态参数。

---

## SSI-COV：协方差驱动算法

### 步骤 1：数据准备与协方差估计

**输入：** 多通道输出时间序列 \(\mathbf{y}[k] \in \mathbb{R}^p\)，\(k = 1,\ldots,N\)，采样间隔 \(\Delta t\)，\(f_s = 1/\Delta t\)。通常 \(p \geq 2\)（建议更多）个传感器，\(N\) 较大（如数万至百万）以得到稳定协方差估计。

**预处理：**

1.1 去直流：\(\mathbf{y}[k] \leftarrow \mathbf{y}[k] - \frac{1}{N}\sum_{j=1}^N \mathbf{y}[j]\)（按通道）。

1.2 （可选）去趋势和/或带通滤波至关心频段。

**协方差估计：**

1.3 选取 **块行数** \(i\)（如 20～80）。所需滞后数至少 \(2i-1\)。对滞后 \(\ell = 0,1,\ldots,2i-1\)，估计 **输出协方差**：

\[\mathbf{R}_\ell = \frac{1}{N-\ell} \sum_{k=1}^{N-\ell} \mathbf{y}[k+\ell] \mathbf{y}[k]^T. \tag{3}\]

其中 \(\mathbf{R}_0\) 为滞后 0 的（样本）输出协方差；\(\mathbf{R}_\ell\) 为 \(\mathbf{y}[k+\ell]\) 与 \(\mathbf{y}[k]\) 的互协方差。这些协方差是下一步构造 Toeplitz 矩阵的基础。

**数学依据**：在模型 (1) 且过程/量测噪声为白噪声、系统渐近稳定时，输出为宽平稳过程，\(\mathbf{R}_\ell\) 收敛到理论值 \(\mathbb{E}[\mathbf{y}[k+\ell]\mathbf{y}[k]^T]\)。该理论协方差由 \((\mathbf{A},\mathbf{C})\) 及噪声协方差唯一决定；因此用样本协方差 (3) 估计 \(\mathbf{R}_\ell\)，即可为后续「从协方差恢复可观性」提供一致估计（\(N \to \infty\) 时）。

**实践要点：**

- \(i\) 需足够大，使 \(\mathcal{O}_i\) 对目标阶次 \(n\) 列满秩（一般 \(i \cdot p \geq n\)，且 \(\mathbf{A}^{i-1}\) 已充分衰减）；常用 \(i = 20\)～\(50\)。
- 协方差数据量：\(N\) 需远大于最大滞后（\(N \gg 2i\)）以降低估计方差；土木结构（如 1～10 Hz 模态、100 Hz \(f_s\)）常用 2～10 分钟数据。

### 步骤 2：构造块 Toeplitz 矩阵

2.1 由协方差构成 **块 Toeplitz 矩阵** \(\mathbf{T}_{1|i} \in \mathbb{R}^{pi \times pi}\)：

\[\mathbf{T}_{1|i} = \begin{bmatrix} \mathbf{R}_1 & \mathbf{R}_2 & \cdots & \mathbf{R}_i \\ \mathbf{R}_2 & \mathbf{R}_3 & \cdots & \mathbf{R}_{i+1} \\ \vdots & \vdots & \ddots & \vdots \\ \mathbf{R}_i & \mathbf{R}_{i+1} & \cdots & \mathbf{R}_{2i-1} \end{bmatrix}. \tag{4}\]

（有的文献在首块列用 \(\mathbf{R}_0\)，下标会整体平移一位；关键是保持 Toeplitz 结构且可分解为可观性 × 可控性。）

2.2 **为什么对 T 做 SVD 就能得到可观性？—— 块 Toeplitz 的分解推导**

在状态空间模型 (1) 下，输出可写为 \(\mathbf{y}[k] = \mathbf{C}\mathbf{A}^{k-k_0}\mathbf{x}[k_0] +\)（噪声引起的项）。在平稳且白噪声假设下，理论协方差可写成（通过 **新息形式** 或 **Lyapunov/Stein 方程**）：
\[
\mathbf{R}_\ell = \mathbf{C} \mathbf{A}^{\ell-1} \mathbf{G}, \quad \ell \geq 1,
\]
其中 \(\mathbf{G} = \mathbb{E}[\mathbf{x}[k+1]\mathbf{y}[k]^T]\) 为状态–输出互协方差，由 \((\mathbf{A},\mathbf{C})\) 与噪声协方差唯一确定。将 \(\mathbf{R}_1,\ldots,\mathbf{R}_{2i-1}\) 按块 Toeplitz (4) 排成 \(\mathbf{T}_{1|i}\)，并利用 \(\mathbf{R}_\ell\) 的上述形式，逐块行写出：
- 第 1 块行：\([\mathbf{C}\mathbf{G},\, \mathbf{C}\mathbf{A}\mathbf{G},\, \ldots,\, \mathbf{C}\mathbf{A}^{i-1}\mathbf{G}] = \mathbf{C} \cdot [\mathbf{G},\, \mathbf{A}\mathbf{G},\, \ldots,\, \mathbf{A}^{i-1}\mathbf{G}]\);
- 第 2 块行：\([\mathbf{C}\mathbf{A}\mathbf{G},\, \ldots] = \mathbf{C}\mathbf{A} \cdot [\mathbf{G},\, \mathbf{A}\mathbf{G},\, \ldots,\, \mathbf{A}^{i-1}\mathbf{G}]\);
- ……
- 第 \(i\) 块行：\(\mathbf{C}\mathbf{A}^{i-1} \cdot [\mathbf{G},\, \mathbf{A}\mathbf{G},\, \ldots,\, \mathbf{A}^{i-1}\mathbf{G}]\).

因此 **按行分块** 可见 \(\mathbf{T}_{1|i}\) 的每一块行都是「\(\mathcal{O}_i\) 的对应块行 × 同一块矩阵」。把 \(\mathcal{O}_i\) 定义为 (2) 的 \(i\) 块行堆叠，把 \(\mathcal{C}_i\) 定义为 \([\mathbf{G},\, \mathbf{A}\mathbf{G},\, \ldots,\, \mathbf{A}^{i-1}\mathbf{G}]\) 及其在 Toeplitz 中按列重复的块结构，即有
\[
\mathbf{T}_{1|i} = \mathcal{O}_i \, \mathcal{C}_i.
\]
故 \(\mathbf{T}_{1|i}\) 的 **行空间** 由 \(\mathcal{O}_i\) 的列张成（当 \(\mathcal{C}_i\) 行满秩时即为 \(\mathcal{O}_i\) 的列空间）。对 \(\mathbf{T}_{1|i}\) 做 SVD 时，**左奇异向量** 即张成该行空间，因此前 \(n\) 个左奇异向量（或 \(\mathbf{U}_n \mathbf{\Sigma}_n^{1/2}\)）给出 **可观性矩阵 \(\mathcal{O}_i\) 列空间的一组基**。这就是「对 T 做 SVD 得到可观性」的完整数学依据。

### 步骤 3：SVD 与可观性

3.1 对块 Toeplitz 做 SVD：

\[\mathbf{T}_{1|i} = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^T. \tag{5}\]

3.2 **截断到阶次 \(n\)：** 选取 \(n\)（如 2× 物理模态数，或由稳定图确定）。保留前 \(n\) 个奇异值与向量：\(\mathbf{U}_n \in \mathbb{R}^{pi \times n}\)，\(\mathbf{\Sigma}_n \in \mathbb{R}^{n \times n}\)，\(\mathbf{V}_n \in \mathbb{R}^{pi \times n}\)。

**数学依据**：理想情况下 \(\mathbf{T}_{1|i} = \mathcal{O}_i \mathcal{C}_i\)，且 \(\mathcal{O}_i\) 列秩、\(\mathcal{C}_i\) 行秩均为 \(n\)（系统阶次），故 \(\mathbf{T}_{1|i}\) 的秩为 \(n\)，仅有 \(n\) 个非零奇异值。实际中噪声与估计误差会使小奇异值非零；截断到 \(n\) 即假定「信号子空间」维数为 \(n\)，对应 \(n\) 阶状态空间模型。

3.3 **可观性矩阵：** 令

\[\mathcal{O}_i = \mathbf{U}_n \mathbf{\Sigma}_n^{1/2}. \tag{6}\]

则 \(\mathcal{O}_i\) 有 \(i\) 个块行，每块 \(p \times n\)。其 **第一块行** 即为输出矩阵 \(\mathbf{C}\)；其余块行的移位关系将在步骤 4 中用于恢复 \(\mathbf{A}\)。

**为何取 \(\mathbf{U}_n \mathbf{\Sigma}_n^{1/2}\)？** SVD 给出 \(\mathbf{T}_{1|i} = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^T\)，行空间由 \(\mathbf{U}\) 的列张成。因子分解 \(\mathbf{T}_{1|i} = (\mathbf{U}_n \mathbf{\Sigma}_n^{1/2})(\mathbf{\Sigma}_n^{1/2} \mathbf{V}_n^T)\) 使左边因子与 \(\mathcal{O}_i\) 张成同一列空间且列数为 \(n\)，便于按块行取「第一块行 = \(\mathbf{C}\)」及「移位求 \(\mathbf{A}\)」；取 \(\mathbf{U}_n \mathbf{\Sigma}_n^{1/2}\) 而非仅 \(\mathbf{U}_n\) 是为了在相似变换意义下与真实 \(\mathcal{O}_i\) 的尺度一致，使恢复的 \(\mathbf{A}\) 特征值不变。

### 步骤 4：恢复 A 与 C

4.1 定义 **移位可观性** 矩阵：

- \(\mathcal{O}_{i,\mathrm{up}}\)：\(\mathcal{O}_i\) **去掉最后一块行**（\((i-1)p \times n\)）。
- \(\mathcal{O}_{i,\mathrm{down}}\)：\(\mathcal{O}_i\) **去掉第一块行**（\((i-1)p \times n\)）。

**移位关系的数学依据**：按可观性矩阵定义 (2)，\(\mathcal{O}_i\) 的第 \(\alpha\) 块行为 \(\mathbf{C}\mathbf{A}^{\alpha-1}\)（\(\alpha = 1,\ldots,i\)）。故 \(\mathcal{O}_{i,\mathrm{up}}\) 的块行为 \(\mathbf{C},\mathbf{C}\mathbf{A},\ldots,\mathbf{C}\mathbf{A}^{i-2}\)，\(\mathcal{O}_{i,\mathrm{down}}\) 的块行为 \(\mathbf{C}\mathbf{A},\mathbf{C}\mathbf{A}^2,\ldots,\mathbf{C}\mathbf{A}^{i-1}\)。因此 \(\mathcal{O}_{i,\mathrm{down}}\) 的每一块行都等于 \(\mathcal{O}_{i,\mathrm{up}}\) 的「下一块行 × \(\mathbf{A}\)」，即
\[
\mathcal{O}_{i,\mathrm{down}} = \mathcal{O}_{i,\mathrm{up}} \, \mathbf{A}.
\]
当 \(\mathcal{O}_{i,\mathrm{up}}\) 列满秩（可观性成立且 \(i\) 足够大）时，\(\mathbf{A}\) 可唯一恢复为
\[\mathbf{A} = \mathcal{O}_{i,\mathrm{up}}^\dagger \, \mathcal{O}_{i,\mathrm{down}}. \tag{7}\]

其中 \(\dagger\) 表示 Moore–Penrose 伪逆。

4.2 **输出矩阵：** \(\mathbf{C}\) 为 \(\mathcal{O}_i\) 的 **第一块行**：

\[\mathbf{C} = \mathcal{O}_i(1:p, 1:n). \tag{8}\]

### 步骤 5：模态参数提取（COV 与 DATA 共用）

5.1 离散特征值：求 \(\mathbf{A}\) 的特征值 \(\mu_j\)，\(j=1,\ldots,n\)。剔除单位圆外或明显非物理的（如实负根）。

5.2 **连续时间极点：** 离散状态矩阵 \(\mathbf{A}\) 的特征值 \(\mu_j\) 对应连续时间极点 \(\lambda_j\)。对采样间隔 \(\Delta t\) 的离散化，标准关系为 \(\mu_j = e^{\lambda_j \Delta t}\)，反解即得

\[\lambda_j = \frac{\ln(\mu_j)}{\Delta t} \quad \text{（复数取主值）}. \tag{9}\]

**数学依据**：连续时间 LTI 系统 \(\dot{\mathbf{x}} = \mathbf{A}_c \mathbf{x} + \ldots\) 在采样间隔 \(\Delta t\) 下的精确离散化（或零阶保持）为 \(\mathbf{x}[k+1] = \mathbf{A} \mathbf{x}[k] + \ldots\)，其中 \(\mathbf{A} = e^{\mathbf{A}_c \Delta t}\)。故 \(\mathbf{A}\) 的特征值 \(\mu_j\) 与 \(\mathbf{A}_c\) 的特征值 \(\lambda_j\) 满足 \(\mu_j = e^{\lambda_j \Delta t}\)。对单自由度振荡 \(s^2 + 2\zeta\omega_n s + \omega_n^2 = 0\)，极点 \(\lambda = -\zeta\omega_n \pm j\omega_n\sqrt{1-\zeta^2}\)，故 \(|\lambda| = \omega_n\)，\(\mathrm{Re}(\lambda) = -\zeta\omega_n\)，即 (10a)(10b)。

对对应一个物理模态的 **共轭对** \(\lambda_r, \lambda_r^*\)，固有圆频率（rad/s）与阻尼比为

\[\omega_r = |\lambda_r|, \tag{10a}\]

\[\zeta_r = -\frac{\mathrm{Re}(\lambda_r)}{|\lambda_r|}. \tag{10b}\]

固有频率（Hz）为 \(f_r = \omega_r/(2\pi)\)。

5.3 **振型：** 设 \(\mathbf{\psi}_r \in \mathbb{C}^n\) 为 \(\mathbf{A}\) 对应 \(\mu_r\) 的右特征向量。传感器处 **振型** 为

\[\mathbf{\phi}_r = \mathbf{C} \mathbf{\psi}_r. \tag{11}\]

做归一化（如单位范数或最大分量置 1）。对实结构可按需取实部或模。

**数学依据**：状态空间模型中，输出 \(\mathbf{y} = \mathbf{C}\mathbf{x}\)。若状态在模态 \(r\) 下按 \(\mathbf{x} \propto \mathbf{\psi}_r\) 演化（\(\mathbf{A}\mathbf{\psi}_r = \mu_r \mathbf{\psi}_r\)），则各传感器测到的幅值比由 \(\mathbf{C}\mathbf{\psi}_r\) 给出，即该模态在 \(p\) 个测点上的 **振型向量** \(\mathbf{\phi}_r\)。因此 (11) 是「状态特征向量 → 输出空间振型」的标准映射。

---

## SSI-DATA：数据驱动算法

SSI-DATA 直接对原始输出数据操作而非协方差：先构造块 Hankel 矩阵，再通过投影步得到与 SSI-COV 相同的可观性子空间，之后恢复 \(\mathbf{A}\)、\(\mathbf{C}\) 及模态参数的步骤与 COV 一致。

### 步骤 1：数据准备与块 Hankel 矩阵

**输入：** 与 COV 相同：\(\mathbf{y}[k] \in \mathbb{R}^p\)，\(k = 1,\ldots,N\)，\(\Delta t\)。预处理：去直流；可选去趋势与带通滤波。

**块 Hankel 矩阵：**

1.1 选取 **块行数** \(i\)（过去）与 \(i\)（未来），以及 **列数** \(j\)。通常 \(j = N - 2i + 1\)（用满可用数据）。要求 \(j \geq 2i\) 且足够大以保证统计稳定。

1.2 构造 **块 Hankel 矩阵** \(\mathbf{Y}_{0|2i-1} \in \mathbb{R}^{2pi \times j}\)：

\[\mathbf{Y}_{0|2i-1} = \frac{1}{\sqrt{j}} \begin{bmatrix} \mathbf{y}[1] & \mathbf{y}[2] & \cdots & \mathbf{y}[j] \\ \mathbf{y}[2] & \mathbf{y}[3] & \cdots & \mathbf{y}[j+1] \\ \vdots & \vdots & \ddots & \vdots \\ \mathbf{y}[2i] & \mathbf{y}[2i+1] & \cdots & \mathbf{y}[2i+j-1] \end{bmatrix}. \tag{12}\]

因子 \(1/\sqrt{j}\) 用于归一化。拆分为：

- **过去：** \(\mathbf{Y}_{0|i-1}\) — 前 \(i\) 块行（\(pi \times j\)）。
- **未来：** \(\mathbf{Y}_{i|2i-1}\) — 后 \(i\) 块行（\(pi \times j\)）。

### 步骤 2：投影

2.1 「未来」行空间在「过去」行空间上的 **斜投影** 为

\[\mathbf{O}_i = \mathbf{Y}_{i|2i-1} \, \mathbf{Y}_{0|i-1}^T \bigl( \mathbf{Y}_{0|i-1} \mathbf{Y}_{0|i-1}^T \bigr)^{-1} \mathbf{Y}_{0|i-1}. \tag{13}\]

该投影刻画了由过去可预测的那部分未来，其列空间与可观性子空间相关。

**数学依据（为何投影给出可观性？）** 在状态空间模型 (1) 下，「过去」输出块 \(\mathbf{Y}_{0|i-1}\) 由过去状态与噪声生成，「未来」输出块 \(\mathbf{Y}_{i|2i-1}\) 由未来状态与噪声生成。在平稳性与白噪声假设下，**未来对过去的条件期望**（即由过去线性预测未来的最佳估计）落在由 **可观性矩阵 \(\mathcal{O}_i\) 张成的子空间** 内：因为「可预测部分」对应状态的确定性演化 \(\mathbf{x} \mapsto \mathbf{A}\mathbf{x}\)，而观测为 \(\mathbf{y} = \mathbf{C}\mathbf{x}\)，故可预测的未来输出由 \(\mathcal{O}_i\) 的列张成。斜投影 (13) 正是「未来在过去的（线性）回归」，其列空间因此与 \(\mathcal{O}_i\) 的列空间一致（差一可逆变换）。对投影做 SVD 即得到该子空间的一组基，从而得到可观性矩阵的估计。实际实现中一般不显式形成该矩阵，而是对过去/未来堆叠矩阵做 **RQ** 或 **LQ** 分解，或对加权后的组合做 **SVD**，以因子形式表示投影。

2.2 **加权投影**（如 CVA — 典型变量分析）：定义过去与未来的权重矩阵

\[\mathbf{W}_p = \bigl( \mathbf{Y}_{0|i-1} \mathbf{Y}_{0|i-1}^T \bigr)^{-1/2}, \tag{14a}\]

\[\mathbf{W}_f = \bigl( \mathbf{Y}_{i|2i-1} \mathbf{Y}_{i|2i-1}^T \bigr)^{-1/2}. \tag{14b}\]

对 **加权投影** 做 SVD：

\[\mathbf{W}_f \, \mathbf{Y}_{i|2i-1} \, \mathbf{Y}_{0|i-1}^T \, \mathbf{W}_p^T = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^T. \tag{15}\]

（也可采用其他加权；无加权即 \(\mathbf{W}_p = \mathbf{I}\)，\(\mathbf{W}_f = \mathbf{I}\)。）

2.3 可观性：截断到阶次 \(n\)：\(\mathbf{U}_n, \mathbf{\Sigma}_n\)。则

\[\mathcal{O}_i = \mathbf{W}_f^{-1} \mathbf{U}_n \mathbf{\Sigma}_n^{1/2}. \tag{16}\]

即由加权投影的左奇异向量（及奇异值）得到可观性矩阵，再用 \(\mathbf{W}_f^{-1}\)「去权」。

### 步骤 3：恢复 A 与 C

3.1 与 SSI-COV 相同：由 \(\mathcal{O}_i\) 得到 \(\mathcal{O}_{i,\mathrm{up}}\)、\(\mathcal{O}_{i,\mathrm{down}}\)。则

\[\mathbf{A} = \mathcal{O}_{i,\mathrm{up}}^\dagger \, \mathcal{O}_{i,\mathrm{down}}, \tag{17a}\]

\[\mathbf{C} = \mathcal{O}_i(1:p, 1:n). \tag{17b}\]

### 步骤 4：模态参数提取

4.1 与 COV 步骤 5 相同：\(\mathbf{A}\) 的特征值 → 连续极点 (9)、(10a)–(10b) → \(f_r, \zeta_r\)；振型由 (11) 得到。

---

## 完整操作步骤（逐步）

**输入：** 多通道输出 \(\mathbf{y}[k]\)，\(k=1,\ldots,N\)；采样率 \(f_s\)，\(\Delta t = 1/f_s\)。

**通用预处理：**

1. 去直流（每通道减均值）。
2. （可选）带通滤波至关心频段。
3. 选取块行数 \(i\)（如 20～50）与模型阶次 \(n\)（如 2× 模态数；可用稳定图细化）。

**SSI-COV：**

4. 用 (3) 估计 \(\ell = 0,\ldots,2i-1\) 的协方差 \(\mathbf{R}_\ell\)。
5. 构造块 Toeplitz \(\mathbf{T}_{1|i}\) (4)；SVD (5)；截断到 \(n\)；由 (6) 得到 \(\mathcal{O}_i\)。
6. 按 (7)(8) 恢复 \(\mathbf{A}\)、\(\mathbf{C}\)。
7. \(\mathbf{A}\) 的特征值 → 极点 (9)、(10a)–(10b) → \(f_r, \zeta_r\)；振型 (11)。

**SSI-DATA：**

4. 构造块 Hankel \(\mathbf{Y}_{0|2i-1}\) (12)；划分过去/未来。
5. 形成（加权）投影并做 SVD (15)；截断到 \(n\)；由 (16) 得到 \(\mathcal{O}_i\)。
6. 按 (17a)(17b) 恢复 \(\mathbf{A}\)、\(\mathbf{C}\)。
7. 与 COV 相同的模态提取（特征值 → 极点 → \(f_r, \zeta_r\)；振型 (11)）。

**稳定图（两种共用）：** 对若干阶次 \(n = n_{\min}, \ldots, n_{\max}\) 重复识别，绘制极点（如频率–阶次）。**稳定** 极点（频率与阻尼随 \(n\) 变化很小）对应物理模态；漂移或散点多为虚假，仅保留稳定极点（如相邻阶次间 \(\Delta f/f < 1\%\)，\(\Delta\zeta/\zeta < 5\%\)）。

**输出：** 模态参数 — 固有频率 \(f_r\)，阻尼比 \(\zeta_r\)，振型 \(\mathbf{\phi}_r\)。

---

## 适用情形与局限

适用：

- **仅输出**（环境振动）或输入–输出数据；无需测量激励。
- 需要 **稳健、高精度** 的模态识别；SSI 是土木与机械 OMA 的常用标准。
- 需区分 **多阶模态**（含 **密模态**）；参数化状态空间与稳定图有助于分离。
- 需要 **阻尼**（SSI 直接给出，而 FDD 需 EFDD 补充）。
- 具备 **长时、多通道** 数据（SSI 受益于大 \(N\) 与 \(p\)）。

局限：

- **阶次选择：** 模型阶次 \(n\) 需人为选取；过小漏模态，过大引入虚假极点。应对多个 \(n\) 做 **稳定图**，只保留稳定极点。
- 计算与内存：块矩阵规模大（\(pi \times pi\) 或 \(2pi \times j\)）；SVD 与分解成本高；在资源受限的边缘设备上较难直接运行，常需简化或上云。
- 激励：假设持续激励（宽频或至少覆盖关心模态）；极窄带或非平稳激励会带来偏差。
- 噪声：噪声大会增加协方差/投影估计的方差；长数据与多通道可提高稳健性。

---

## 工程实践要点与可复现参数

| 方面 | 建议与典型取值 |
|------|----------------|
| 块行数 \(i\) | 足够大使 \(\mathcal{O}_i\) 对阶次 \(n\) 列满秩；典型 \(i = 20\)～\(50\)；\(i \cdot p \geq n\)。 |
| 模型阶次 \(n\) | 用稳定图：试 \(n = 8, 10, 12, \ldots\) 至约 2× 预期模态数；保留随阶次稳定的极点（如 \(\Delta f/f < 1\%\)，\(\Delta\zeta/\zeta < 5\%\)）。 |
| 数据长度 \(N\) | \(N \gg 2i\)；土木（如 1～10 Hz）在 100 Hz 下常用 2～10 分钟；低阻尼或低信噪比时更长。 |
| COV 与 DATA | COV：先压缩为协方差，一次 Toeplitz + 一次 SVD；DATA：完整数据 Hankel、投影 + SVD，可用 CVA 加权改善条件；结果类型相同。 |
| 数值稳定性 | 避免奇异或病态矩阵（如 DATA 中 \(\mathbf{Y}_{0|i-1}\mathbf{Y}_{0|i-1}^T\)）；用正则化或伪逆；剔除过小奇异值。 |
| 验证 | 与 FDD/EFDD 或 ERA 对比；用识别模型重算协方差或拟合；检查振型一致性。 |

---

## 边缘与在线计算

是否适用 — 对低功耗边缘 **较不适用**：大矩阵（Hankel/Toeplitz、SVD）、多阶次稳定图与长数据缓冲导致内存与算力需求高。更适合在服务器/云端或较强网关上运行；边缘可做轻量筛查（如 PP 或低阶 ARX），再上传片段做 SSI。

实现策略：

1. 缩小规模：限制 \(i\)、\(n\)；减少通道或协方差滞后数（COV）；限制 \(j\)（DATA）。
2. 边缘 COV：若在线（滑动窗或递推）计算协方差，可用较小 \(i\)、固定 \(n\) 在边缘跑 SSI-COV，完整稳定图离线。
3. 边缘 DATA：块 Hankel 与投影占内存大；可上传原始或降采样数据，在云端做 SSI。
4. 分级：边缘做快速检查（如 PP、简单相关）；可疑或定期片段上传；云端做完整 SSI 与稳定图。

挑战：内存（Hankel/Toeplitz）、CPU（SVD）、阶次选择（稳定图）以及长记录的数据量。

---

## 与其他方法比较

SSI vs FDD/EFDD：

- **SSI：** 时域、参数化、状态空间；在一套框架内给出频率、阻尼、振型；需阶次选择与稳定图；对密模态更有利。
- **FDD：** 频域、非参数；无阻尼（EFDD 补充）；无需阶次；实现简单但对阻尼与密模态精度较低。

SSI vs ERA / NExT-ERA：

- **SSI：** 由 **协方差（COV）** 或 **原始数据（DATA）** 构造 Hankel/Toeplitz，经投影 + SVD；稳健性强，为 OMA 常用标准；矩阵规模大。
- **ERA：** 由 **Markov 参数**（脉冲响应或 NExT 相关）构造 Hankel；一次 SVD、形式紧凑；NExT-ERA 通过相关实现仅输出。仅输出场景下 SSI 一般更稳健、灵活；当已有脉冲/相关时 ERA 更轻量。

SSI-COV vs SSI-DATA：

- **COV：** 一次压缩为协方差；一个 Toeplitz + SVD；若 \(N\) 很大可节省内存；协方差估计方差可能较大。
- **DATA：** 使用完整数据块与投影；可加 CVA 等加权以改善数值/统计性质；更灵活，软件中常首选；内存占用更高。

何时选择：

- 当仅输出 OMA 需要 **精度、阻尼与密模态分离** 的最佳折中时，优先 **SSI**；实现简单、内存固定可选 **SSI-COV**；在资源允许且追求稳健与质量时用 **SSI-DATA**（含 CVA）。
