# 随机子空间识别（SSI）

本页对 [SHM roadmap](../shm.zh.md) 中的 **2.8 SSI** 做展开： **随机子空间识别（Stochastic Subspace Identification, SSI）** 是结构健康监测中最重要、应用最广的仅输出（或输入–输出）模态识别方法之一。从输出时间序列出发，SSI 构造块 Hankel/Toeplitz 矩阵，通过 **投影** 与 **SVD** 得到 **可观性子空间**（及可控性子空间），恢复离散状态空间模型，并提取 **固有频率**、**阻尼比** 与 **振型**。两种主要实现形式为 **SSI-COV**（协方差驱动）：用输出协方差滞后构成块 Toeplitz 矩阵； **SSI-DATA** （数据驱动）：用原始数据块 Hankel 与斜投影。两者最终得到同类型的状态空间模型与模态参数；选择哪种主要取决于实现方式、数值性质与数据使用上的权衡。

---

## 概念

**核心思路** — 在 **环境激励**（未知或未测量）下，结构响应可视为由白噪声驱动的 **线性时不变（LTI）离散状态空间系统** 的输出。该系统的 **可观性矩阵** 张成的子空间可从仅输出数据中恢复。**SSI** 的做法是：

1. **SSI-COV：** 估计多滞后输出的 **协方差矩阵**，排成 **块 Toeplitz** 矩阵，对其做 SVD 得到可观性矩阵的列空间；再由可观性矩阵恢复状态矩阵 \(\mathbf{A}\) 与输出矩阵 \(\mathbf{C}\)。
2. **SSI-DATA：** 用输出数据构造 **块 Hankel 矩阵**，划分为「过去」与「未来」块，计算未来行空间对过去行空间的 **斜投影**，对（加权形式的）该投影做 SVD 得到可观性矩阵；与 COV 相同方式恢复 \(\mathbf{A}\)、\(\mathbf{C}\)。

由离散状态矩阵 \(\mathbf{A}\) 的特征值得到 **离散极点**，再换算为 **连续时间极点**（固有频率 \(\omega_r\)、阻尼比 \(\zeta_r\)）。由 \(\mathbf{C}\) 与 \(\mathbf{A}\) 的特征向量得到传感器处的 **振型** 。无需已知输入力（仅输出运行模态分析）。

**状态空间模型（离散时间）** — 底层 LTI 模型为

\[\mathbf{x}[k+1] = \mathbf{A} \mathbf{x}[k] + \mathbf{w}[k], \qquad \mathbf{y}[k] = \mathbf{C} \mathbf{x}[k] + \mathbf{v}[k]. \tag{1}\]

式 (1) 符号含义：

- \(\mathbf{x}[k] \in \mathbb{R}^n\)：状态向量；\(n\) 为 **模型阶次**（状态维数）。
- \(\mathbf{A} \in \mathbb{R}^{n \times n}\)：状态矩阵；其特征值为离散极点，决定模态频率与阻尼。
- \(\mathbf{C} \in \mathbb{R}^{p \times n}\)：输出矩阵；将状态映射到测量输出，用于求振型。
- \(\mathbf{y}[k] \in \mathbb{R}^p\)：输出向量（如 \(p\) 个传感器加速度）。
- \(\mathbf{w}[k], \mathbf{v}[k]\)：过程噪声与量测噪声（常假设为白噪声、零均值）；输入隐含在 \(\mathbf{w}\) 中。

**可观性** — **可观性矩阵**（块行）为

\[\mathcal{O}_i = \begin{bmatrix} \mathbf{C} \\ \mathbf{C}\mathbf{A} \\ \vdots \\ \mathbf{C}\mathbf{A}^{i-1} \end{bmatrix}. \tag{2}\]

SSI 恢复的是 \(\mathcal{O}_i\) **列空间**的一组基（差一个相似变换）。由 \(\mathcal{O}_i\) 可得到 \(\mathbf{A}\) 和 \(\mathbf{C}\)：第一块行为 \(\mathbf{C}\)；\(\mathbf{A}\) 由移位结构恢复（如 \(\mathcal{O}_{i,\mathrm{up}}^\dagger \mathcal{O}_{i,\mathrm{down}}\)，其中「up」去掉最后一块行，「down」去掉第一块行）。

**SSI-COV 与 SSI-DATA** — **SSI-COV** 基于 **协方差估计**：用输出协方差 \(\mathbf{R}_\ell = \mathbb{E}[\mathbf{y}_{k+\ell}\mathbf{y}_k^T]\)（或样本估计）构成块 Toeplitz 矩阵，对该矩阵（或其因子）做一次 SVD 即得可观性。**SSI-DATA** 基于 **原始数据**：用 \(\mathbf{y}[k]\) 构造块 Hankel 矩阵，计算「未来」输出对「过去」输出的正交/斜投影，对该投影做 SVD 得到可观性。COV 先将数据压缩为协方差（数据量小、规模固定）；DATA 使用完整数据块，并可加权重（如 CVA）以改善数值与统计性质。两者终点一致：可观性 → \(\mathbf{A},\mathbf{C}\) → 模态参数。

---

## SSI-COV：协方差驱动算法

### 步骤 1：数据准备与协方差估计

输入：多通道输出时间序列 \(\mathbf{y}[k] \in \mathbb{R}^p\)，\(k = 1,\ldots,N\)，采样间隔 \(\Delta t\)，\(f_s = 1/\Delta t\)。通常 \(p \geq 2\)（建议更多）个传感器，\(N\) 较大（如数万至百万）以得到稳定协方差估计。

预处理：

1.1 去直流：\(\mathbf{y}[k] \leftarrow \mathbf{y}[k] - \frac{1}{N}\sum_{j=1}^N \mathbf{y}[j]\)（按通道）。

1.2 （可选）去趋势和/或带通滤波至关心频段。

**协方差估计：**

1.3 选取 **块行数** \(i\)（如 20～80）。所需滞后数至少 \(2i-1\)。对滞后 \(\ell = 0,1,\ldots,2i-1\)，估计 **输出协方差**：

\[\mathbf{R}_\ell = \frac{1}{N-\ell} \sum_{k=1}^{N-\ell} \mathbf{y}[k+\ell] \mathbf{y}[k]^T. \tag{3}\]

即 \(\mathbf{R}_0\) 为滞后 0 的（样本）输出协方差；\(\mathbf{R}_\ell\) 为 \(\mathbf{y}[k+\ell]\) 与 \(\mathbf{y}[k]\) 的互协方差。

实践要点：

- \(i\) 需足够大，使 \(\mathcal{O}_i\) 对目标阶次 \(n\) 列满秩（一般 \(i \cdot p \geq n\)，且 \(\mathbf{A}^{i-1}\) 已充分衰减）；常用 \(i = 20\)～\(50\)。
- 协方差数据量：\(N\) 需远大于最大滞后（\(N \gg 2i\)）以降低估计方差；土木结构（如 1～10 Hz 模态、100 Hz \(f_s\)）常用 2～10 分钟数据。

### 步骤 2：构造块 Toeplitz 矩阵

2.1 由协方差构成 **块 Toeplitz 矩阵** \(\mathbf{T}_{1|i} \in \mathbb{R}^{pi \times pi}\)：

\[\mathbf{T}_{1|i} = \begin{bmatrix} \mathbf{R}_1 & \mathbf{R}_2 & \cdots & \mathbf{R}_i \\ \mathbf{R}_2 & \mathbf{R}_3 & \cdots & \mathbf{R}_{i+1} \\ \vdots & \vdots & \ddots & \vdots \\ \mathbf{R}_i & \mathbf{R}_{i+1} & \cdots & \mathbf{R}_{2i-1} \end{bmatrix}. \tag{4}\]

（有的文献在首块列用 \(\mathbf{R}_0\)，下标会整体平移一位；关键是保持 Toeplitz 结构且可分解为可观性 × 可控性。）

2.2 在状态空间模型 (1) 下，\(\mathbf{T}_{1|i} = \mathcal{O}_i \, \mathcal{C}_i\)，其中 \(\mathcal{O}_i\) 为可观性矩阵 (2)，\(\mathcal{C}_i\) 为同结构的可控性矩阵。故 \(\mathbf{T}_{1|i}\) 的 **列空间** 与 \(\mathcal{O}_i\) 的列空间相同。

### 步骤 3：SVD 与可观性

3.1 对块 Toeplitz 做 SVD：

\[\mathbf{T}_{1|i} = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^T. \tag{5}\]

3.2 **截断到阶次 \(n\)：** 选取 \(n\)（如 2× 物理模态数，或由稳定图确定）。保留前 \(n\) 个奇异值与向量：\(\mathbf{U}_n \in \mathbb{R}^{pi \times n}\)，\(\mathbf{\Sigma}_n \in \mathbb{R}^{n \times n}\)，\(\mathbf{V}_n \in \mathbb{R}^{pi \times n}\)。

3.3 可观性矩阵：令

\[\mathcal{O}_i = \mathbf{U}_n \mathbf{\Sigma}_n^{1/2}. \tag{6}\]

则 \(\mathcal{O}_i\) 有 \(i\) 个块行，每块 \(p \times n\)；**第一块行** 即为输出矩阵 \(\mathbf{C}\)。

### 步骤 4：恢复 A 与 C

4.1 定义 **移位可观性** 矩阵：

- \(\mathcal{O}_{i,\mathrm{up}}\)：\(\mathcal{O}_i\) **去掉最后一块行**（\((i-1)p \times n\)）。
- \(\mathcal{O}_{i,\mathrm{down}}\)：\(\mathcal{O}_i\) **去掉第一块行**（\((i-1)p \times n\)）。

由状态空间结构有 \(\mathcal{O}_{i,\mathrm{down}} = \mathcal{O}_{i-1,\mathrm{up}} \mathbf{A}\)（\(\mathcal{O}_{i-1,\mathrm{up}}\) 为 \(\mathcal{O}_i\) 的前 \(i-1\) 块行）。故

\[\mathbf{A} = \mathcal{O}_{i,\mathrm{up}}^\dagger \, \mathcal{O}_{i,\mathrm{down}}, \tag{7}\]

其中 \(\dagger\) 表示 Moore–Penrose 伪逆。

4.2 输出矩阵：\(\mathbf{C}\) 为 \(\mathcal{O}_i\) 的 **第一块行**：

\[\mathbf{C} = \mathcal{O}_i(1:p, 1:n). \tag{8}\]

### 步骤 5：模态参数提取（COV 与 DATA 共用）

5.1 离散特征值：求 \(\mathbf{A}\) 的特征值 \(\mu_j\)，\(j=1,\ldots,n\)。剔除单位圆外或明显非物理的（如实负根）。

5.2 **连续时间极点：** 对每个特征值 \(\mu_j\)，

\[\lambda_j = \frac{\ln(\mu_j)}{\Delta t} \quad \text{（取主值）}. \tag{9}\]

对对应一个物理模态的 **共轭对** \(\lambda_r, \lambda_r^*\)：

\[\omega_r = |\lambda_r|, \qquad \zeta_r = -\frac{\mathrm{Re}(\lambda_r)}{|\lambda_r|}. \tag{10}\]

固有频率 \(f_r = \omega_r/(2\pi)\)。

5.3 振型：设 \(\mathbf{\psi}_r \in \mathbb{C}^n\) 为 \(\mathbf{A}\) 对应 \(\mu_r\) 的右特征向量。传感器处 **振型** 为

\[\mathbf{\phi}_r = \mathbf{C} \mathbf{\psi}_r. \tag{11}\]

做归一化（如单位范数或最大分量置 1）。对实结构可按需取实部或模。

---

## SSI-DATA：数据驱动算法

### 步骤 1：数据准备与块 Hankel 矩阵

输入：与 COV 相同：\(\mathbf{y}[k] \in \mathbb{R}^p\)，\(k = 1,\ldots,N\)，\(\Delta t\)。预处理：去直流；可选去趋势与带通滤波。

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

实际实现中一般不显式形成该矩阵，而是对过去/未来堆叠矩阵做 **RQ** 或 **LQ** 分解，或对加权后的组合做 **SVD**，以因子形式表示投影。

2.2 **加权投影（如 CVA — 典型变量分析）：** 定义

\[\mathbf{W}_p = \bigl( \mathbf{Y}_{0|i-1} \mathbf{Y}_{0|i-1}^T \bigr)^{-1/2}, \qquad \mathbf{W}_f = \bigl( \mathbf{Y}_{i|2i-1} \mathbf{Y}_{i|2i-1}^T \bigr)^{-1/2}. \tag{14}\]

对 **加权投影** 做 SVD：

\[\mathbf{W}_f \, \mathbf{Y}_{i|2i-1} \, \mathbf{Y}_{0|i-1}^T \, \mathbf{W}_p^T = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^T. \tag{15}\]

（也可采用其他加权；无加权即 \(\mathbf{W}_p = \mathbf{I}\)，\(\mathbf{W}_f = \mathbf{I}\)。）

2.3 可观性：截断到阶次 \(n\)：\(\mathbf{U}_n, \mathbf{\Sigma}_n\)。则

\[\mathcal{O}_i = \mathbf{W}_f^{-1} \mathbf{U}_n \mathbf{\Sigma}_n^{1/2}. \tag{16}\]

即由加权投影的左奇异向量（及奇异值）得到可观性矩阵，再用 \(\mathbf{W}_f^{-1}\)「去权」。

### 步骤 3：恢复 A 与 C

3.1 与 SSI-COV 相同：由 \(\mathcal{O}_i\) 得到 \(\mathcal{O}_{i,\mathrm{up}}\)、\(\mathcal{O}_{i,\mathrm{down}}\)；则

\[\mathbf{A} = \mathcal{O}_{i,\mathrm{up}}^\dagger \, \mathcal{O}_{i,\mathrm{down}}, \qquad \mathbf{C} = \mathcal{O}_i(1:p, 1:n). \tag{17}\]

### 步骤 4：模态参数提取

4.1 与 COV 步骤 5 相同：\(\mathbf{A}\) 的特征值 → 连续极点 (9)–(10) → \(f_r, \zeta_r\)；振型由 (11) 得到。

---

## 完整操作步骤（逐步）

输入：多通道输出 \(\mathbf{y}[k]\)，\(k=1,\ldots,N\)；采样率 \(f_s\)，\(\Delta t = 1/f_s\)。

通用预处理：

1. 去直流（每通道减均值）。
2. （可选）带通滤波至关心频段。
3. 选取块行数 \(i\)（如 20～50）与模型阶次 \(n\)（如 2× 模态数；可用稳定图细化）。

**SSI-COV：**

4. 用 (3) 估计 \(\ell = 0,\ldots,2i-1\) 的协方差 \(\mathbf{R}_\ell\)。
5. 构造块 Toeplitz \(\mathbf{T}_{1|i}\) (4)；SVD (5)；截断到 \(n\)；由 (6) 得到 \(\mathcal{O}_i\)。
6. 按 (7)(8) 恢复 \(\mathbf{A}\)、\(\mathbf{C}\)。
7. \(\mathbf{A}\) 的特征值 → 极点 → (9)–(10) 得 \(f_r, \zeta_r\)；振型 (11)。

SSI-DATA：

4. 构造块 Hankel \(\mathbf{Y}_{0|2i-1}\) (12)；划分过去/未来。
5. 形成（加权）投影并做 SVD (15)；截断到 \(n\)；由 (16) 得到 \(\mathcal{O}_i\)。
6. 按 (17) 恢复 \(\mathbf{A}\)、\(\mathbf{C}\)。
7. 与 COV 相同的模态提取。

稳定图（两种共用）：对若干阶次 \(n = n_{\min}, \ldots, n_{\max}\) 重复识别，绘制极点（如频率–阶次）。**稳定** 极点（频率与阻尼随 \(n\) 变化很小）对应物理模态；漂移或散点多为虚假，仅保留稳定极点（如相邻阶次间 \(\Delta f/f < 1\%\)，\(\Delta\zeta/\zeta < 5\%\)）。

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
