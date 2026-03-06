# 特征系统实现算法（ERA）/ NExT-ERA

本页对 [SHM roadmap](../shm.zh.md) 中的 **2.7 ERA / NExT-ERA** 做展开： **ERA** 从脉冲响应或自由衰减构造 **Markov 参数** 与 **Hankel 矩阵**，经 **SVD** 与最小实现得到状态空间并提取模态参数； **NExT-ERA** 在环境激励下用输出的自/互相关构造等效脉冲响应，再对之应用 ERA，实现仅输出数据的运行模态分析。

---

## 教学视频

<iframe width="800" height="450" src="https://www.youtube-nocookie.com/embed/5TM_JrW5qDo" frameborder="0" allowfullscreen></iframe>

---

## 概念

**核心思路** — **ERA**（Eigensystem Realization Algorithm）是时域辨识方法：从脉冲响应或自由衰减序列得到 Markov 参数，排成 Hankel 矩阵，做 SVD 并截断得到最小实现，即离散状态空间 \((\mathbf{A},\mathbf{B},\mathbf{C})\)；由 \(\mathbf{A}\) 的特征值得连续时间极点（固有频率、阻尼比），由 \(\mathbf{C}\) 与特征向量得振型。 **NExT**（Natural Excitation Technique）指出：在白噪声或宽频环境激励下，两通道响应的互相关与脉冲响应满足同一齐次微分方程，故仅输出数据估计的自/互相关可当作“等效脉冲响应”代入 ERA，即 **NExT-ERA**。流程：ERA 为 脉冲/自由衰减 → Markov → Hankel → SVD → 状态空间 → 模态；NExT-ERA 为 仅输出 → 相关估计 → 等效脉冲 → ERA。

**原理（简要）** — 离散 LTI 状态空间模型：

\[\mathbf{x}[k+1] = \mathbf{A} \mathbf{x}[k] + \mathbf{B} u[k], \qquad \mathbf{y}[k] = \mathbf{C} \mathbf{x}[k] + \mathbf{D} u[k]. \tag{1}\]

脉冲响应（Markov 参数）为 \(\mathbf{h}_0 = \mathbf{D}\)，\(\mathbf{h}_k = \mathbf{C} \mathbf{A}^{k-1} \mathbf{B}\)（\(k \geq 1\)）。将 \(\mathbf{h}_k\) 按块组成 Hankel 矩阵 \(\mathbf{H}\)，其秩与可观测/可控子空间相关；对 \(\mathbf{H}\) 做 SVD 并截断到目标阶次 \(n\)，可从 \(\mathbf{U},\mathbf{\Sigma},\mathbf{V}\) 恢复 \(\mathbf{A},\mathbf{B},\mathbf{C}\)（最小实现）。\(\mathbf{A}\) 的特征值 \(\mu_i\) 对应离散极点；变换到连续时间 \(\lambda_i = \ln(\mu_i)/\Delta t\)，即得模态极点，进而得到 \(\omega_r\)、\(\zeta_r\) 与振型 \(\mathbf{\phi}_r\)。

**式 (1) 符号含义：**

- \(\mathbf{x}[k] \in \mathbb{R}^n\)：状态向量；\(n\) 为状态维数（模型阶次）。

- \(\mathbf{A} \in \mathbb{R}^{n \times n}\)：状态矩阵；其特征值决定模态频率与阻尼。

- \(\mathbf{B},\mathbf{C},\mathbf{D}\)：输入、输出与直通矩阵；脉冲响应由 \(\mathbf{h}_k = \mathbf{C} \mathbf{A}^{k-1} \mathbf{B}\) 给出。

- \(u[k], \mathbf{y}[k]\)：输入与输出（多通道时 \(\mathbf{y}\) 为向量）。

**ERA vs NExT-ERA** — **ERA** 需要脉冲响应或自由衰减（来自冲击试验或已知输入）； **NExT-ERA** 只需环境激励下的响应，用相关函数代替脉冲响应后同样做 Hankel + SVD + 最小实现，故适用于仅输出的运行模态分析（如桥梁、楼宇等）。

---

## 详细算法与实现

### 步骤 1：数据准备与“脉冲/等效脉冲”序列

**ERA（有脉冲或自由衰减）：**

1.1  **输入：** 脉冲响应序列 \(\mathbf{h}_0, \mathbf{h}_1, \ldots, \mathbf{h}_L\)（多通道时为矩阵序列），或自由衰减响应经识别得到的等效 Markov 参数；采样间隔 \(\Delta t\)。实践中 \(L\) 建议覆盖 2～3 倍最关心模态的衰减时间（例如阻尼比 1% 时衰减时间约 \(1/(\zeta\omega_n)\) 量级），按 \(L \approx 2\sim3 \times (1/(\zeta \omega_n)) \times f_s\) 取整。

1.2  **预处理：** 去直流（每通道减均值）；可选去趋势与带通滤波（保留关心频段）；若为自由衰减，可先做 RDT 等得到更干净衰减再转为 Markov 形式。

**NExT-ERA（仅输出）：**

1.3  **输入：** 多通道响应时间序列 \(\{y_i[k]\}\)，\(i=1,\ldots,n\)，\(k=1,\ldots,N\)，采样率 \(f_s\)，\(\Delta t = 1/f_s\)。建议 \(N\) 足够大以使相关估计稳定（例如至少数万样本，或 2～5 分钟以上数据）。

1.4  **估计相关函数：** 选定参考通道（如 \(y_{\mathrm{ref}}\)，建议选信噪比高且不在关心模态节点上的通道）。对每通道 \(i\) 估计与参考的互相关（或自相关）：

\[R_{i,\mathrm{ref}}(\tau) = \frac{1}{N-\tau} \sum_{k=1}^{N-\tau} y_i[k] \, y_{\mathrm{ref}}[k+\tau], \qquad \tau = 0,1,\ldots,L. \tag{2}\]

在白噪声激励下，\(R_{i,\mathrm{ref}}(\tau)\) 与“参考到通道 \(i\) 的脉冲响应”成比例，故可将 \(R_{i,\mathrm{ref}}(0), R_{i,\mathrm{ref}}(1), \ldots\) 视为等效 Markov 参数用于后续 ERA。

**可复现细节（NExT-ERA）：** 单参考、多输出时，第 \(k\) 个 Markov 块可取为列向量 \(\mathbf{Y}_k = [R_{1,\mathrm{ref}}(k), \ldots, R_{n,\mathrm{ref}}(k)]^T\)（\(n\) 为通道数）；若需矩阵形式，可将 \(\mathbf{Y}_k\) 视为 \(n\times 1\) 块，Hankel 中每块为上述列向量。

**实践要点：**

- 参考通道：选信号质量好、且不在关心模态节点上的通道；可试多个参考比较 **稳定图** 质量。

- 相关长度 \(L\)：需覆盖关心模态的衰减时间；典型为数百至数千个采样点（如 \(f_s=100\) Hz、阻尼 1%、最低模态约 1 Hz 时，衰减时间约 15 s，\(L \approx 1500\)）。

- NExT 假设：激励近似白噪声或宽频；若激励有色，相关与脉冲响应比例关系会受影响。

### 步骤 2：构造 Hankel 矩阵

2.1  **由 Markov 参数组块：** 设 \(\mathbf{Y}_k\) 为第 \(k\) 个 Markov 参数块（多通道时为一矩阵或列向量）。取块行数 \(r\)、块列数 \(s\)（\(r,s\) 足够大以覆盖模态；典型取 \(r=s=20\sim50\)，且 \(r+s-1 \leq L+1\) 以保证 \(\mathbf{Y}_{r+s-1}\) 可用），构造：

\[\mathbf{H}(0) = \begin{bmatrix} \mathbf{Y}_1 & \mathbf{Y}_2 & \cdots & \mathbf{Y}_s \\ \mathbf{Y}_2 & \mathbf{Y}_3 & \cdots & \mathbf{Y}_{s+1} \\ \vdots & \vdots & \ddots & \vdots \\ \mathbf{Y}_r & \mathbf{Y}_{r+1} & \cdots & \mathbf{Y}_{r+s-1} \end{bmatrix}, \qquad \mathbf{H}(1) = \begin{bmatrix} \mathbf{Y}_2 & \mathbf{Y}_3 & \cdots & \mathbf{Y}_{s+1} \\ \mathbf{Y}_3 & \mathbf{Y}_4 & \cdots & \mathbf{Y}_{s+2} \\ \vdots & \vdots & \ddots & \vdots \\ \mathbf{Y}_{r+1} & \mathbf{Y}_{r+2} & \cdots & \mathbf{Y}_{r+s} \end{bmatrix}. \tag{3}\]

\(\mathbf{H}(1)\) 为 \(\mathbf{H}(0)\) 向下平移一块行得到。

2.2  **单通道/标量情形：** 若为标量脉冲响应 \(h_k\)，则 \(\mathbf{Y}_k = h_k\)，\(\mathbf{H}(0),\mathbf{H}(1)\) 为普通 Hankel 矩阵。

### 步骤 3：SVD 与最小实现

3.1  **对 \(\mathbf{H}(0)\) 做 SVD：** \(\mathbf{H}(0) = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^T\)（或 \(\mathbf{V}^H\)），\(\mathbf{\Sigma} = \mathrm{diag}(\sigma_1,\ldots,\sigma_{\min})\)。

3.2  **截断阶次 \(n\)：** 根据奇异值分布或稳定图，保留前 \(n\) 个奇异值（\(n\) 约为 2×模态数），得 \(\mathbf{U}_n, \mathbf{\Sigma}_n, \mathbf{V}_n\)。

3.3  **恢复状态空间：** 定义 \(\mathbf{O} = \mathbf{U}_n \mathbf{\Sigma}_n^{1/2}\)（可观性块），\(\mathbf{C}_c = \mathbf{\Sigma}_n^{1/2} \mathbf{V}_n^T\)（可控性块）。则

\[\mathbf{A} = \mathbf{\Sigma}_n^{-1/2} \mathbf{U}_n^T \mathbf{H}(1) \mathbf{V}_n \mathbf{\Sigma}_n^{-1/2}, \tag{4}\]

\(\mathbf{B}\) 取为 \(\mathbf{\Sigma}_n^{-1/2} \mathbf{U}_n^T\) 乘以 \(\mathbf{H}(0)\) 的首块列；\(\mathbf{C}\) 取为 \(\mathbf{H}(0)\) 的首块行乘以 \(\mathbf{V}_n \mathbf{\Sigma}_n^{-1/2}\)（标准 ERA 中，若 \(\mathbf{H}(0)\) 首列为 \([\mathbf{Y}_1;\mathbf{Y}_2;\ldots;\mathbf{Y}_r]\)，则 \(\mathbf{B} = \mathbf{\Sigma}_n^{-1/2}\mathbf{U}_n^T [\mathbf{Y}_1;\ldots;\mathbf{Y}_r]\)；首行为 \([\mathbf{Y}_1,\mathbf{Y}_2,\ldots,\mathbf{Y}_s]\) 则 \(\mathbf{C}\) 为该行右乘 \(\mathbf{V}_n\mathbf{\Sigma}_n^{-1/2}\)）。实现时注意块尺寸与矩阵维度一致。

**实现要点：**

- 阶次选择：多试几个 \(n\)（如 \(n=8,10,12,\ldots,2\times\) 预期模态数），做 **稳定图**（极点随 \(n\) 变化）；稳定极点对应真模态，漂移或散点多为虚假。稳定判据常用：频率相对变化 \(<\Delta f/f\)（如 1%）且阻尼相对变化 \(<\Delta\zeta/\zeta\)（如 5%）。

- 数值：\(\mathbf{H}(0)\) 可能病态；SVD 截断有正则化效果，必要时可略增大 \(n\) 再按稳定图筛模态；避免 \(\mathbf{\Sigma}_n\) 中过小奇异值（可设最小奇异值阈值或忽略过小奇异值）。

### 步骤 4：从状态空间提取模态参数

4.1  **离散特征值：** 对 \(\mathbf{A}\) 求特征值 \(\mu_i\)，\(i=1,\ldots,n\)。剔除不在单位圆内或明显为噪声的极点。

4.2  **连续极点：** \(\lambda_i = \ln(\mu_i) / \Delta t\)（主值）。对共轭对 \(\lambda_{r},\lambda_{r}^*\)：

\[\omega_r = |\lambda_r|, \qquad \zeta_r = -\mathrm{Re}(\lambda_r) / |\lambda_r|. \tag{5}\]

固有频率 \(f_r = \omega_r / (2\pi)\)，阻尼比 \(\zeta_r\)。

4.3  **振型：** 由 \(\mathbf{C}\) 与 \(\mathbf{A}\) 的右特征向量 \(\mathbf{\psi}_r\) 得到输出空间振型 \(\mathbf{\phi}_r \propto \mathbf{C} \mathbf{\psi}_r\)；按通道归一化（如单位范数或最大分量置 1）。

**实践要点：**

- 共轭对：物理模态对应成对共轭极点；单实极点多为数值或噪声，可剔除。

- 稳定图：随 \(n\) 增加，真模态频率/阻尼应稳定，虚假模态会漂移或突变；仅保留在连续多个阶次上稳定的极点。

---

## 完整操作步骤（逐步）

**输入（ERA）：** 脉冲响应或自由衰减序列；或 **（NExT-ERA）** 多通道响应 \(\{y_i[k]\}\)，采样率 \(f_s\)。

**步骤 1：预处理与“脉冲”序列**

1.1 若 ERA：确认脉冲/自由衰减序列与 \(\Delta t\)；可选去直流、滤波。

1.2 若 NExT-ERA：选参考通道；按 (2) 估计互相关 \(R_{i,\mathrm{ref}}(\tau)\)，\(\tau=0,\ldots,L\)；将 \(R\) 视为等效 Markov 序列。

**步骤 2：构造 Hankel**

2.1 由 Markov 参数（或等效）组块，形成 \(\mathbf{H}(0)\) 与 \(\mathbf{H}(1)\)。

**步骤 3：SVD 与实现**

3.1 对 \(\mathbf{H}(0)\) 做 SVD；按奇异值或稳定图选取阶次 \(n\)，截断得 \(\mathbf{U}_n, \mathbf{\Sigma}_n, \mathbf{V}_n\)。

3.2 按 (4) 及标准 ERA 公式得到 \(\mathbf{A},\mathbf{B},\mathbf{C}\)。

**步骤 4：模态提取**

4.1 求 \(\mathbf{A}\) 的特征值 \(\mu_i\)；转为连续极点 \(\lambda_i\)，按 (5) 得 \(\omega_r, \zeta_r\)。

4.2 由 \(\mathbf{C}\) 与 \(\mathbf{A}\) 特征向量得振型 \(\mathbf{\phi}_r\)，归一化。

**输出：** 模态参数 — 固有频率 \(f_r\)，阻尼比 \(\zeta_r\)，振型 \(\mathbf{\phi}_r\)。

**可复现操作清单（NExT-ERA 示例）：**

1. 准备：\(n_y\) 通道响应 \(y_i[k]\)，\(k=1,\ldots,N\)，采样率 \(f_s\)；去直流（每通道减均值）；可选带通滤波。
2. 选参考通道 \(y_{\mathrm{ref}}\)；取 \(L = \lfloor 2.5 \cdot f_s / (\zeta \omega_n) \rfloor\)（\(\zeta\approx 0.01\)，\(\omega_n\) 取关心频段中心弧度频率的粗估）。
3. 按式 (2) 计算 \(R_{i,\mathrm{ref}}(\tau)\)，\(\tau=0,\ldots,L\)；组成 \(\mathbf{Y}_k = [R_{1,\mathrm{ref}}(k),\ldots,R_{n_y,\mathrm{ref}}(k)]^T\)。
4. 取 \(r=s=30\)（或 20～50），构造 \(\mathbf{H}(0),\mathbf{H}(1)\)；对 \(\mathbf{H}(0)\) 做 SVD。
5. 试 \(n=8,10,12,\ldots,24\)，对每个 \(n\) 用式 (4) 得 \(\mathbf{A}\)，求特征值并转为 \(f_r,\zeta_r\)；画稳定图，保留频率/阻尼随 \(n\) 稳定的极点（如 \(\Delta f/f<1\%\)）。
6. 用式 (5) 与 \(\mathbf{C}\)、特征向量得振型并归一化；可选：用识别模型重算相关与原始 \(R\) 对比验证。

---

## 适用情形与局限

**适用：**

- **ERA**：有脉冲响应或自由衰减（冲击试验、阶跃松弛等）；需要时域、一次 SVD、无迭代的辨识。
- **NExT-ERA**：仅输出、环境激励（桥梁、楼宇、大型结构）；激励近似宽频/白噪声时效果较好。
- 多通道：多测点可得振型；参考通道选择影响质量。

**局限：**

- 阶次敏感：\(n\) 过小漏模态，过大引入虚假极点；需稳定图或多次试算。
- NExT 假设：白噪声/宽频激励；窄带或有色激励会偏置相关与模态估计。
- 相关估计：需要足够数据与滞后长度；长相关增加 Hankel 规模与计算量。
- 密模态：与 FDD/SSI 类似，模态非常接近时分离难度增加。

---

## 工程实践要点与可复现参数

| 方面 | 建议与典型取值 |
|------|----------------|
| 参考通道（NExT-ERA） | 选信噪比高、且不在关心模态节点上的通道；可试多个参考比较稳定图。 |
| 相关长度 \(L\) | 覆盖 2–3 倍最关心模态的衰减时间；典型 \(L = \lfloor 2.5/(\zeta \omega_n) \cdot f_s \rfloor\)，\(\zeta\) 用先验或 0.01；过长则 Hankel 大、计算贵。 |
| 块行/块列 \(r,s\) | 取 \(r,s\) 足够大使 Hankel 秩能反映系统阶次；典型 \(r=s=20\sim50\)，且 \(r+s-1\le L+1\)。 |
| 阶次 \(n\) | 用稳定图：试 \(n=8,10,12,\ldots\) 至约 2×预期模态数；保留频率/阻尼随 \(n\) 稳定的极点（如 \(\Delta f/f<1\%\)，\(\Delta\zeta/\zeta<5\%\)）。 |
| 数值稳定性 | SVD 截断即正则化；忽略过小奇异值（如 \(\sigma_i/\sigma_1 < 10^{-6}\)），避免 \(\mathbf{A}\) 病态。 |
| 验证 | 用识别出的 \(\mathbf{A},\mathbf{B},\mathbf{C}\) 重算脉冲响应或相关，与原始数据对比；或与 FDD/SSI 结果交叉核对。 |

---

## 边缘与在线计算

**是否适用** — 较适用。仅输出（ **NExT-ERA** ）无需输入测量；单次 Hankel 构建 + 一次 SVD，无迭代；可通过限制 Hankel 规模与阶次 \(n\) 控制算力。

**实现策略：**

1. 相关估计：在线滑动窗或分块计算互相关；需缓冲与乘加运算，长度 \(L\) 不宜过大。
2. Hankel 与 SVD：限制 \(r,s\) 与 \(n\)；可只计算前若干奇异值（精简 SVD）。
3. 稳定图：边缘可只算少数几个 \(n\) 做简单稳定判据；完整稳定图可放云端或离线。

**挑战：**

- 内存：Hankel 与 SVD 矩阵规模随 \(r,s,n\) 增长；需限制以适配 MCU。
- 阶次选择：边缘上难以做完整稳定图；多依赖先验或固定 \(n\)。
- 相关估计：长数据与长滞后需要缓冲与计算量。

---

## 与其他方法比较

**ERA vs FDD：**

- **ERA**：时域、参数化、一次 SVD + 特征值；直接给出频率、阻尼、振型；需阶次选择。
- **FDD**：频域、非参数、各频率线 SVD；无阻尼（需 EFDD）；无需阶次，但密模态时分离依赖奇异值曲线。

**ERA vs SSI：**

- **ERA**：从 Markov/相关构造 Hankel，实现简洁；单次 SVD，计算量相对小。
- **SSI**：从数据块直接构造 Hankel/Toeplitz，投影 + SVD；通常更稳健、工程应用更广，但矩阵与计算更大。

**NExT-ERA vs 其他仅输出方法：**

- **NExT-ERA**：相关 → 等效脉冲 → ERA；结构清晰，适合有 RDT 等预处理时与 ERA 衔接。
- **FDD/EFDD**：频域 PSD + SVD；无需阶次，EFDD 可补阻尼。
- **SSI**：直接从数据构造、稳健性好；计算与内存需求更高。

**何时选择：**

- 有脉冲或自由衰减时优先考虑 **ERA**。

- 仅输出、环境激励且希望时域、一次 SVD、实现简单时，用 **NExT-ERA**；若需最强稳健性与完整稳定图，可考虑 SSI。
