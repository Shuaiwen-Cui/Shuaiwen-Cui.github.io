# Stochastic Subspace Identification (SSI)

This page expands **2.8 SSI** from the [SHM roadmap](../shm.en.md): **Stochastic Subspace Identification (SSI)** is one of the most important and widely used methods in structural health monitoring for output-only (or input–output) modal identification. From output time series, SSI builds block Hankel/Toeplitz matrices, obtains the **observability** (and controllability) subspace via **projection** and **SVD**, recovers a discrete state-space model, and extracts **natural frequency**, **damping ratio**, and **mode shapes**. The two main variants are **SSI-COV** (covariance-driven), which uses output covariance lags to form a block Toeplitz matrix, and **SSI-DATA** (data-driven), which uses the raw data block Hankel and oblique projection. Both yield the same type of state-space model and modal parameters; the choice between them is mainly implementation, numerical, and data-usage trade-offs.

---

## Concept

**Core idea** — Under **ambient** (unknown or unmeasured) excitation, the structure’s response can be modelled as the output of a **linear time-invariant (LTI) discrete state-space system** driven by white noise. The **observability matrix** of this system spans a subspace that can be recovered from output-only data. **SSI** does this by:

1. **SSI-COV:** Estimating output **covariance matrices** at several lags, arranging them into a **block Toeplitz** matrix, and taking an SVD to obtain the observability range space; the state matrix \(\mathbf{A}\) and output matrix \(\mathbf{C}\) are then recovered from the observability matrix.
2. **SSI-DATA:** Building a **block Hankel matrix** from the output data, splitting it into “past” and “future” blocks, computing the **oblique projection** of the future row space onto the past row space, and taking an SVD of (a weighted form of) this projection to obtain the observability matrix; \(\mathbf{A}\) and \(\mathbf{C}\) are recovered as in COV.

From the discrete state matrix \(\mathbf{A}\), eigenvalues give **discrete poles**; these are converted to **continuous-time poles** (natural frequency \(\omega_r\), damping ratio \(\zeta_r\)). From \(\mathbf{C}\) and the eigenvectors of \(\mathbf{A}\), **mode shapes** at the sensors are obtained. No knowledge of the input force is required (output-only operational modal analysis).

**State-space model (discrete time)** — The underlying LTI model is

\[\mathbf{x}[k+1] = \mathbf{A} \mathbf{x}[k] + \mathbf{w}[k], \qquad \mathbf{y}[k] = \mathbf{C} \mathbf{x}[k] + \mathbf{v}[k]. \tag{1}\]

**Notation for (1):**

- \(\mathbf{x}[k] \in \mathbb{R}^n\): state vector; \(n\) is the **model order** (state dimension).
- \(\mathbf{A} \in \mathbb{R}^{n \times n}\): state matrix; its eigenvalues are the discrete poles; determines modal frequency and damping.
- \(\mathbf{C} \in \mathbb{R}^{p \times n}\): output matrix; maps state to measured outputs; used for mode shapes.
- \(\mathbf{y}[k] \in \mathbb{R}^p\): output vector (e.g. accelerations at \(p\) sensors).
- \(\mathbf{w}[k], \mathbf{v}[k]\): process and measurement noise (often assumed white, zero mean); input is implicitly in \(\mathbf{w}\).

**Observability** — The **observability matrix** (block rows) is

\[\mathcal{O}_i = \begin{bmatrix} \mathbf{C} \\ \mathbf{C}\mathbf{A} \\ \vdots \\ \mathbf{C}\mathbf{A}^{i-1} \end{bmatrix}. \tag{2}\]

SSI recovers a basis for the **range of \(\mathcal{O}_i\)** (up to a similarity transform). From \(\mathcal{O}_i\) one can obtain \(\mathbf{A}\) and \(\mathbf{C}\): the first block row is \(\mathbf{C}\); \(\mathbf{A}\) is recovered from the shift structure (e.g. \(\mathcal{O}_{i,\mathrm{up}}^\dagger \mathcal{O}_{i,\mathrm{down}}\), where “up” drops the last block row and “down” drops the first).

**SSI-COV vs SSI-DATA** — **SSI-COV** works on **covariance estimates**: it forms a block Toeplitz matrix from output covariances \(\mathbf{R}_\ell = \mathbb{E}[\mathbf{y}_{k+\ell}\mathbf{y}_k^T]\) (or sample estimates). One SVD of this matrix (or a factor of it) yields the observability. **SSI-DATA** works on **raw data**: it builds a block Hankel matrix from \(\mathbf{y}[k]\), computes the orthogonal/oblique projection of “future” outputs onto “past” outputs, and takes an SVD of this projection to get the observability. COV compresses data into covariances first (fewer numbers, fixed size); DATA uses the full data block and can use weighting (e.g. CVA) for better numerical and statistical properties. Both end at the same place: observability → \(\mathbf{A},\mathbf{C}\) → modal parameters.

---

## SSI-COV: Covariance-driven algorithm

### Step 1: Data preparation and covariance estimation

Input: Multi-channel output time series \(\mathbf{y}[k] \in \mathbb{R}^p\), \(k = 1,\ldots,N\), sampling interval \(\Delta t\), \(f_s = 1/\Delta t\). Typically \(p \geq 2\) (preferably more) sensors and \(N\) large (e.g. tens of thousands to millions) for stable covariance estimates.

Preprocessing:

1.1 Remove DC: \(\mathbf{y}[k] \leftarrow \mathbf{y}[k] - \frac{1}{N}\sum_{j=1}^N \mathbf{y}[j]\) per channel.

1.2 (Optional) Detrend and/or band-pass filter to the frequency band of interest.

**Covariance estimation:**

1.3 Choose the **number of block rows** \(i\) (e.g. 20–80). The number of lags needed is at least \(2i-1\). For lag \(\ell = 0,1,\ldots,2i-1\), estimate the **output covariance**:

\[\mathbf{R}_\ell = \frac{1}{N-\ell} \sum_{k=1}^{N-\ell} \mathbf{y}[k+\ell] \mathbf{y}[k]^T. \tag{3}\]

So \(\mathbf{R}_0\) is the (sample) output covariance at lag 0; \(\mathbf{R}_\ell\) is the cross-covariance between \(\mathbf{y}[k+\ell]\) and \(\mathbf{y}[k]\).

Practical notes:

- \(i\) should be large enough so that \(\mathcal{O}_i\) has full column rank for the desired order \(n\) (typically \(i \cdot p \geq n\), and \(i\) such that \(\mathbf{A}^{i-1}\) has decayed). Often \(i = 20\)–\(50\).
- Covariance length: \(N\) must be much larger than the maximum lag (\(N \gg 2i\)) for low-variance estimates; for civil structures (e.g. 1–10 Hz modes, 100 Hz \(f_s\)), 2–10 minutes of data is common.

### Step 2: Build the block Toeplitz matrix

2.1 Form the **block Toeplitz matrix** \(\mathbf{T}_{1|i} \in \mathbb{R}^{pi \times pi}\) from the covariances:

\[\mathbf{T}_{1|i} = \begin{bmatrix} \mathbf{R}_1 & \mathbf{R}_2 & \cdots & \mathbf{R}_i \\ \mathbf{R}_2 & \mathbf{R}_3 & \cdots & \mathbf{R}_{i+1} \\ \vdots & \vdots & \ddots & \vdots \\ \mathbf{R}_i & \mathbf{R}_{i+1} & \cdots & \mathbf{R}_{2i-1} \end{bmatrix}. \tag{4}\]

(Some conventions use \(\mathbf{R}_0\) in the first block column; the exact indexing may shift by one. The key is that the matrix has a Toeplitz structure and factorizes as observability × controllability.)

2.2 Under the state-space model (1), \(\mathbf{T}_{1|i} = \mathcal{O}_i \, \mathcal{C}_i\), where \(\mathcal{O}_i\) is the observability matrix (2) and \(\mathcal{C}_i\) is the controllability matrix (in the same block structure). So the **column space** of \(\mathbf{T}_{1|i}\) equals the column space of \(\mathcal{O}_i\).

### Step 3: SVD and observability

3.1 Compute the SVD of the block Toeplitz matrix:

\[\mathbf{T}_{1|i} = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^T. \tag{5}\]

3.2 Truncate to order \(n\): Choose \(n\) (e.g. 2× number of physical modes, or from a stabilisation diagram). Keep the first \(n\) singular values and vectors: \(\mathbf{U}_n \in \mathbb{R}^{pi \times n}\), \(\mathbf{\Sigma}_n \in \mathbb{R}^{n \times n}\), \(\mathbf{V}_n \in \mathbb{R}^{pi \times n}\).

3.3 **Observability matrix:** Set

\[\mathcal{O}_i = \mathbf{U}_n \mathbf{\Sigma}_n^{1/2}. \tag{6}\]

Then \(\mathcal{O}_i\) has \(i\) block rows of size \(p \times n\) each; the **first block row** is the output matrix \(\mathbf{C}\).

### Step 4: Recover A and C

4.1 Define shifted observability matrices:

- \(\mathcal{O}_{i,\mathrm{up}}\): \(\mathcal{O}_i\) **without the last block row** (\((i-1)p \times n\)).
- \(\mathcal{O}_{i,\mathrm{down}}\): \(\mathcal{O}_i\) **without the first block row** (\((i-1)p \times n\)).

From the state-space structure, \(\mathcal{O}_{i,\mathrm{down}} = \mathcal{O}_{i-1,\mathrm{up}} \mathbf{A}\) (with \(\mathcal{O}_{i-1,\mathrm{up}}\) the first \(i-1\) block rows of \(\mathcal{O}_i\)). So

\[\mathbf{A} = \mathcal{O}_{i,\mathrm{up}}^\dagger \, \mathcal{O}_{i,\mathrm{down}}, \tag{7}\]

where \(\dagger\) denotes the Moore–Penrose pseudo-inverse.

4.2 Output matrix: \(\mathbf{C}\) is the **first block row** of \(\mathcal{O}_i\):

\[\mathbf{C} = \mathcal{O}_i(1:p, 1:n). \tag{8}\]

### Step 5: Modal parameter extraction (common to COV and DATA)

5.1 Discrete eigenvalues: Compute eigenvalues \(\mu_j\) of \(\mathbf{A}\), \(j=1,\ldots,n\). Discard those outside the unit circle or clearly non-physical (e.g. real negative).

5.2 Continuous-time poles: For each eigenvalue \(\mu_j\),

\[\lambda_j = \frac{\ln(\mu_j)}{\Delta t} \quad \text{(principal value)}. \tag{9}\]

For a **conjugate pair** \(\lambda_r, \lambda_r^*\) corresponding to one physical mode:

\[\omega_r = |\lambda_r|, \qquad \zeta_r = -\frac{\mathrm{Re}(\lambda_r)}{|\lambda_r|}. \tag{10}\]

Natural frequency \(f_r = \omega_r/(2\pi)\).

5.3 Mode shapes: Let \(\mathbf{\psi}_r \in \mathbb{C}^n\) be the right eigenvector of \(\mathbf{A}\) for \(\mu_r\). The **mode shape** at the sensors is

\[\mathbf{\phi}_r = \mathbf{C} \mathbf{\psi}_r. \tag{11}\]

Normalise (e.g. unit norm or max component 1). For real structures, take real part or magnitude if needed.

---

## SSI-DATA: Data-driven algorithm

### Step 1: Data preparation and block Hankel matrix

Input: Same as COV: \(\mathbf{y}[k] \in \mathbb{R}^p\), \(k = 1,\ldots,N\), \(\Delta t\). Preprocessing: remove DC; optional detrend and band-pass filter.

Block Hankel matrix:

1.1 Choose **number of block rows** \(i\) (past) and \(i\) (future), and **number of columns** \(j\). Typically \(j = N - 2i + 1\) (use all available data). Require \(j \geq 2i\) and large enough for statistical stability.

1.2 Build the **block Hankel matrix** \(\mathbf{Y}_{0|2i-1} \in \mathbb{R}^{2pi \times j}\):

\[\mathbf{Y}_{0|2i-1} = \frac{1}{\sqrt{j}} \begin{bmatrix} \mathbf{y}[1] & \mathbf{y}[2] & \cdots & \mathbf{y}[j] \\ \mathbf{y}[2] & \mathbf{y}[3] & \cdots & \mathbf{y}[j+1] \\ \vdots & \vdots & \ddots & \vdots \\ \mathbf{y}[2i] & \mathbf{y}[2i+1] & \cdots & \mathbf{y}[2i+j-1] \end{bmatrix}. \tag{12}\]

The factor \(1/\sqrt{j}\) is for normalisation. Split into:

- **Past:** \(\mathbf{Y}_{0|i-1}\) — first \(i\) block rows (\(pi \times j\)).
- **Future:** \(\mathbf{Y}_{i|2i-1}\) — next \(i\) block rows (\(pi \times j\)).

### Step 2: Projection

2.1 The **oblique projection** of the row space of the future onto the row space of the past is

\[\mathbf{O}_i = \mathbf{Y}_{i|2i-1} \, \mathbf{Y}_{0|i-1}^T \bigl( \mathbf{Y}_{0|i-1} \mathbf{Y}_{0|i-1}^T \bigr)^{-1} \mathbf{Y}_{0|i-1}. \tag{13}\]

In practice one does **not** form this explicitly. Instead, use the **RQ** or **LQ** factorisation of the stacked past/future matrix, or the **SVD** of a weighted combination, so that the projection is represented in factored form.

2.2 Weighted projection (e.g. CVA — Canonical Variate Analysis): Define

\[\mathbf{W}_p = \bigl( \mathbf{Y}_{0|i-1} \mathbf{Y}_{0|i-1}^T \bigr)^{-1/2}, \qquad \mathbf{W}_f = \bigl( \mathbf{Y}_{i|2i-1} \mathbf{Y}_{i|2i-1}^T \bigr)^{-1/2}. \tag{14}\]

Then compute the SVD of the **weighted projection**:

\[\mathbf{W}_f \, \mathbf{Y}_{i|2i-1} \, \mathbf{Y}_{0|i-1}^T \, \mathbf{W}_p^T = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^T. \tag{15}\]

(Other weightings exist; unweighted corresponds to \(\mathbf{W}_p = \mathbf{I}\), \(\mathbf{W}_f = \mathbf{I}\).)

2.3 **Observability:** Truncate to order \(n\): \(\mathbf{U}_n, \mathbf{\Sigma}_n\). Then

\[\mathcal{O}_i = \mathbf{W}_f^{-1} \mathbf{U}_n \mathbf{\Sigma}_n^{1/2}. \tag{16}\]

So the observability matrix is recovered from the left singular vectors (and singular values) of the weighted projection, then “un-weighted” by \(\mathbf{W}_f^{-1}\).

### Step 3: Recover A and C

3.1 Same as SSI-COV: from \(\mathcal{O}_i\), form \(\mathcal{O}_{i,\mathrm{up}}\) and \(\mathcal{O}_{i,\mathrm{down}}\); then

\[\mathbf{A} = \mathcal{O}_{i,\mathrm{up}}^\dagger \, \mathcal{O}_{i,\mathrm{down}}, \qquad \mathbf{C} = \mathcal{O}_i(1:p, 1:n). \tag{17}\]

### Step 4: Modal parameter extraction

4.1 Identical to COV Step 5: eigenvalues of \(\mathbf{A}\) → continuous poles (9)–(10) → \(f_r, \zeta_r\); mode shapes from (11).

---

## Complete procedure (step-by-step)

Input: Multi-channel output \(\mathbf{y}[k]\), \(k=1,\ldots,N\); sampling rate \(f_s\), \(\Delta t = 1/f_s\).

Common preprocessing:

1. Remove DC (subtract mean per channel).
2. (Optional) Band-pass filter to the band of interest.
3. Choose block row number \(i\) (e.g. 20–50) and model order \(n\) (e.g. 2× number of modes; refine with stabilisation diagram).

**SSI-COV:**

4. Estimate covariances \(\mathbf{R}_\ell\) for \(\ell = 0,\ldots,2i-1\) using (3).
5. Build block Toeplitz \(\mathbf{T}_{1|i}\) (4); SVD (5); truncate to \(n\); form \(\mathcal{O}_i\) (6).
6. Recover \(\mathbf{A}\) (7) and \(\mathbf{C}\) (8).
7. Eigenvalues of \(\mathbf{A}\) → poles → \(f_r, \zeta_r\) (9)–(10); mode shapes (11).

**SSI-DATA:**

4. Build block Hankel \(\mathbf{Y}_{0|2i-1}\) (12); split past/future.
5. Form (weighted) projection and compute SVD (15); truncate to \(n\); form \(\mathcal{O}_i\) (16).
6. Recover \(\mathbf{A}\) and \(\mathbf{C}\) (17).
7. Same modal extraction as COV.

Stabilisation diagram (both): For several orders \(n = n_{\min}, \ldots, n_{\max}\), repeat the identification. Plot poles (e.g. frequency vs order). **Stable** poles (frequency and damping change little with \(n\)) correspond to physical modes; drifting or scattered poles are often spurious. Keep only stable poles (e.g. \(\Delta f/f < 1\%\), \(\Delta\zeta/\zeta < 5\%\) between consecutive orders).

Output: Modal parameters — natural frequencies \(f_r\), damping ratios \(\zeta_r\), mode shapes \(\mathbf{\phi}_r\).

---

## When to use and limitations

Use when:

- **Output-only** (ambient vibration) or input–output data; no need to measure excitation.
- **Robust, high-quality** modal identification is required; SSI is a standard in civil and mechanical OMA.
- **Multiple modes**, including **close modes**, need to be separated (parametric state-space model and stabilisation diagram help).
- **Damping** is required (SSI gives it directly, unlike FDD which needs EFDD).
- **Long, multi-channel** data are available (SSI benefits from large \(N\) and \(p\)).

Limitations:

- **Order selection:** Model order \(n\) must be chosen; too small misses modes, too large adds spurious poles. Use **stabilisation diagram** over a range of \(n\) and keep only stable poles.
- **Computation and memory:** Block matrices are large (\(pi \times pi\) or \(2pi \times j\)); SVD and factorisations are costly; less suited to very resource-limited edge devices without simplification or offload.
- **Excitation:** Assumes persistent excitation (broadband or at least covering the modes of interest); very narrowband or non-stationary excitation can bias results.
- **Noise:** Heavy noise increases variance of covariance/projection estimates; long data and multiple channels improve robustness.

---

## Engineering practice and reproducible parameters

| Aspect | Guideline and typical values |
|--------|-------------------------------|
| Block rows \(i\) | Large enough so \(\mathcal{O}_i\) has full column rank for order \(n\); typically \(i = 20\)–\(50\); \(i \cdot p \geq n\). |
| Model order \(n\) | Use stabilisation diagram: try \(n = 8, 10, 12, \ldots\) up to ~2× expected number of modes; keep poles stable across orders (e.g. \(\Delta f/f < 1\%\), \(\Delta\zeta/\zeta < 5\%\)). |
| Data length \(N\) | \(N \gg 2i\); for civil (e.g. 1–10 Hz), 2–10 min at 100 Hz is common; longer for low damping or low SNR. |
| COV vs DATA | COV: compresses to covariances, one Toeplitz + one SVD; DATA: full data Hankel, projection + SVD; DATA can use CVA weighting for better conditioning; both give same type of result. |
| Numerical stability | Avoid singular or ill-conditioned matrices (e.g. \(\mathbf{Y}_{0|i-1}\mathbf{Y}_{0|i-1}^T\) in DATA); use regularisation or pseudo-inverse; drop tiny singular values. |
| Validation | Compare with FDD/EFDD or ERA; recompute covariances or fit from identified model; check mode shapes for consistency. |

---

## Edge and online computing

Suitability — **Moderately to poorly suited** for low-power edge: large matrices (Hankel/Toeplitz, SVD), multiple orders for stabilisation diagram, and long data buffers increase memory and compute. Best used on server/cloud or powerful gateways; edge can do lightweight screening (e.g. PP or low-order ARX) and upload segments for SSI.

Implementation strategy:

1. Reduce problem size: Limit \(i\) and \(n\); use fewer channels or a subset of lags (COV); limit \(j\) (DATA).
2. COV on edge: If covariances are computed online (sliding window or recursive), one can run SSI-COV with fixed \(n\) and small \(i\) to reduce SVD size; full stabilisation diagram offline.
3. DATA on edge: Block Hankel and projection are memory-heavy; consider uploading raw or downsampled data and running SSI in the cloud.
4. Tiered: Edge: quick checks (e.g. PP, simple correlation); upload suspect or periodic segments; cloud: full SSI and stabilisation diagram.

Challenges: Memory (Hankel/Toeplitz), CPU (SVD), order selection (stabilisation diagram), and data volume for long records.

---

## Comparison with other methods

SSI vs FDD/EFDD:

- **SSI:** Time-domain, parametric, state-space; gives frequency, damping, mode shape in one framework; needs order selection and stabilisation diagram; better for close modes.
- **FDD:** Frequency-domain, non-parametric; no damping (EFDD adds it); no order choice; simpler but less accurate for damping and close modes.

SSI vs ERA / NExT-ERA:

- **SSI:** Builds Hankel/Toeplitz from **covariances (COV)** or **raw data (DATA)** and uses projection + SVD; very robust, standard in OMA; larger matrices.
- **ERA:** Builds Hankel from **Markov parameters** (impulse response or NExT correlation); one SVD, compact; NExT-ERA is output-only via correlation. SSI is generally more robust and flexible for output-only; ERA is lighter when impulse/correlation is already available.

SSI-COV vs SSI-DATA:

- **COV:** Single compression to covariances; one Toeplitz + SVD; lower memory if \(N\) is huge; covariance estimates can be noisier.
- **DATA:** Uses full data block; projection step; can use CVA or other weighting for better numerical/statistical properties; more flexible, often preferred in software; higher memory.

When to choose:

- Prefer **SSI** when you need the best balance of accuracy, damping, and close-mode separation for output-only OMA; use **SSI-COV** for simpler implementation and fixed memory, **SSI-DATA** (with CVA) for maximum robustness and quality when resources allow.
