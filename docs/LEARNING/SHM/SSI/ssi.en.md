# Stochastic Subspace Identification (SSI)

This page expands **2.8 SSI** from the [SHM roadmap](../shm.en.md). **Stochastic Subspace Identification (SSI)** is one of the most widely used **output-only** modal identification methods in structural health monitoring: it estimates **natural frequencies, damping ratios, and mode shapes** from multi-channel response time series **without measuring excitation**. SSI is the **stochastic, output-only** branch of **subspace methods**. The two common implementations are **SSI-COV** (covariance-driven) and **SSI-DATA** (data-driven); both yield the same type of state-space model and modal parameters, and differ mainly in how data are used and in numerical implementation.

---

## 1. Motivation and context

### 1.1 Why output-only identification?

- **Ambient excitation is unmeasured**: Excitation of large structures (bridges, buildings, turbines) comes from wind, traffic, ground motion, etc., and is rarely known or measured. Traditional modal analysis assumes known inputs (hammer, shaker), which is often impractical in the field.
- **Operational conditions**: We care about dynamic behaviour under **real service**; when only response sensors (acceleration, strain, etc.) are available, identification must rely **only on output time series**.
- **Goal**: Estimate **natural frequency \(f_r\), damping ratio \(\zeta_r\), and mode shapes \(\mathbf{\phi}_r\)** from multi-channel output \(\mathbf{y}[k]\) without knowing or measuring the input—i.e. **operational modal analysis (OMA)**.

### 1.2 Why a state-space model?

- **Unified framework**: The discrete LTI state-space \((\mathbf{A},\mathbf{C})\) describes both dynamics and observation; modal parameters follow from eigenvalues of \(\mathbf{A}\) and from \(\mathbf{C}\) and eigenvectors.
- **Compatibility with random excitation**: Modelling unknown excitation and measurement noise as **white noise**, the “white-noise driven, output-only” setup fits the state equation (1); the system remains LTI and tractable.
- **Identifiability**: Under observability, the **second-order statistics of the output** (covariances at several lags) are **uniquely determined** by \((\mathbf{A},\mathbf{C})\) **up to a similarity transform**, so the observability subspace can be **recovered** from covariances or data blocks, and then \(\mathbf{A},\mathbf{C}\). In short:
  - **Second-order statistics**: the multi-lag covariances \(\mathbf{R}_\ell = \mathbb{E}[\mathbf{y}[k+\ell]\mathbf{y}[k]^T]\); under white-noise drive, the full sequence is uniquely fixed by \((\mathbf{A},\mathbf{C})\) and noise covariances.
  - **Up to similarity**: State-space representations are not unique (change of coordinates \(\mathbf{x}'=\mathbf{T}\mathbf{x}\) leaves input–output behaviour unchanged), so we can only identify an equivalence class; **modal parameters** (eigenvalues, mode-shape directions) are invariant under similarity and are therefore uniquely determined.
  - **Why observability**: If the system were not observable, some state directions would never appear in the output and could not be recovered from output alone; observability ensures that the state is “visible” over multiple output steps, so the observability matrix has full column rank and its column space can be recovered from data.
  - **Chain**: \(\mathbf{R}_\ell = \mathbf{C}\mathbf{A}^{\ell-1}\mathbf{G}\) → block Toeplitz \(\mathbf{T}=\mathcal{O}_i\mathcal{C}_i\) → row space of \(\mathbf{T}\) = column space of \(\mathcal{O}_i\) → SVD gives observability subspace → shift relation gives \(\mathbf{A},\mathbf{C}\).

---

## 2. Core mathematical principles

### 2.1 One sentence and two basic principles

**One sentence**: **The observability subspace can be recovered from output data.**

**Two basic principles**:

1. **Data have a low-rank subspace structure**: The block matrix built from the output (block Toeplitz of covariances, or block Hankel/projection of data) has **rank \(n\)** (the system order) in the ideal case, and its row/column space is spanned by the observability matrix; an SVD of this matrix extracts a basis for the observability subspace.
2. **Shift invariance**: The block rows of the observability matrix satisfy “next block row = current block row × \(\mathbf{A}\)”. Once a basis for the observability subspace is available, this shift relation **uniquely recovers** \(\mathbf{A}\) and \(\mathbf{C}\).

Together: **first obtain the observability subspace from data (principle 1), then obtain \(\mathbf{A},\mathbf{C}\) from the shift relation (principle 2)**.

### 2.2 Why are “output covariances uniquely determined by the observability matrix”?

Under model (1) with stationarity and white-noise assumptions, the multi-lag covariances can be written

$$
\mathbf{R}_\ell = \mathbf{C} \mathbf{A}^{\ell-1} \mathbf{G}, \quad \ell \geq 1,
$$

where \(\mathbf{G}=\mathbb{E}[\mathbf{x}[k+1]\mathbf{y}[k]^T]\) is the state–output cross-covariance. So \(\mathbf{R}_\ell\) depends on the system only through “\(\mathbf{C}\mathbf{A}^{\ell-1}\)”, and \(\mathbf{C},\mathbf{C}\mathbf{A},\ldots\) are exactly the block rows of the observability matrix \(\mathcal{O}_i\). Arranging \(\mathbf{R}_1,\mathbf{R}_2,\ldots\) into the block Toeplitz \(\mathbf{T}_{1|i}\), each block row is a linear combination of the corresponding block rows of \(\mathcal{O}_i\) (coefficients from \(\mathbf{G}\)), hence **\(\mathbf{T}_{1|i} = \mathcal{O}_i \mathcal{C}_i\)** (\(\mathcal{C}_i\) built from \(\mathbf{G}\), etc.). So the **row space** of \(\mathbf{T}_{1|i}\) equals the **column space** of \(\mathcal{O}_i\) (when \(\mathcal{C}_i\) has full row rank). Thus: the observability matrix uniquely determines this row space; conversely, an SVD of \(\mathbf{T}_{1|i}\) uniquely recovers the observability subspace. The data-driven form (block Hankel + oblique projection) gives the same observability subspace in theory.

### 2.3 Essence of subspace methods

**Subspace methods**: Instead of fitting \((\mathbf{A},\mathbf{C})\) directly, one builds a block matrix from data (Hankel/Toeplitz, etc.), uses **SVD or projection + SVD** to extract **a basis for the observability subspace**, then uses the **shift structure** of the observability matrix to recover \(\mathbf{A},\mathbf{C}\) from that basis. In short: **first find a basis for the subspace, then recover the state-space parameters from that basis**.

**Deterministic vs stochastic (SSI) subspace methods**:

| Aspect | Deterministic (e.g. ERA) | Stochastic (SSI) |
|--------|--------------------------|------------------|
| Input | Known, measured | Unknown, treated as white noise |
| Data | Input–output or impulse response | Output only or output covariances |
| Block matrix | Often block Hankel from Markov parameters | Block Toeplitz (covariances) or block Hankel + projection |
| “Drive” side | Controllability \([\mathbf{B},\mathbf{A}\mathbf{B},\ldots]\) | Controllability-like \([\mathbf{G},\mathbf{A}\mathbf{G},\ldots]\) (no \(\mathbf{B}\)) |
| Typical algorithms | ERA, Ho–Kalman, deterministic N4SID | SSI-COV, SSI-DATA, stochastic N4SID |
| Scenario | Lab, known input | Field, output-only, OMA |

Both families share the same observability matrix and shift structure; they differ only in the source of data and how the block matrix is built. SSI is designed for “no known input, output-only” OMA.

---

## 3. Model and notation

### 3.1 State-space model

The discrete-time LTI model is

$$
\mathbf{x}[k+1] = \mathbf{A} \mathbf{x}[k] + \mathbf{w}[k], \qquad \mathbf{y}[k] = \mathbf{C} \mathbf{x}[k] + \mathbf{v}[k]. \tag{1}
$$

The state evolves as above; we only observe \(\mathbf{y}[k]\). The goal of SSI is to recover \(\mathbf{A},\mathbf{C}\) from output data and then obtain modal parameters.

**Notation**: \(\mathbf{x}[k]\in\mathbb{R}^n\) state, \(n\) is the **model order**; \(\mathbf{A}\in\mathbb{R}^{n\times n}\) state matrix; \(\mathbf{C}\in\mathbb{R}^{p\times n}\) output matrix; \(\mathbf{y}[k]\in\mathbb{R}^p\) output (\(p\) sensors); \(\mathbf{w}[k],\mathbf{v}[k]\) process and measurement noise (often white, zero mean); \(i\) is the **number of block rows** (user choice, e.g. 20–50).

### 3.2 Observability matrix

The observability matrix (stack of \(\mathbf{C}\) and its products with powers of \(\mathbf{A}\)) is

$$
\mathcal{O}_i = \begin{bmatrix} \mathbf{C} \\ \mathbf{C}\mathbf{A} \\ \vdots \\ \mathbf{C}\mathbf{A}^{i-1} \end{bmatrix} \in \mathbb{R}^{pi \times n}. \tag{2}
$$

**Meaning**: The 1st block row is \(\mathbf{C}\) (output from current state), the 2nd is \(\mathbf{C}\mathbf{A}\) (output from next-step state), …, the \(i\)-th is \(\mathbf{C}\mathbf{A}^{i-1}\). The **columns** of \(\mathcal{O}_i\) span the “state directions visible in \(i\) output steps”—the observability subspace. **Size**: \(i\) block rows of size \(p\times n\) each, so \(\mathcal{O}_i\) is **\(pi\) rows × \(n\) columns**.

### 3.3 Shift and the upper / lower part of the observability matrix

To recover \(\mathbf{A}\) from \(\mathcal{O}_i\), define two matrices obtained by dropping one block row (both \((i-1)p \times n\)): \(\mathcal{O}_{i,\mathrm{up}}\) (upper part) and \(\mathcal{O}_{i,\mathrm{down}}\) (lower part).

- **\(\mathcal{O}_{i,\mathrm{up}}\)**: \(\mathcal{O}_i\) **without the last block row** (keep block rows 1 to \(i-1\)).
- **\(\mathcal{O}_{i,\mathrm{down}}\)**: \(\mathcal{O}_i\) **without the first block row** (keep block rows 2 to \(i\)).

Then each block row of \(\mathcal{O}_{i,\mathrm{down}}\) equals the corresponding block row of \(\mathcal{O}_{i,\mathrm{up}}\) right-multiplied by \(\mathbf{A}\), i.e. **\(\mathcal{O}_{i,\mathrm{down}} = \mathcal{O}_{i,\mathrm{up}} \mathbf{A}\)**. Hence

$$
\mathbf{A} = \mathcal{O}_{i,\mathrm{up}}^\dagger \mathcal{O}_{i,\mathrm{down}}, \qquad \mathbf{C} = \text{first block row of } \mathcal{O}_i. \tag{3}
$$

\(\dagger\) denotes the Moore–Penrose pseudo-inverse. So: once \(\mathcal{O}_i\) (or a basis for its column space) is obtained from data, \(\mathbf{A}\) and \(\mathbf{C}\) are uniquely recovered by the above.

---

## 4. Algorithms (concrete steps)

The following is written as "what you have, what to do at each step, what you get", so you can implement it directly or check it against code.

### 4.1 SSI-COV (covariance-driven)

**What you have**: Response time series from \(p\) channels (e.g. \(p\) accelerometers), \(N\) samples per channel, sampling interval \(\Delta t\) (e.g. 100 Hz ⇒ \(\Delta t=0.01\) s). Data can be stored as a \(p \times N\) matrix; column \(k\) is the \(p\)-dimensional output \(\mathbf{y}[k]\) at time \(k\).

**Preprocessing (required)**: Subtract the mean of each channel (remove DC); otherwise covariances are biased. Optional: band-pass filter to the band of interest (e.g. 1–20 Hz).

**Step 1 — Compute covariances**  
- Choose block row number \(i\), typically 20–50 (larger uses more memory but gives a richer observability matrix).  
- For lags \(\ell = 0, 1, 2, \ldots, 2i-1\), estimate a \(p \times p\) matrix \(\mathbf{R}_\ell\) with

$$
\mathbf{R}_\ell = \frac{1}{N-\ell} \sum_{k=1}^{N-\ell} \mathbf{y}[k+\ell]\mathbf{y}[k]^T. \tag{4}
$$

- In plain terms: entry \((a,b)\) of \(\mathbf{R}_\ell\) is the (sample) correlation between channel \(a\) and channel \(b\) with lag \(\ell\). You get \(2i\) matrices \(\mathbf{R}_0,\mathbf{R}_1,\ldots,\mathbf{R}_{2i-1}\), each \(p\times p\).

**Step 2 — Build the block Toeplitz matrix**  
- Arrange \(\mathbf{R}_1\) through \(\mathbf{R}_{2i-1}\) into a block Toeplitz matrix \(\mathbf{T}_{1|i}\): the block in row block \(\alpha\), column block \(\beta\) is \(\mathbf{R}_{\alpha+\beta-1}\) (so blocks on the same anti-diagonal are equal). Total size \(pi \times pi\) (e.g. \(p=4,\,i=30\) ⇒ \(120\times 120\)):

$$
\mathbf{T}_{1|i} = \begin{bmatrix} \mathbf{R}_1 & \mathbf{R}_2 & \cdots & \mathbf{R}_i \\ \mathbf{R}_2 & \mathbf{R}_3 & \cdots & \mathbf{R}_{i+1} \\ \vdots & \vdots & \ddots & \vdots \\ \mathbf{R}_i & \mathbf{R}_{i+1} & \cdots & \mathbf{R}_{2i-1} \end{bmatrix}. \tag{5}
$$

- Theory: the row space of \(\mathbf{T}_{1|i}\) equals the column space of the observability matrix, so the next step uses SVD to extract the main directions.

**Step 3 — SVD of T and extract observability**  
- Compute SVD \(\mathbf{T}_{1|i} = \mathbf{U}\mathbf{\Sigma}\mathbf{V}^T\). \(\mathbf{U}\) is \(pi \times pi\), \(\mathbf{\Sigma}\) has singular values in descending order.  
- Choose model order \(n\). Often unknown at first; try several values (e.g. \(n=8,10,12,\ldots,40\)) and use the stabilisation diagram to decide which poles are physical. For now fix some \(n\) (e.g. 20).  
- Keep only the first \(n\) singular values and corresponding left/right singular vectors: \(\mathbf{U}_n\) (\(pi \times n\)), \(\mathbf{\Sigma}_n\) (\(n \times n\)), \(\mathbf{V}_n\) (\(pi \times n\)). Set

$$
\mathcal{O}_i = \mathbf{U}_n \mathbf{\Sigma}_n^{1/2}. \tag{6}
$$

- \(\mathcal{O}_i\) is \(pi \times n\). **Split it into \(i\) blocks of \(p\) rows each**: block 1 (rows 1–\(p\)) is the output matrix \(\mathbf{C}\); blocks 2–\(i\) are related to block 1 by a "shift by \(\mathbf{A}\)", which is used to recover \(\mathbf{A}\).

**Step 4 — Get A and C from O_i**  
- \(\mathbf{C}\): take the **first \(p\) rows and first \(n\) columns** of \(\mathcal{O}_i\), i.e. \(\mathbf{C} = \mathcal{O}_i(1:p,\,1:n)\) (\(p \times n\)).  
- \(\mathbf{A}\):  
  - \(\mathcal{O}_{i,\mathrm{up}}\) = \(\mathcal{O}_i\) **with the last \(p\) rows removed** (keep the first \((i-1)p\) rows), size \((i-1)p \times n\);  
  - \(\mathcal{O}_{i,\mathrm{down}}\) = \(\mathcal{O}_i\) **with the first \(p\) rows removed** (keep the last \((i-1)p\) rows), size \((i-1)p \times n\).  
  - Compute \(\mathbf{A} = \mathcal{O}_{i,\mathrm{up}}^\dagger \, \mathcal{O}_{i,\mathrm{down}}\) (\(\dagger\) = pseudo-inverse) to get \(\mathbf{A}\) of size \(n \times n\).

**Step 5 — Modal parameters from A and C**  
See 4.3 (how to get frequency and damping from eigenvalues of \(\mathbf{A}\), and mode shapes from \(\mathbf{C}\) and eigenvectors).

### 4.2 SSI-DATA (data-driven)

**What you have**: Same as COV—\(p\) channels, \(N\) samples, \(\Delta t\); same preprocessing: remove DC, optional band-pass.

**Step 1 — Build block Hankel and split into "past" and "future"**  
- Choose block row number \(i\) (e.g. 20–50) and number of columns \(j = N - 2i + 1\) (so each column is a segment of length \(2i\)).  
- Build the **block Hankel** matrix: \(2i\) block rows, each of \(p\) rows; block row \(t\) contains data at times \(t,\,t+1,\,\ldots,\,t+j-1\) for all \(p\) channels (stacked into \(p\) rows). Total size \(2pi \times j\).  
- **Past block** \(\mathbf{Y}_p = \mathbf{Y}_{0|i-1}\): first \(i\) block rows (first \(pi\) rows), size \(pi \times j\). Each column = "first \(i\) steps" of output for some starting time.  
- **Future block** \(\mathbf{Y}_f = \mathbf{Y}_{i|2i-1}\): last \(i\) block rows (last \(pi\) rows), size \(pi \times j\). Each column = "next \(i\) steps" for the same starting time.  
- In plain terms: each column is a window "past \(i\) steps + future \(i\) steps"; we use the past to predict the future; the best predictable part spans the observability subspace.

**Projection and oblique projection (in depth)**  
Why does the oblique projection of future onto past yield the observability subspace, and how does it differ from an ordinary (orthogonal) projection?

- **Geometry of past and future**: \(\mathbf{Y}_{0|i-1}\) and \(\mathbf{Y}_{i|2i-1}\) are both \(pi \times j\). Each **column** is a time window of length \(2i\); the past block uses the first \(i\) time steps, the future block the last \(i\). In terms of **rows**, the **row space** of past/future is a subspace of \(\mathbb{R}^j\); in terms of **columns**, the **column space** is a subspace of \(\mathbb{R}^{pi}\). What matters here is column space: the subspace spanned by the columns of the past block, of the future block, and the subspace spanned by the part of the future that is linearly predictable from the past.

- **Orthogonal vs oblique projection**: If we view each **column** of the future block as a vector in \(\mathbb{R}^{pi}\), the **orthogonal projection** onto the column space of the past would give only the part of the future that lies in the past column space. In general past and future column spaces differ, so that projection loses the structure we need. The **oblique projection** instead performs **linear regression**: use the columns of the past to least-squares fit the columns of the future. The result is the matrix that, as prediction of future from past, is optimal in the least-squares sense. That is the oblique projection of future onto past: it projects the future along a direction complementary to the past onto the space spanned by the past, and thereby extracts the future that can be explained by the past.

- **Formula**: Let \(\mathbf{Y}_p = \mathbf{Y}_{0|i-1}\) (past) and \(\mathbf{Y}_f = \mathbf{Y}_{i|2i-1}\) (future). Least-squares prediction of the future from the past (column-wise regression) gives

$$\mathbf{O}_i = \mathbf{Y}_f \mathbf{Y}_p^T (\mathbf{Y}_p \mathbf{Y}_p^T)^{-1} \mathbf{Y}_p.$$

  This is the **oblique projection of future onto past** (\(\mathbf{O}_i\) is \(pi \times j\)). The **column space** of \(\mathbf{O}_i\) is exactly the subspace spanned by the part of the future that is linearly predictable from the past.

- **Why this column space equals the observability subspace**: Under the LTI plus white-noise model, the future output is determined by the current state (evolving via \(\mathbf{A},\mathbf{C}\)) plus future noise. The part of the future that is predictable from the past is exactly the **deterministic evolution** of the state (without future noise), and the map from state to output is given by \(\mathbf{C}\) and powers of \(\mathbf{A}\)—the block rows of the observability matrix \(\mathcal{O}_i\). So the subspace spanned by that predictable part is the **column space** of \(\mathcal{O}_i\). Thus: take the SVD of \(\mathbf{O}_i\); its **left** singular vectors span the observability subspace; keep the first \(n\) and combine with singular values to get an estimate of \(\mathcal{O}_i\) (up to an invertible right factor).

- **Weighting (CVA)**: In practice one often whitens or weights the past and future blocks first (e.g. CVA: left-multiply by \((\mathbf{Y}_p \mathbf{Y}_p^T)^{-1/2}\) and \((\mathbf{Y}_f \mathbf{Y}_f^T)^{-1/2}\)), then take the SVD of the weighted oblique projection and finally un-weight to get \(\mathcal{O}_i\). This improves conditioning and emphasizes the main correlated directions.

**Step 2 — Oblique projection and SVD for observability**  
- Compute \(\mathbf{O}_i = \mathbf{Y}_f \mathbf{Y}_p^T (\mathbf{Y}_p \mathbf{Y}_p^T)^{-1} \mathbf{Y}_p\) (\(pi \times j\)). If using CVA: whiten \(\mathbf{Y}_p,\mathbf{Y}_f\) first, form the weighted oblique projection, take SVD, then un-weight to get \(\mathcal{O}_i\).  
- SVD of \(\mathbf{O}_i\): \(\mathbf{O}_i = \mathbf{U}\mathbf{\Sigma}\mathbf{V}^T\). The **column space** is spanned by the left columns of \(\mathbf{U}\), so take the first \(n\) left singular vectors and singular values and set \(\mathcal{O}_i = \mathbf{U}_n \mathbf{\Sigma}_n^{1/2}\), giving \(pi \times n\) (interpret the rows in blocks of \(p\) as in COV).

**Step 3 — Get A and C from O_i**  
Same as COV: \(\mathbf{C} = \mathcal{O}_i(1:p,\,1:n)\); \(\mathcal{O}_{i,\mathrm{up}}\) = \(\mathcal{O}_i\) with last \(p\) rows removed, \(\mathcal{O}_{i,\mathrm{down}}\) = \(\mathcal{O}_i\) with first \(p\) rows removed, \(\mathbf{A} = \mathcal{O}_{i,\mathrm{up}}^\dagger \mathcal{O}_{i,\mathrm{down}}\).

**Step 4 — Modal parameters**  
See 4.3.

### 4.3 Modal parameter extraction (common to COV and DATA)

Once you have \(\mathbf{A}\) (\(n\times n\)) and \(\mathbf{C}\) (\(p\times n\)), use the following to get frequency, damping and mode shapes (suitable to code directly).

**1. Eigenvalues of A (discrete poles)**  
- Eigen-decompose \(\mathbf{A}\) to get \(n\) eigenvalues \(\mu_1,\ldots,\mu_n\) (complex).  
- **Keep only those inside the unit circle and in conjugate pairs**: physical modes correspond to conjugate pairs with \(|\mu|<1\); those on the unit circle or real are often numerical/model artefacts and can be dropped. With large \(n\) you typically get several conjugate pairs; each pair gives one mode (frequency + damping + mode shape).

**2. Discrete → continuous poles**  
- In the discrete model \(\mu = e^{\lambda \Delta t}\), so \(\lambda = \ln(\mu)/\Delta t\) (principal value for complex \(\mu\)).  
- For each conjugate pair \(\mu_r,\,\mu_r^*\), use only one (e.g. the one with positive imaginary part) to get \(\lambda_r\); the other root is \(\lambda_r^*\).

**3. Frequency and damping ratio**  
- For \(\lambda_r = a_r + \mathrm{i}\,b_r\) (\(a_r\leq 0\)):  
  - Natural angular frequency \(\omega_r = |\lambda_r| = \sqrt{a_r^2 + b_r^2}\);  
  - Natural frequency \(f_r = \omega_r/(2\pi)\) (Hz);  
  - Damping ratio \(\zeta_r = -a_r / |\lambda_r|\) (\(0 \leq \zeta_r \leq 1\) for underdamped).  
- In plain terms: \(|\lambda_r|\) gives frequency; the (negative) real part gives decay rate; damping ratio = |real part| / \(|\lambda_r|\).

**4. Mode shapes**  
- Take the **right eigenvector** \(\mathbf{\psi}_r\) of \(\mathbf{A}\) for \(\mu_r\) (\(n \times 1\)).  
- Mode shape (amplitude and phase at each channel) is \(\mathbf{\phi}_r = \mathbf{C} \mathbf{\psi}_r\) (\(p \times 1\)).  
- **Normalisation**: common choices—(a) normalise \(\mathbf{\phi}_r\) to unit norm (e.g. 2-norm 1), or (b) set one channel (e.g. the first) to 1 and scale the rest. State which convention you use when reporting or plotting.

**Output**: For each physical mode you get three quantities—natural frequency \(f_r\) (Hz), damping ratio \(\zeta_r\), and mode shape \(\mathbf{\phi}_r\) (\(p\)-vector).

### 4.4 Stabilisation diagram and order selection

**Issue**: Model order \(n\) is unknown; too small \(n\) misses modes, too large \(n\) gives many spurious poles. The stabilisation diagram helps separate "true" and "false" poles.

**What to do**  
- Run SSI (COV or DATA) for \(n = n_{\min},\,n_{\min}+2,\,\ldots,\,n_{\max}\) (e.g. \(n=8,10,12,\ldots,40\)). Each run gives a set of poles (frequency + damping).  
- Plot: horizontal axis = order \(n\), vertical axis = frequency (or use colour/size for damping). Each point \((n,\,f)\) is one pole at that order.

**How to read it**  
- **Stable poles**: As \(n\) increases, some frequency–damping points hardly move and form **vertical lines** or **vertical bands**—these are physical modes.  
- **Spurious poles**: Points that jump with \(n\) or have unrealistic damping are usually numerical/overfitting; ignore them.

**Typical stability criteria** (compare two consecutive orders \(n,\,n+2\)):  
- Frequency: \(\Delta f / f < 1\%\) (or stricter 0.5%);  
- Damping: \(\Delta\zeta / \zeta < 5\%\) (or 10%).  
Poles satisfying these are marked "stable"; keep only the \(f_r,\,\zeta_r,\,\mathbf{\phi}_r\) corresponding to stable poles as your final result.

---

## 5. Practice

### 5.1 When to use and limitations

**Use when**: Output-only (ambient) or input–output data; robust, accurate modal identification; multiple or close modes; damping and mode shapes needed; long, multi-channel data available.

**Limitations**: Order \(n\) must be chosen (manually or via stabilisation diagram); block matrices and SVD are costly; method assumes broadband/persistent excitation and approximately white noise; strong noise or non-stationarity can degrade estimates.

### 5.2 Typical parameters and validation

| Parameter | Guideline |
|-----------|-----------|
| Block rows \(i\) | So that \(\mathcal{O}_i\) has full column rank for order \(n\); typically \(i=20\)–50, \(i\cdot p \geq n\) |
| Model order \(n\) | Use stabilisation diagram: try \(n=8,10,12,\ldots\) up to ~2× expected number of modes; keep stable poles |
| Data length \(N\) | \(N \gg 2i\); for civil 1–10 Hz at 100 Hz sampling, 2–10 minutes is common |
| Validation | Compare with FDD/EFDD or ERA; recompute covariances from identified model; check mode-shape consistency |

### 5.3 Edge and online use; comparison with other methods

**Edge/online**: SSI is demanding for low-power edge devices due to large matrices and multi-order stabilisation diagrams; better suited to server/cloud, or lightweight screening on the edge with segments uploaded for SSI.

**Comparison**:  
- **SSI vs FDD**: SSI is time-domain, parametric, gives damping and mode shapes directly, but needs order and stabilisation diagram; FDD is frequency-domain, non-parametric, no damping (EFDD adds it), simpler to implement.  
- **SSI vs ERA**: SSI builds from covariances or raw data via projection+SVD, often preferred in OMA for robustness; ERA builds Hankel from Markov parameters/impulse response and is more compact; when known input is available, ERA can be lighter.  
- **SSI-COV vs SSI-DATA**: COV compresses to covariances first, one Toeplitz+SVD, fixed memory; DATA uses full data and projection, can add CVA-type weighting, often the default in software, higher memory.

---

**Output**: Modal parameters — natural frequencies \(f_r\), damping ratios \(\zeta_r\), mode shapes \(\mathbf{\phi}_r\).
