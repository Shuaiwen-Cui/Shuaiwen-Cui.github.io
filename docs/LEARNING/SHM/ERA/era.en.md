# Eigensystem Realization Algorithm (ERA) / NExT-ERA

This page expands **2.7 ERA / NExT-ERA** from the [SHM roadmap](../shm.en.md): **ERA** builds **Markov parameters** and a **Hankel matrix** from impulse or free-decay response, then uses **SVD** and minimal realisation to obtain a state-space model and modal parameters; **NExT-ERA** under ambient excitation builds equivalent impulse responses from output auto- and cross-correlations, then applies ERA for output-only operational modal analysis.

---

## Tutorial video

<iframe width="800" height="450" src="https://www.youtube-nocookie.com/embed/5TM_JrW5qDo" frameborder="0" allowfullscreen></iframe>。

---

## Concept

**Core idea** — **ERA** (Eigensystem Realization Algorithm) is a time-domain identification method: from impulse or free-decay sequences we form Markov parameters, assemble a Hankel matrix, perform SVD and truncation to get a minimal realisation, yielding a discrete state-space \((\mathbf{A},\mathbf{B},\mathbf{C})\); the eigenvalues of \(\mathbf{A}\) give continuous-time poles (natural frequency and damping ratio), and \(\mathbf{C}\) with eigenvectors give mode shapes. **NExT** (Natural Excitation Technique) states that under white or broadband ambient excitation, the cross-correlation of two response channels satisfies the same homogeneous differential equation as the impulse response, so we can estimate auto/cross-correlations from output-only data and use them as “equivalent impulse response” in ERA, i.e. **NExT-ERA**. In short: ERA is impulse/free decay → Markov → Hankel → SVD → state-space → modes; NExT-ERA is output-only → correlation estimate → equivalent impulse → ERA.

**Principle (brief)** — Discrete LTI state-space:

\[\mathbf{x}[k+1] = \mathbf{A} \mathbf{x}[k] + \mathbf{B} u[k], \qquad \mathbf{y}[k] = \mathbf{C} \mathbf{x}[k] + \mathbf{D} u[k]. \tag{1}\]

The pulse response (Markov parameters) are \(\mathbf{h}_0 = \mathbf{D}\), \(\mathbf{h}_k = \mathbf{C} \mathbf{A}^{k-1} \mathbf{B}\) (\(k \geq 1\)). Stacking \(\mathbf{h}_k\) into a Hankel matrix \(\mathbf{H}\) links its rank to observability/controllability; SVD of \(\mathbf{H}\) truncated to order \(n\) yields \(\mathbf{U},\mathbf{\Sigma},\mathbf{V}\) from which \(\mathbf{A},\mathbf{B},\mathbf{C}\) are recovered (minimal realisation). Eigenvalues \(\mu_i\) of \(\mathbf{A}\) are discrete poles; map to continuous time by \(\lambda_i = \ln(\mu_i)/\Delta t\) to get modal poles, hence \(\omega_r\), \(\zeta_r\) and mode shapes \(\mathbf{\phi}_r\).

**Notation for (1):**

- \(\mathbf{x}[k] \in \mathbb{R}^n\): state vector; \(n\) is the state dimension (model order).

- \(\mathbf{A} \in \mathbb{R}^{n \times n}\): state matrix; its eigenvalues determine modal frequency and damping.

- \(\mathbf{B},\mathbf{C},\mathbf{D}\): input, output, and direct matrices; pulse response is \(\mathbf{h}_k = \mathbf{C} \mathbf{A}^{k-1} \mathbf{B}\).

- \(u[k], \mathbf{y}[k]\): input and output (vector for multi-channel).

**ERA vs NExT-ERA** — **ERA** requires impulse or free-decay response (e.g. from impact tests or known input); **NExT-ERA** uses only response under ambient excitation, replacing the impulse response with correlation functions before running the same Hankel + SVD + realisation, so it is used for output-only operational modal analysis (e.g. bridges, buildings).

---

## Detailed algorithm and implementation

### Step 1: Data preparation and “pulse / equivalent pulse” sequence

**ERA (with impulse or free decay):**

1.1  **Input:** Impulse response sequence \(\mathbf{h}_0, \mathbf{h}_1, \ldots, \mathbf{h}_L\) (matrix sequence for multi-channel), or equivalent Markov parameters from free decay; sampling interval \(\Delta t\). In practice, \(L\) should cover 2–3 times the decay time of the modes of interest (e.g. for 1% damping, decay time \(\sim 1/(\zeta\omega_n)\)), so take \(L \approx \lfloor (2\sim3) \times (1/(\zeta\omega_n)) \times f_s \rfloor\).

1.2  **Preprocessing:** Remove DC (subtract mean per channel); optional detrend and band-pass filter (keep the band of interest); if using free decay, apply RDT or similar first for cleaner decay, then convert to Markov form.

**NExT-ERA (output-only):**

1.3  **Input:** Multi-channel response time series \(\{y_i[k]\}\), \(i=1,\ldots,n\), \(k=1,\ldots,N\), sampling rate \(f_s\), \(\Delta t = 1/f_s\). Use enough samples \(N\) for stable correlation estimates (e.g. at least tens of thousands, or 2–5 minutes of data).

1.4  **Estimate correlation:** Choose a reference channel (e.g. \(y_{\mathrm{ref}}\), preferably high SNR and not at a node of the modes of interest). For each channel \(i\) estimate cross-correlation with reference (or auto-correlation):

\[R_{i,\mathrm{ref}}(\tau) = \frac{1}{N-\tau} \sum_{k=1}^{N-\tau} y_i[k] \, y_{\mathrm{ref}}[k+\tau], \qquad \tau = 0,1,\ldots,L. \tag{2}\]

Under white-noise excitation, \(R_{i,\mathrm{ref}}(\tau)\) is proportional to the “impulse response from reference to channel \(i\)”, so \(R_{i,\mathrm{ref}}(0), R_{i,\mathrm{ref}}(1), \ldots\) can be used as equivalent Markov parameters in the following ERA steps.

**Reproducible detail (NExT-ERA):** For single reference and multiple outputs, the \(k\)-th Markov block can be the column vector \(\mathbf{Y}_k = [R_{1,\mathrm{ref}}(k), \ldots, R_{n,\mathrm{ref}}(k)]^T\) (\(n\) = number of channels); for matrix form, treat \(\mathbf{Y}_k\) as an \(n\times 1\) block so each block in the Hankel is this column.

**Practical notes:**

- Reference channel: choose one with good SNR and not at a node of the modes of interest; try several references and compare **stabilisation diagrams**.
- Correlation length \(L\): cover the decay time of the modes of interest; typically hundreds to a few thousand samples (e.g. \(f_s=100\) Hz, 1% damping, lowest mode ~1 Hz → decay time ~15 s, \(L \approx 1500\)).
- NExT assumption: excitation approximately white or broadband; coloured excitation can bias correlation and modal estimates.

### Step 2: Build the Hankel matrix

2.1  **Block Markov parameters:** Let \(\mathbf{Y}_k\) be the \(k\)-th Markov parameter block (a matrix or column vector for multi-channel). Choose block row count \(r\) and block column count \(s\) (large enough to capture the dynamics; typically \(r=s=20\sim50\), and \(r+s-1 \le L+1\) so that \(\mathbf{Y}_{r+s-1}\) is available), and form:

\[\mathbf{H}(0) = \begin{bmatrix} \mathbf{Y}_1 & \mathbf{Y}_2 & \cdots & \mathbf{Y}_s \\ \mathbf{Y}_2 & \mathbf{Y}_3 & \cdots & \mathbf{Y}_{s+1} \\ \vdots & \vdots & \ddots & \vdots \\ \mathbf{Y}_r & \mathbf{Y}_{r+1} & \cdots & \mathbf{Y}_{r+s-1} \end{bmatrix}, \qquad \mathbf{H}(1) = \begin{bmatrix} \mathbf{Y}_2 & \mathbf{Y}_3 & \cdots & \mathbf{Y}_{s+1} \\ \mathbf{Y}_3 & \mathbf{Y}_4 & \cdots & \mathbf{Y}_{s+2} \\ \vdots & \vdots & \ddots & \vdots \\ \mathbf{Y}_{r+1} & \mathbf{Y}_{r+2} & \cdots & \mathbf{Y}_{r+s} \end{bmatrix}. \tag{3}\]

\(\mathbf{H}(1)\) is \(\mathbf{H}(0)\) shifted down by one block row.

2.2  **Scalar case:** If the pulse response is scalar \(h_k\), then \(\mathbf{Y}_k = h_k\) and \(\mathbf{H}(0),\mathbf{H}(1)\) are standard Hankel matrices.

### Step 3: SVD and minimal realisation

3.1  **SVD of \(\mathbf{H}(0)\):** \(\mathbf{H}(0) = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^T\) (or \(\mathbf{V}^H\)), \(\mathbf{\Sigma} = \mathrm{diag}(\sigma_1,\ldots,\sigma_{\min})\).

3.2  **Truncate to order \(n\):** From the singular value drop or stabilisation diagram, keep the first \(n\) singular values (\(n\) on the order of 2× number of modes), giving \(\mathbf{U}_n, \mathbf{\Sigma}_n, \mathbf{V}_n\).

3.3  **Recover state-space:** Define \(\mathbf{O} = \mathbf{U}_n \mathbf{\Sigma}_n^{1/2}\) (observability block), \(\mathbf{C}_c = \mathbf{\Sigma}_n^{1/2} \mathbf{V}_n^T\) (controllability block). Then

\[\mathbf{A} = \mathbf{\Sigma}_n^{-1/2} \mathbf{U}_n^T \mathbf{H}(1) \mathbf{V}_n \mathbf{\Sigma}_n^{-1/2}, \tag{4}\]

\(\mathbf{B}\) is \(\mathbf{\Sigma}_n^{-1/2}\mathbf{U}_n^T\) times the first block column of \(\mathbf{H}(0)\); \(\mathbf{C}\) is the first block row of \(\mathbf{H}(0)\) times \(\mathbf{V}_n\mathbf{\Sigma}_n^{-1/2}\) (standard ERA; ensure block dimensions match in code).

**Implementation notes:**

- Order selection: try several \(n\) (e.g. \(n=8,10,12,\ldots,\) up to about 2× expected number of modes) and build a **stabilisation diagram** (poles vs order); stable poles correspond to physical modes, drifting or scattered ones are often spurious. Common stability criteria: relative frequency change \(<\Delta f/f\) (e.g. 1%) and relative damping change \(<\Delta\zeta/\zeta\) (e.g. 5%).
- Numerics: \(\mathbf{H}(0)\) can be ill-conditioned; SVD truncation acts as regularisation; if needed use a slightly larger \(n\) and filter modes via the stabilisation diagram; avoid very small singular values in \(\mathbf{\Sigma}_n\) (set a minimum threshold or drop them).

### Step 4: Extract modal parameters from state-space

4.1  **Discrete eigenvalues:** Compute eigenvalues \(\mu_i\) of \(\mathbf{A}\), \(i=1,\ldots,n\). Drop poles outside the unit circle or clearly due to noise.

4.2  **Continuous poles:** \(\lambda_i = \ln(\mu_i) / \Delta t\) (principal value). For a conjugate pair \(\lambda_{r},\lambda_{r}^*\):

\[\omega_r = |\lambda_r|, \qquad \zeta_r = -\mathrm{Re}(\lambda_r) / |\lambda_r|. \tag{5}\]

Natural frequency \(f_r = \omega_r / (2\pi)\), damping ratio \(\zeta_r\).

4.3  **Mode shapes:** From \(\mathbf{C}\) and the right eigenvectors \(\mathbf{\psi}_r\) of \(\mathbf{A}\), the output-space mode shape is \(\mathbf{\phi}_r \propto \mathbf{C} \mathbf{\psi}_r\); normalise by channel (e.g. unit norm or max component 1).

**Practical notes:**

- Conjugate pairs: physical modes correspond to conjugate pole pairs; single real poles are often numerical or noise and can be dropped.
- Stabilisation diagram: as \(n\) increases, true modal frequency/damping should stabilise; spurious modes tend to drift or jump; keep only poles that are stable over several consecutive orders.

---

## Complete procedure (step-by-step)

**Input (ERA):** Impulse or free-decay sequence; or **(NExT-ERA)** multi-channel response \(\{y_i[k]\}\), sampling rate \(f_s\).

**Step 1: Preprocessing and “pulse” sequence**

1.1 For ERA: Have impulse/free-decay sequence and \(\Delta t\); optional remove DC and filter.

1.2 For NExT-ERA: Choose reference channel; estimate cross-correlation \(R_{i,\mathrm{ref}}(\tau)\) from (2), \(\tau=0,\ldots,L\); treat \(R\) as equivalent Markov sequence.

**Step 2: Build Hankel**

2.1 Block the Markov parameters (or equivalent) to form \(\mathbf{H}(0)\) and \(\mathbf{H}(1)\).

**Step 3: SVD and realisation**

3.1 SVD of \(\mathbf{H}(0)\); choose order \(n\) from singular values or stabilisation diagram; truncate to \(\mathbf{U}_n, \mathbf{\Sigma}_n, \mathbf{V}_n\).

3.2 Recover \(\mathbf{A},\mathbf{B},\mathbf{C}\) from (4) and standard ERA formulas.

**Step 4: Modal extraction**

4.1 Eigenvalues \(\mu_i\) of \(\mathbf{A}\); convert to continuous poles \(\lambda_i\), then \(\omega_r, \zeta_r\) from (5).

4.2 Mode shapes \(\mathbf{\phi}_r\) from \(\mathbf{C}\) and eigenvectors of \(\mathbf{A}\); normalise.

**Output:** Modal parameters — natural frequencies \(f_r\), damping ratios \(\zeta_r\), mode shapes \(\mathbf{\phi}_r\).

**Reproducible checklist (NExT-ERA example):**

1. Prepare: \(n_y\)-channel response \(y_i[k]\), \(k=1,\ldots,N\), sampling rate \(f_s\); remove DC (subtract mean per channel); optional band-pass filter.
2. Choose reference channel \(y_{\mathrm{ref}}\); set \(L = \lfloor 2.5 \cdot f_s / (\zeta \omega_n) \rfloor\) (\(\zeta\approx 0.01\), \(\omega_n\) = rough centre frequency in rad/s of the band of interest).
3. Compute \(R_{i,\mathrm{ref}}(\tau)\) from (2), \(\tau=0,\ldots,L\); form \(\mathbf{Y}_k = [R_{1,\mathrm{ref}}(k),\ldots,R_{n_y,\mathrm{ref}}(k)]^T\).
4. Set \(r=s=30\) (or 20–50), build \(\mathbf{H}(0),\mathbf{H}(1)\); run SVD on \(\mathbf{H}(0)\).
5. Try \(n=8,10,12,\ldots,24\); for each \(n\) get \(\mathbf{A}\) from (4), compute eigenvalues and convert to \(f_r,\zeta_r\); plot stabilisation diagram and keep poles whose frequency/damping stabilise with \(n\) (e.g. \(\Delta f/f<1\%\)).
6. Get mode shapes from (5) and \(\mathbf{C}\), eigenvectors; normalise; optionally validate by recomputing correlation from the identified model and comparing with original \(R\).

---

## When to use and limitations

**Use when:**

- **ERA:** Impulse or free-decay response is available (impact test, step relaxation, etc.); time-domain, single SVD, no iteration is desired.
- **NExT-ERA:** Output-only, ambient excitation (bridges, buildings, large structures); works best when excitation is approximately broadband/white.
- Multi-channel: multiple sensors give mode shapes; reference channel choice affects quality.

**Limitations:**

- Order sensitivity: too small \(n\) misses modes; too large introduces spurious poles; use stabilisation diagram or multiple trial orders.
- NExT assumption: white/broadband excitation; narrowband or coloured excitation can bias correlation and modal estimates.
- Correlation estimation: requires sufficient data and lag length; long correlation increases Hankel size and compute.
- Close modes: as with FDD/SSI, separating very close modes is harder.

---

## Engineering practice and reproducible parameters

| Aspect | Guideline and typical values |
|--------|-------------------------------|
| Reference channel (NExT-ERA) | Choose high SNR and avoid nodes of modes of interest; try several references and compare stabilisation diagrams. |
| Correlation length \(L\) | Cover 2–3 times the decay time of the modes of interest; e.g. \(L = \lfloor 2.5/(\zeta\omega_n) \cdot f_s \rfloor\) with \(\zeta\approx 0.01\); too long makes Hankel large and costly. |
| Block rows/columns \(r,s\) | Take \(r,s\) large enough so Hankel rank reflects system order; typically \(r=s=20\sim50\), with \(r+s-1\le L+1\). |
| Order \(n\) | Use stabilisation diagram: try \(n=8,10,12,\ldots\) up to ~2× expected number of modes; keep poles whose frequency/damping stabilise (e.g. \(\Delta f/f<1\%\), \(\Delta\zeta/\zeta<5\%\)). |
| Numerical stability | SVD truncation regularises; drop very small singular values (e.g. \(\sigma_i/\sigma_1 < 10^{-6}\)) to avoid ill-conditioned \(\mathbf{A}\). |
| Validation | Recompute impulse response or correlation from identified \(\mathbf{A},\mathbf{B},\mathbf{C}\) and compare with raw data; or cross-check with FDD/SSI. |

---

## Edge and online computing

**Suitability** — Moderately suited. Output-only ( **NExT-ERA** ) needs no input measurement; one Hankel build + one SVD, no iteration; compute can be controlled by limiting Hankel size and order \(n\).

**Implementation strategy:**

1. Correlation estimation: sliding window or block-wise correlation online; needs buffer and multiply-add; keep length \(L\) moderate.
2. Hankel and SVD: limit \(r,s\) and \(n\); economy SVD (first few singular values only) if acceptable.
3. Stabilisation diagram: on the edge, compute only a few orders for a simple stability check; full diagram can be done in the cloud or offline.

**Challenges:**

- Memory: Hankel and SVD scale with \(r,s,n\); limits needed for MCUs.
- Order selection: full stabilisation diagram is hard on the edge; often rely on prior or fixed \(n\).
- Correlation estimation: long data and long lags need buffer and compute.

---

## Comparison with other methods

**ERA vs FDD:**

- **ERA:** Time-domain, parametric, single SVD + eigenvalue solve; gives frequency, damping, mode shapes directly; requires order selection.
- **FDD:** Frequency-domain, non-parametric, SVD per frequency line; no damping (need EFDD); no order choice, but close-mode separation relies on singular value curves.

**ERA vs SSI:**

- **ERA:** Hankel from Markov/correlation; compact formulation; single SVD, relatively low compute.
- **SSI:** Hankel/Toeplitz built directly from data, projection + SVD; generally more robust and widely used in practice, but larger matrices and more compute.

**NExT-ERA vs other output-only methods:**

- **NExT-ERA:** Correlation → equivalent impulse → ERA; clear structure; fits well with RDT or similar preprocessing feeding ERA.
- **FDD/EFDD:** Frequency-domain PSD + SVD; no order; EFDD adds damping.
- **SSI:** Direct from data, robust; higher compute and memory.

**When to choose:**

- Prefer **ERA** when impulse or free-decay data are available.
- For output-only, ambient excitation with a simple time-domain, single-SVD pipeline, use **NExT-ERA**; for maximum robustness and full stabilisation diagram, consider SSI.
