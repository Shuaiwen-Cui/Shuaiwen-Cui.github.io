# Frequency Domain Decomposition (FDD) / Enhanced FDD (EFDD)

This page expands **2.4 FDD / EFDD** from the [SHM roadmap](../shm.en.md): decompose the output power spectral density matrix via SVD to extract mode shapes and frequencies; EFDD refines damping estimates by single-mode fitting in frequency bands.

---

## Tutorial video

<iframe width="800" height="450" src="https://www.youtube-nocookie.com/embed/placeholder" frameborder="0" allowfullscreen></iframe>

---

## Concept

**Core idea** — For output-only data (ambient excitation), the **power spectral density (PSD) matrix** of the responses contains modal information. At each frequency, perform **singular value decomposition (SVD)** of the PSD matrix; the **dominant singular vector** (first left singular vector) approximates the **mode shape**, and the **first singular value** (largest) peaks at natural frequencies. **EFDD** improves damping estimates by fitting a single-mode model to the singular value curve around each peak. So FDD is: **estimate PSD matrix** → **SVD at each frequency** → **extract mode shapes and frequencies** → **(EFDD) fit damping**.

**Why it works** — Under white or broadband ambient excitation, the output PSD matrix \(\mathbf{S}_{yy}(\omega)\) can be decomposed as

\[\mathbf{S}_{yy}(\omega) = \sum_{r=1}^{N_m} \frac{d_r \boldsymbol{\phi}_r \boldsymbol{\phi}_r^H}{(\omega_r^2 - \omega^2)^2 + (2\zeta_r \omega_r \omega)^2} + \mathbf{S}_n(\omega), \tag{1}\]

**Notation for (1):**

- \(\mathbf{S}_{yy}(\omega) \in \mathbb{C}^{n \times n}\): output PSD matrix at frequency \(\omega\) (\(n\) = number of sensors).
- \(N_m\): number of modes.
- \(\boldsymbol{\phi}_r \in \mathbb{C}^n\): mode shape vector for mode \(r\).
- \(\omega_r, \zeta_r\): undamped natural frequency and damping ratio of mode \(r\).
- \(d_r\): modal participation factor (depends on excitation and mode).
- \(\mathbf{S}_n(\omega)\): noise PSD matrix.
- \((\cdot)^H\): Hermitian transpose.

At a natural frequency \(\omega \approx \omega_r\), the \(r\)-th mode dominates, so the **first singular vector** of \(\mathbf{S}_{yy}(\omega)\) approximates \(\boldsymbol{\phi}_r\), and the **first singular value** is large. FDD exploits this by doing SVD at each frequency line.

**FDD vs EFDD** — **FDD** gives frequencies (from peaks of the first singular value curve) and mode shapes (from the first singular vector at each peak), but **no damping**. **EFDD** adds a step: around each peak, fit a single-degree-of-freedom (SDOF) model to the singular value curve to extract damping ratio \(\zeta_r\).

---

## Detailed algorithm and implementation

### Step 1: Data preparation and PSD estimation

**Input:** Multi-channel output time series \(\{y_i[k]\}\), \(i = 1, \ldots, n\), \(k = 1, \ldots, N\), sampled at rate \(f_s\).

**PSD estimation (Welch method):**

1. **Segment the data:** Divide each channel into overlapping segments (e.g. 50% overlap). Typical segment length: 1024–4096 samples (adjust for desired frequency resolution \(\Delta f = f_s / N_{\text{seg}}\)).

2. **Window each segment:** Apply a window (e.g. Hanning, Hamming) to reduce spectral leakage. For segment \(m\) of channel \(i\):
   \[y_{i,m}^{\text{win}}[k] = w[k] \cdot y_{i,m}[k],\]
   where \(w[k]\) is the window function.

3. **Compute FFT:** For each windowed segment, compute the FFT:
   \[Y_{i,m}(\omega_k) = \text{FFT}\{y_{i,m}^{\text{win}}[k]\},\]
   where \(\omega_k = 2\pi k f_s / N_{\text{seg}}\), \(k = 0, \ldots, N_{\text{seg}}/2\).

4. **Estimate PSD matrix:** At each frequency line \(\omega_k\), form the PSD matrix:
   \[\mathbf{S}_{yy}(\omega_k) = \frac{1}{M} \sum_{m=1}^{M} \mathbf{Y}_m(\omega_k) \mathbf{Y}_m^H(\omega_k), \tag{2}\]
   where \(\mathbf{Y}_m(\omega_k) = [Y_{1,m}(\omega_k), \ldots, Y_{n,m}(\omega_k)]^T\) is the vector of FFT values at frequency \(\omega_k\) for segment \(m\), and \(M\) is the number of segments.

**Practical notes:**

- **Segment length:** Longer segments → better frequency resolution but fewer averages (higher variance). Typical: 2048 samples at \(f_s = 100\) Hz gives \(\Delta f \approx 0.05\) Hz.
- **Overlap:** 50% overlap is common; 75% gives more averages but more computation.
- **Window:** Hanning is standard; use flat-top for amplitude accuracy, rectangular for transient signals.
- **Number of averages:** \(M \geq 20\) is recommended for stable estimates; more averages reduce variance but require longer data.

### Step 2: SVD at each frequency line

**For each frequency \(\omega_k\) in the range of interest:**

1. **Extract PSD matrix:** Get \(\mathbf{S}_{yy}(\omega_k) \in \mathbb{C}^{n \times n}\) (or \(\mathbb{R}^{n \times n}\) if using one-sided PSD).

2. **Perform SVD:**
   \[\mathbf{S}_{yy}(\omega_k) = \mathbf{U}_k \boldsymbol{\Sigma}_k \mathbf{V}_k^H, \tag{3}\]
   where:
   - \(\mathbf{U}_k = [\mathbf{u}_{1,k}, \ldots, \mathbf{u}_{n,k}]\): left singular vectors (columns).
   - \(\boldsymbol{\Sigma}_k = \text{diag}(\sigma_{1,k}, \ldots, \sigma_{n,k})\): singular values (\(\sigma_{1,k} \geq \sigma_{2,k} \geq \cdots \geq \sigma_{n,k} \geq 0\)).
   - \(\mathbf{V}_k\): right singular vectors (not needed for FDD).

3. **Store results:**
   - **First singular value:** \(\sigma_{1,k}\) (used to find frequency peaks).
   - **First singular vector:** \(\mathbf{u}_{1,k}\) (approximates mode shape at \(\omega_k\)).

**Implementation notes:**

- Use **economy SVD** (compute only the first few singular vectors) if only the first mode is of interest.
- For real PSD matrices (one-sided), use real SVD; for complex (two-sided), use complex SVD.
- **Numerical stability:** If \(\mathbf{S}_{yy}(\omega_k)\) is ill-conditioned (e.g. near noise floor), the first singular vector may be unreliable; check condition number or use regularization.

### Step 3: Extract frequencies and mode shapes (FDD)

**Frequency extraction:**

1. **Plot the first singular value curve:** \(\sigma_1(\omega_k)\) vs \(\omega_k\) (or \(f_k = \omega_k / (2\pi)\)).

2. **Identify peaks:** Find local maxima in \(\sigma_1(\omega_k)\). Use peak detection:
   - **Simple:** Find bins where \(\sigma_{1,k} > \sigma_{1,k-1}\) and \(\sigma_{1,k} > \sigma_{1,k+1}\).
   - **Robust:** Use peak prominence (minimum height above neighbors) to reject noise.
   - **Sub-bin accuracy:** Fit a parabola around each peak to get sub-bin frequency estimate.

3. **Record peak frequencies:** For each peak \(r\), record the frequency \(\omega_r\) (or \(f_r\)).

**Mode shape extraction:**

1. **At each peak frequency \(\omega_r\):**
   - Extract the first singular vector: \(\boldsymbol{\phi}_r \approx \mathbf{u}_{1,k_r}\), where \(k_r\) is the frequency bin index of peak \(r\).
   - **Normalize:** Typically normalize to unit norm or set the largest component to 1:
     \[\boldsymbol{\phi}_r \leftarrow \frac{\boldsymbol{\phi}_r}{\|\boldsymbol{\phi}_r\|} \quad \text{or} \quad \boldsymbol{\phi}_r \leftarrow \frac{\boldsymbol{\phi}_r}{\max_i |\phi_{r,i}|}.\]

2. **Handle complex mode shapes:** If \(\boldsymbol{\phi}_r\) is complex (from complex SVD), you can:
   - Use the **magnitude** \(|\boldsymbol{\phi}_r|\) for real mode shapes (common for lightly damped structures).
   - Use the **real part** \(\text{Re}(\boldsymbol{\phi}_r)\) if the structure is real-valued.
   - Keep complex if the structure has complex modes (e.g. gyroscopic effects, non-proportional damping).

**Practical notes:**

- **Peak selection:** Only select peaks that are clearly above the noise floor. A rule of thumb: \(\sigma_{1,\text{peak}} > 2 \cdot \sigma_{1,\text{noise}}\), where \(\sigma_{1,\text{noise}}\) is the average of \(\sigma_1\) in non-resonant regions.
- **Close modes:** If two peaks are very close (within a few frequency bins), the first singular vector may mix both modes. Use the second singular vector or apply EFDD/SSI for better separation.
- **Mode shape consistency:** Check that mode shapes at nearby frequencies (within the same peak region) are similar; large variations indicate noise or mode mixing.

### Step 4: Damping estimation (EFDD)

**EFDD adds damping estimation by fitting a single-mode model around each peak.**

**For each identified peak at frequency \(\omega_r\):**

1. **Select frequency band:** Choose a frequency band around \(\omega_r\): \([\omega_r - \Delta\omega, \omega_r + \Delta\omega]\), where \(\Delta\omega\) is typically 2–5 times the expected bandwidth (e.g. \(\Delta\omega \approx 5 \zeta_r \omega_r\) if damping is known approximately).

2. **Extract singular value curve:** In this band, extract \(\sigma_1(\omega_k)\) for all frequency bins \(k\) in the band.

3. **Fit SDOF model:** Fit the theoretical SDOF response to the singular value curve. The model is:
   \[\sigma_1(\omega) \approx \frac{A}{(\omega_r^2 - \omega^2)^2 + (2\zeta_r \omega_r \omega)^2}, \tag{4}\]
   where \(A\) is an amplitude factor.

4. **Optimization:** Minimize the error between measured \(\sigma_1(\omega_k)\) and model (4) w.r.t. \(\omega_r\) and \(\zeta_r\) (and optionally \(A\)). Use nonlinear optimization (e.g. Levenberg–Marquardt, Gauss–Newton):
   \[\min_{\omega_r, \zeta_r, A} \sum_{k \in \text{band}} \left| \sigma_1(\omega_k) - \frac{A}{(\omega_r^2 - \omega_k^2)^2 + (2\zeta_r \omega_r \omega_k)^2} \right|^2. \tag{5}\]

5. **Extract damping:** The optimized \(\zeta_r\) is the damping ratio estimate.

**Practical notes:**

- **Initial guess:** Use the peak frequency from Step 3 as initial \(\omega_r\); use a typical damping (e.g. \(\zeta = 0.01\) for civil structures, \(0.05\) for mechanical systems) as initial \(\zeta_r\).
- **Bandwidth selection:** Too narrow → insufficient data; too wide → includes neighboring modes. Start with \(\Delta\omega \approx 3 \zeta_r \omega_r\) and adjust.
- **Convergence:** If optimization fails to converge, the mode may be too close to another mode or too lightly damped; try a narrower band or use FDD frequency only (no damping).
- **Validation:** Check that the fitted curve matches the measured \(\sigma_1\) visually; large residuals indicate poor fit (possibly due to mode mixing or noise).

---

## Complete procedure (step-by-step)

**Input:** Multi-channel time series \(\{y_i[k]\}\), \(i = 1, \ldots, n\), \(k = 1, \ldots, N\), sampling rate \(f_s\).

**Step 1: Preprocessing**
- Remove DC offset: \(y_i[k] \leftarrow y_i[k] - \text{mean}(y_i)\).
- (Optional) Detrend: remove linear trend if present.
- (Optional) Filter: apply band-pass filter to focus on frequency range of interest.

**Step 2: PSD matrix estimation**
- Choose segment length \(N_{\text{seg}}\) (e.g. 2048) and overlap (e.g. 50%).
- For each channel \(i\), compute windowed FFTs for all segments.
- At each frequency bin \(k\), form \(\mathbf{S}_{yy}(\omega_k)\) using (2).
- Result: PSD matrices \(\mathbf{S}_{yy}(\omega_k)\) for \(k = 0, \ldots, N_{\text{seg}}/2\).

**Step 3: SVD at each frequency**
- For each \(\omega_k\):
  - Compute SVD: \(\mathbf{S}_{yy}(\omega_k) = \mathbf{U}_k \boldsymbol{\Sigma}_k \mathbf{V}_k^H\).
  - Store \(\sigma_{1,k}\) (first singular value) and \(\mathbf{u}_{1,k}\) (first singular vector).
- Result: \(\sigma_1(\omega_k)\) curve and \(\mathbf{u}_1(\omega_k)\) vectors.

**Step 4: Frequency identification (FDD)**
- Plot \(\sigma_1(\omega_k)\) vs \(f_k = \omega_k / (2\pi)\).
- Identify peaks using peak detection (with prominence threshold).
- For each peak \(r\), record frequency \(f_r\) (with sub-bin accuracy if fitted).

**Step 5: Mode shape extraction (FDD)**
- For each peak \(r\) at frequency bin \(k_r\):
  - Extract mode shape: \(\boldsymbol{\phi}_r = \mathbf{u}_{1,k_r}\).
  - Normalize (unit norm or max component = 1).
  - (Optional) Check consistency with nearby bins.
- Result: Mode shapes \(\boldsymbol{\phi}_r\) for \(r = 1, \ldots, N_m\).

**Step 6: Damping estimation (EFDD, optional)**
- For each peak \(r\):
  - Select frequency band around \(f_r\).
  - Extract \(\sigma_1(\omega_k)\) in the band.
  - Fit SDOF model (4) using optimization (5).
  - Extract damping ratio \(\zeta_r\).
- Result: Damping ratios \(\zeta_r\) for \(r = 1, \ldots, N_m\).

**Output:** Modal parameters: frequencies \(f_r\), mode shapes \(\boldsymbol{\phi}_r\), damping ratios \(\zeta_r\) (if EFDD).

---

## When to use and limitations

**Use when:**
- **Output-only data** (ambient excitation, operational conditions).
- **Multiple sensors** (at least 2, preferably 4+ for reliable mode shapes).
- **Broadband or white excitation** (ensures all modes are excited).
- **Well-separated modes** (FDD works best when modes are not too close).
- **Moderate damping** (very light damping: peaks are sharp, EFDD fitting may be difficult; very heavy damping: peaks are broad, hard to separate).

**Limitations:**
- **No damping from FDD:** FDD gives frequencies and mode shapes only; need EFDD for damping.
- **Close modes:** When modes are very close (within a few frequency bins), the first singular vector may mix modes; use second/third singular vectors or switch to SSI.
- **Mode shape accuracy:** Singular vectors are approximations; accuracy depends on SNR, number of sensors, and mode separation.
- **Excitation assumption:** Assumes white or broadband excitation; narrowband excitation may bias results.
- **Computational cost:** SVD at each frequency line; for many frequencies and sensors, cost grows as \(O(n^3 \cdot N_f)\), where \(N_f\) is the number of frequency bins.

---

## Engineering practice: detailed guidelines

### PSD estimation parameters

| Parameter | Typical values | Notes |
|-----------|----------------|-------|
| Segment length \(N_{\text{seg}}\) | 1024–4096 samples | Longer → better frequency resolution; shorter → more averages. Choose based on \(f_s\) and desired \(\Delta f\). |
| Overlap | 50% (standard), 75% (more averages) | More overlap → more averages but more computation. |
| Window | Hanning (standard), Hamming, flat-top | Hanning: good balance; flat-top: amplitude accuracy; rectangular: transient signals. |
| Number of averages \(M\) | \(\geq 20\) (recommended) | More averages → lower variance but need longer data. |
| Frequency resolution \(\Delta f\) | \(f_s / N_{\text{seg}}\) | Typical: 0.05–0.1 Hz for civil structures (e.g. \(f_s = 100\) Hz, \(N_{\text{seg}} = 2048\) → \(\Delta f \approx 0.05\) Hz). |

### Peak detection and mode selection

| Aspect | Guidelines |
|--------|------------|
| Peak prominence | Set threshold: \(\sigma_{1,\text{peak}} > 2 \cdot \sigma_{1,\text{noise}}\) (or use automatic prominence). |
| Sub-bin frequency | Fit parabola around peak: \(f_{\text{peak}} = f_k + \delta\), where \(\delta\) is from parabola vertex. |
| Close modes | If peaks within 3–5 bins, check second singular vector; if still mixed, use SSI. |
| Spurious peaks | Check mode shape consistency; spurious peaks have inconsistent or noisy mode shapes. |
| Frequency range | Focus on range where modes are expected (e.g. 0.1–10 Hz for bridges, 1–50 Hz for buildings). |

### EFDD damping fitting

| Aspect | Guidelines |
|--------|------------|
| Frequency band width | Start with \(\Delta\omega \approx 3 \zeta_r \omega_r\); adjust if fit is poor (narrower if close modes, wider if noisy). |
| Initial guess | \(\omega_r\): from FDD peak; \(\zeta_r\): 0.01 (civil), 0.05 (mechanical), or from prior knowledge. |
| Optimization method | Levenberg–Marquardt is robust; Gauss–Newton is faster if near solution. |
| Convergence | If fails, try narrower band, different initial guess, or skip damping (use FDD frequency only). |
| Validation | Plot fitted curve vs measured \(\sigma_1\); large residuals indicate poor fit (mode mixing or noise). |

### Common issues and solutions

| Issue | Symptoms | Solutions |
|-------|----------|----------|
| Noisy singular value curve | Many small peaks, unstable mode shapes | Increase number of averages \(M\); use longer segments; apply smoothing to \(\sigma_1\) curve. |
| Mode mixing | Mode shapes change rapidly near peak | Check second/third singular vectors; use narrower frequency band for EFDD; switch to SSI. |
| Missing modes | Expected mode not found in \(\sigma_1\) curve | Check excitation covers that frequency; increase frequency resolution; check sensor placement (avoid nodes). |
| Poor damping fit | EFDD optimization fails or gives unrealistic \(\zeta\) | Try different initial guess; adjust frequency band; check if mode is too close to another; use FDD only. |
| Inconsistent mode shapes | Mode shape varies across nearby frequencies | Check SNR; increase averages; verify sensor synchronization; check for nonlinearity. |

---

## Edge and online computing

**Suitability** — **Moderately suited.** Output-only (no input measurement needed); SVD per frequency line can be done in batches; manageable when channels and lines are limited. However, PSD estimation needs buffer, and EFDD damping fit adds computation.

**Implementation strategy:**

1. **PSD estimation:** Use **sliding-window Welch** or **block processing**:
   - Buffer \(N_{\text{seg}}\) samples per channel.
   - Compute windowed FFT for current block.
   - Update PSD estimate (exponential moving average or block average).
   - Trade-off: more averages → better quality but more latency.

2. **SVD computation:** For each frequency bin:
   - **Economy SVD:** Compute only first 1–2 singular vectors (most modes are in the first).
   - **Batch processing:** Process multiple frequency bins in parallel if hardware allows.
   - **Incremental updates:** For slowly varying systems, update SVD incrementally (e.g. rank-1 updates).

3. **Peak detection:** Lightweight (comparisons and local fits); can run in real time.

4. **EFDD damping:** More expensive (nonlinear optimization); options:
   - Run EFDD **offline** or in **cloud** (upload peak frequencies, download damping).
   - Use **pre-computed damping** from baseline (assume damping doesn't change much).
   - **Skip EFDD** on edge, use FDD frequencies only (damping from other methods or prior).

**Practical edge deployment:**

- **Tiered pipeline:**
  1. **Edge:** PSD estimation + SVD + peak detection → frequencies and mode shapes (FDD).
  2. **Cloud/offline:** EFDD damping fitting for identified peaks.
- **Resource limits:**
  - Limit number of frequency bins (focus on expected modal range).
  - Limit number of sensors (use subset if needed).
  - Use economy SVD (first 1–2 modes only).
- **Real-time screening:** Compare FDD frequencies with baseline; flag significant shifts; upload suspect segments for full analysis.

**Challenges:**
- **Memory:** PSD matrices and SVD results need storage; limit frequency range and sensors.
- **Computation:** SVD at each frequency; for many bins, cost is significant on MCUs.
- **Latency:** PSD needs buffer; real-time updates have delay (e.g. \(N_{\text{seg}} / f_s\) seconds).
- **Quality:** Edge SNR may be lower; fewer averages → noisier results; use smoothing or longer buffers.

---

## Comparison with other methods

**FDD vs PP:**
- **FDD:** Uses SVD of PSD matrix → more robust to noise, better for close modes, gives mode shapes directly.
- **PP:** Simpler (just peak picking on spectra) → faster, less robust, mode shapes from amplitude ratios.

**FDD vs SSI:**
- **FDD:** Frequency-domain, non-parametric → simpler, no order selection, but no damping (need EFDD).
- **SSI:** Time-domain, parametric → more accurate, gives damping, but needs order selection and more computation.

**FDD vs EFDD:**
- **FDD:** Fast, gives frequencies and mode shapes, no damping.
- **EFDD:** Adds damping estimation via SDOF fitting → more complete but slower.

**When to choose:**
- Use **FDD** for quick screening, well-separated modes, when damping is not critical.
- Use **EFDD** when damping is needed and modes are reasonably separated.
- Use **SSI** for highest accuracy, close modes, or when full modal model is needed.
