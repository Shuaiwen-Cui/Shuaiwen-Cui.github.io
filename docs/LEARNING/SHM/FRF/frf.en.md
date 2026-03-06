# FRF curve fitting

This page expands **2.3 FRF** from the [SHM roadmap](../shm.en.md): fit measured Frequency Response Functions (FRFs) in the frequency domain with rational functions (or orthogonal polynomials), and extract modal frequency, damping, and mode shape from poles and residues.

---

## Tutorial video

<iframe width="800" height="450" src="https://www.youtube-nocookie.com/embed/5wpFgU_oAss" frameborder="0" allowfullscreen></iframe>

---

## Concept

**Core idea** — The Frequency Response Function \(H(\omega)\) is the ratio of output to input in the frequency domain; it is a **rational function** (ratio of polynomials) whose **poles** correspond to the system's natural frequencies and damping, and whose **residues** (at each pole) give the mode shapes. By **fitting a rational function** to measured FRF data, we can extract all modal parameters: frequency, damping, and mode shape. So FRF curve fitting is: **measure FRF** → **fit rational function** → **extract poles and residues** → **modal parameters**.

**FRF definition** — For a linear time-invariant system, the FRF \(H_{ij}(\omega)\) from input \(j\) to output \(i\) is

\[H_{ij}(\omega) = \frac{Y_i(\omega)}{U_j(\omega)}, \tag{1}\]

**Notation for (1):**

- \(H_{ij}(\omega)\): FRF from input \(j\) to output \(i\) at frequency \(\omega\).
- \(Y_i(\omega), U_j(\omega)\): Fourier transforms of output \(i\) and input \(j\).
- \(\omega\): frequency (rad/s).

In modal testing, inputs are typically forces (impact hammer, shaker) and outputs are responses (acceleration, velocity, displacement).

**Rational function form** — The FRF can be written as a rational function (ratio of polynomials). The **Rational Fraction Polynomial (RFP)** form is common:

\[H(\omega) = \frac{\sum_{k=0}^{m} a_k (j\omega)^k}{\sum_{k=0}^{n} b_k (j\omega)^k} = \frac{N(\omega)}{D(\omega)}, \tag{2}\]

**Notation for (2):**

- \(a_k, b_k\): numerator and denominator polynomial coefficients (real or complex).
- \(m, n\): orders of numerator and denominator polynomials.
- \(j\): imaginary unit.
- \(N(\omega), D(\omega)\): numerator and denominator polynomials.

Alternatively, in **partial fraction form** (pole–residue form):

\[H(\omega) = \sum_{r=1}^{N_m} \frac{R_r}{j\omega - s_r} + \text{(higher-order terms)}, \tag{3}\]

**Notation for (3):**

- \(N_m\): number of modes in the frequency band.
- \(s_r\): **pole** for mode \(r\) (complex, \(s_r = -\zeta_r \omega_r \pm j\omega_r\sqrt{1-\zeta_r^2}\)).
- \(R_r\): **residue** for mode \(r\) (complex, related to mode shape).
- \(\omega_r, \zeta_r\): undamped natural frequency and damping ratio of mode \(r\).

**Modal parameter extraction** — From the fitted rational function:

1. **Poles** \(s_r\): solve \(D(\omega) = 0\) (or extract from partial fraction form) → get \(\omega_r\) and \(\zeta_r\).
2. **Residues** \(R_r\): from the partial fraction expansion → mode shapes (for MIMO, residues form a matrix whose columns are mode shapes).

---

## Algorithm in brief

**FRF estimation** — First, estimate the FRF from measured input–output data. Common methods:

- **H1 estimator:** \(H_1(\omega) = S_{yu}(\omega) / S_{uu}(\omega)\), where \(S_{yu}\) is the cross-power spectral density (CPSD) between output and input, and \(S_{uu}\) is the auto-PSD of input. Assumes output noise only.
- **H2 estimator:** \(H_2(\omega) = S_{yy}(\omega) / S_{uy}(\omega)\), assumes input noise only.
- **Hv estimator:** Uses SVD of the CPSD matrix for multiple inputs/outputs.

**Rational function fitting** — Fit the measured FRF \(\hat{H}(\omega_k)\) at frequencies \(\omega_k\) with a rational function. The problem is **nonlinear** in the polynomial coefficients (or poles/residues). Common approaches:

1. **Linear least squares (LS) on numerator/denominator:** Rewrite (2) as \(D(\omega) H(\omega) = N(\omega)\), linearize in \(a_k, b_k\), solve LS. Iterate if needed (e.g. Sanathanan–Koerner iteration).
2. **Nonlinear optimization:** Minimize \(\sum_k |H(\omega_k) - \hat{H}(\omega_k)|^2\) w.r.t. poles and residues (or polynomial coefficients) using Gauss–Newton, Levenberg–Marquardt, etc.
3. **Orthogonal polynomials:** Use orthogonal basis (e.g. Forsythe polynomials) to improve numerical conditioning.

**Pole extraction** — From the fitted denominator \(D(\omega)\) (or from partial fraction form), find roots \(s_r\) such that \(D(s_r) = 0\). Map to continuous time if needed (e.g. \(s = j\omega\) for frequency-domain fitting).

**Residue and mode shape extraction** — For each pole \(s_r\), compute the residue \(R_r\) (from partial fraction expansion or from the fitted model). For MIMO systems, residues form a matrix; the columns (or rows) give mode shapes.

---

## Procedure (outline)

1. **Data acquisition:** Measure input (force) and output (response) signals simultaneously; ensure known, controlled excitation (impact hammer, shaker).
2. **FRF estimation:** Compute FRF from input–output data using H1, H2, or Hv estimator; average over multiple measurements if available (e.g. multiple impacts, repeated shaker sweeps).
3. **Frequency band selection:** Choose the frequency band of interest; may fit multiple bands separately or use a global fit.
4. **Model order selection:** Choose numerator and denominator orders (\(m, n\)) or number of modes \(N_m\); use stabilisation diagram (poles vs order) or information criteria (AIC, BIC).
5. **Rational function fitting:** Fit the measured FRF with a rational function using LS, nonlinear optimization, or orthogonal polynomials; iterate if needed (e.g. Sanathanan–Koerner).
6. **Pole extraction:** Find roots of the denominator polynomial → discrete or continuous poles \(s_r\).
7. **Modal parameter extraction:** From poles \(s_r\), extract \(\omega_r\) and \(\zeta_r\); from residues \(R_r\), extract mode shapes (for MIMO).
8. **Validation:** Compare fitted FRF with measured FRF; check pole stability across orders; validate mode shapes (e.g. MAC with reference).

---

## When to use and limitations

**Use when** — **Known excitation** is available (impact hammer, shaker); accurate modal parameters (frequency, damping, mode shape) are required; laboratory modal testing or controlled field tests; input–output data are measured; frequency-domain analysis is preferred.

**Limitations** — **Requires known excitation:** not suitable for output-only (ambient) data; needs controlled input measurement. **Computational cost:** Rational fitting and root finding are non-trivial; full-band, multi-channel fitting can be resource-heavy. **Model order sensitivity:** Too low underfits, too high introduces spurious modes; stabilisation diagram helps but adds cost. **Numerical conditioning:** High-order polynomials can be ill-conditioned; orthogonal polynomials or partial fraction form help. **Frequency resolution:** Limited by measurement resolution; interpolation or sub-bin fitting may be needed.

---

## Engineering practice: practical notes

| Aspect | Notes |
|--------|--------|
| Excitation and measurement | Use impact hammer or shaker with known force; measure input and output simultaneously; multiple averages improve SNR and reduce noise bias. |
| FRF estimator choice | H1 for output noise; H2 for input noise; Hv for multiple inputs/outputs; coherence function indicates quality. |
| Frequency band and order | Fit band-by-band or globally; use stabilisation diagram to select order; true modes stabilize, spurious poles drift. |
| Numerical methods | Linear LS (Sanathanan–Koerner iteration) is fast but may need iteration; nonlinear optimization (Levenberg–Marquardt) is more accurate but slower; orthogonal polynomials improve conditioning. |
| Pole extraction | Root finding can be sensitive; use robust algorithms (e.g. companion matrix eigenvalue); check for stable poles (left half-plane for continuous, inside unit circle for discrete). |
| Mode shape scaling | Residues give relative mode shapes; absolute scaling requires mass matrix or normalization (e.g. unit modal mass, unit max component). |
| Validation | Compare fitted vs measured FRF (magnitude, phase); check coherence; validate poles with stabilisation diagram; compare mode shapes with reference (MAC). |

---

## Edge and online computing

**Suitability** — **Partially suited.** Fitting can be limited to a frequency band; reduced-order fit lowers cost; but needs known excitation and FRF estimation; rational fit and root finding are non-trivial; full-band multi-channel is resource-heavy.

**Potential** — **Band-limited fitting:** Fit only the frequency band of interest (e.g. first few modes) to reduce order and compute. **Reduced-order models:** Use low-order rational functions (few modes) for edge deployment; trade accuracy for speed. **Pre-computed models:** Fit offline, deploy poles/residues on edge for real-time FRF evaluation or modal filtering. **Tiered pipeline:** Estimate FRF and do coarse fit on edge; upload full data for high-order fit in cloud.

**Challenges** — **FRF estimation:** Needs input measurement and CPSD computation; H1/H2 estimators require buffer and FFT. **Rational fitting:** Nonlinear optimization is heavy; linear LS with iteration is lighter but still non-trivial. **Root finding:** Polynomial root finding (e.g. companion matrix eigenvalue) needs matrix operations; high order is costly. **Multi-channel:** Multiple FRFs increase data and compute; may need to fit one FRF at a time or use reduced MIMO. **Real-time:** Full FRF fitting is typically offline; edge can do coarse screening (e.g. "has frequency shifted?") with pre-fitted models, then upload for refinement.

**Practical strategy** — Use **pre-fitted models** on edge: fit FRF and rational function offline, deploy poles/residues; edge evaluates FRF or filters response in real time. For online updates, use **incremental fitting** (e.g. recursive LS on polynomial coefficients) or **band-limited updates** (fit only changed frequency bands). A **tiered approach**: edge does H1/H2 estimation and coarse low-order fit for screening; suspect segments uploaded for high-order fit and validation.
