# Transfer function fitting (TFF)

This page expands **2.4 Transfer function fitting** from the [SHM roadmap](../shm.en.md): fit a rational transfer function (e.g. \(H(s)=N(s)/D(s)\) or polynomial ratio in frequency) to measured FRF or input–output frequency-domain data, and extract modal frequency, damping, and mode shape from the poles and residues of the fitted model.

---

## Tutorial video

<iframe width="800" height="450" src="https://www.youtube-nocookie.com/embed/8-0M2FdRf3Q" frameborder="0" allowfullscreen></iframe>

---

## Concept

**Core idea** — The system **transfer function** \(H(s)\) (with \(s\) the complex frequency) or its frequency-domain form \(H(\omega)\) can be written as a **rational function** (ratio of polynomials) in modal analysis: its **poles** correspond to natural frequency and damping, and its **residues** to mode shapes. **Transfer function fitting** is the process of fitting such a rational function to FRF or frequency-domain input–output data from tests, so as to obtain a parametric model and extract all modal parameters. In short: **obtain FRF or frequency-domain data** → **choose rational form and fit** → **compute poles and residues** → **modal parameters**.

**Transfer function form** — In continuous time the rational form in \(s\) is common; when data are in the frequency domain, \(j\omega\) is used. The **Rational Fraction Polynomial (RFP)** form is:

\[H(\omega) = \frac{N(\omega)}{D(\omega)} = \frac{\sum_{k=0}^{m} a_k (j\omega)^k}{\sum_{k=0}^{n} b_k (j\omega)^k}, \tag{1}\]

**Notation for (1):**

- \(H(\omega)\): transfer function in frequency (i.e. FRF for SISO).
- \(a_k, b_k\): numerator and denominator polynomial coefficients (real or complex).
- \(m, n\): numerator and denominator orders; typically \(n \geq m\), with \(n\) on the order of twice the number of modes.
- \(j\): imaginary unit.

**Pole–residue form** (partial fractions):

\[H(\omega) = \sum_{r=1}^{N_m} \frac{R_r}{j\omega - s_r} + \text{(optional higher-order terms)}, \tag{2}\]

**Notation for (2):**

- \(N_m\): number of modes in the frequency band.
- \(s_r\): **pole** of mode \(r\) (complex), \(s_r = -\zeta_r \omega_r \pm j\omega_r\sqrt{1-\zeta_r^2}\).
- \(R_r\): **residue** of mode \(r\) (complex, related to mode shape).
- \(\omega_r, \zeta_r\): undamped natural frequency and damping ratio.

**Relation to FRF curve fitting** — FRF curve fitting is the rational fitting of **measured FRF** in the frequency domain and is one instance of transfer function fitting; this page emphasises the **transfer function** view (rational \(H(s)/H(\omega)\), poles and residues) and covers the general process of obtaining a parametric modal model from FRF or frequency-domain input–output data.

---

## Algorithm in brief

**Data source** — The fit is usually performed on:

- **Measured FRF:** from input–output via H1/H2/Hv estimators (see [FRF](../FRF/frf.en.md)); or
- **Frequency-domain input–output:** Fourier coefficients at each frequency, fitting \(Y(\omega)/U(\omega)\) or equivalent.

**Rational function fitting** — Fit \(\hat{H}(\omega_k)\) at discrete frequencies \(\omega_k\). The problem is **nonlinear** in the coefficients (or poles/residues). Common approaches:

1. **Linear LS on numerator/denominator:** Write (1) as \(D(\omega)H(\omega)=N(\omega)\), linearise in \(a_k,b_k\), solve LS; iterate with weighting if needed (e.g. Sanathanan–Koerner).
2. **Nonlinear optimisation:** Optimise poles \(s_r\) and residues \(R_r\) (or polynomial coefficients) directly, minimising \(\sum_k |H(\omega_k)-\hat{H}(\omega_k)|^2\) with Gauss–Newton, Levenberg–Marquardt, etc.
3. **Orthogonal polynomials:** Use an orthogonal basis (e.g. Forsythe) instead of monomials to improve conditioning, then find roots for poles.

**Poles and residues** — Poles \(s_r\) are the roots of the fitted denominator \(D\), or read from form (2). Residues \(R_r\) come from partial fraction expansion or the fitted model; for MIMO, residues form a matrix whose columns (or rows) give mode shapes.

**Modal parameters** — \(\omega_r = |s_r|\) (or solve for \(\omega_r,\zeta_r\) from real and imaginary parts of \(s_r\)); mode shapes from the residue matrix and output mapping.

---

## Procedure (outline)

1. **Data preparation:** Obtain FRF or frequency-domain input–output data (known excitation required; see FRF estimation).
2. **Band and order:** Choose frequency band of interest; choose numerator/denominator orders \(m,n\) or number of modes \(N_m\), e.g. via stabilisation diagram (poles vs order) or AIC/BIC.
3. **Rational fit:** Fit (1) or (2) over the chosen band using linear LS (with iteration if needed), nonlinear optimisation, or orthogonal polynomials.
4. **Pole extraction:** Get poles \(s_r\) from denominator roots or from (2); drop unstable or clearly spurious poles.
5. **Residues and mode shapes:** Compute residues \(R_r\) at each pole; for multi-channel, form mode shapes from the residue matrix and normalise.
6. **Validation:** Compare fitted curve with measured FRF; check pole stability with order; compare mode shapes with reference (e.g. MAC) if available.

---

## When to use and limitations

**Use when** — **FRF or frequency-domain input–output data** are available (e.g. impact or shaker test); a **parametric modal model** (poles, residues, or rational coefficients) is needed; classical modal analysis, model updating, or control design; frequency-domain, rational representation is preferred.

**Limitations** — **Requires known excitation:** not for output-only (ambient) data; controlled input and FRF (or frequency-domain I/O) estimation are needed. **Order sensitivity:** too low underfits, too high introduces spurious poles; stabilisation diagram or criteria help. **Compute and numerics:** rational fitting and root finding are non-trivial; high order can be ill-conditioned; orthogonal basis or partial fraction form help. **Multi-channel:** multiple FRFs or MIMO increase data and compute; fit channel-by-channel or use reduced MIMO.

---

## Engineering practice: practical notes

| Aspect | Notes |
|--------|--------|
| Data and FRF | Same as [FRF](../FRF/frf.en.md): known excitation, H1/H2/Hv estimation, multiple averages for SNR; coherence for quality. |
| Rational form | RFP (1) is simple to implement; pole–residue (2) has clear physical meaning and can be numerically better; orthogonal polynomials improve conditioning. |
| Order and stabilisation | Try several orders, plot poles vs order; stable poles are physical modes, drifting or scattered ones are often spurious. |
| Fitting method | Linear LS + iteration is fast and easy; nonlinear optimisation is more accurate but needs initial guess; use linear solution as initial guess if needed. |
| Poles and residues | Use robust root finding (e.g. companion matrix eigenvalue); residues give relative mode shapes, absolute scaling needs mass or normalization. |
| Validation | Compare fitted and measured FRF (magnitude, phase); check poles with stabilisation diagram; use MAC etc. for mode shapes in multi-channel. |

---

## Edge and online computing

**Suitability** — **Partially suited.** Band and order can be limited to reduce compute; but FRF or frequency-domain data are still required (hence input measurement and estimation), and rational fitting plus root finding remain non-trivial on resource-limited devices.

**Potential** — **Band and reduced order:** Fit only the band and number of modes of interest to reduce order and compute. **Pre-fitted deployment:** Fit offline, deploy poles/residues or rational coefficients on edge for real-time FRF evaluation or modal filtering. **Tiered:** Edge does FRF estimation and coarse low-order fit for screening; suspect data uploaded for high-order fit and validation.

**Challenges** — **FRF and data:** Input measurement, buffer, CPSD/FFT. **Fitting and roots:** Nonlinear or iterative fitting and polynomial root finding are not light on MCUs. **Real-time:** Full fitting is usually offline; edge is better with pre-fitted models for coarse checks (e.g. frequency shift?), then upload for refinement.

**Practical strategy** — Use **pre-fitted models** on edge: do transfer function fitting offline, deploy poles/residues; edge only evaluates FRF or response. For online updates, use **recursive or incremental fitting** (e.g. recursive LS on polynomial coefficients) or **band-limited updates**; or a **tiered** setup with coarse fit on edge and refined fit in the cloud.

---

## Relation to FRF curve fitting

- **Transfer function fitting (this page):** Emphasises fitting a rational **transfer function** \(H(s)/H(\omega)\) to data (FRF or frequency-domain I/O) to obtain poles, residues, and modal parameters; applies to classical modal analysis, modelling, and control.
- **FRF curve fitting:** Focuses on rational curve fitting of **measured FRF**; the procedure and algorithms are the same as here, and can be seen as transfer function fitting when the data source is FRF; see [FRF](../FRF/frf.en.md).

The two are the same in mathematics and implementation: rational fitting plus pole/residue extraction; this page stresses the “transfer function” and general data source, the FRF page the “FRF estimation + curve fitting” workflow.
