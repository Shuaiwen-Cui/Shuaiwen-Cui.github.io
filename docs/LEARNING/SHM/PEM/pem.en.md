# ARX / ARMAX and prediction-error methods (PEM)

This page expands **2.2 ARX / ARMAX (PEM)** from the [SHM roadmap](../shm.en.md): describe input–output or output-only dynamics with difference-equation models (ARX, ARMAX, etc.), estimate parameters by minimising one-step-ahead prediction error, and obtain modal parameters from the identified model.

---

## Concept

**Core idea** — Use discrete-time difference equations (ARX, ARMAX, AR, ARMA) to describe the system; treat the **one-step-ahead prediction error** as a function of the model parameters and estimate parameters by **minimising a prediction-error criterion** (e.g. least squares, maximum likelihood); then obtain the **poles** of the identified discrete transfer function or state-space form, and hence **modal frequency and damping ratio** (and **mode shapes** when using MIMO or state-space). So PEM is: **model + prediction-error criterion + optimisation** → parameter estimate → modal parameters.

**Model forms (brief)** — Let \(q^{-1}\) be the delay operator (\(q^{-1} y[k] = y[k-1]\)).

- **ARX:** \(A(q^{-1}) y[k] = B(q^{-1}) u[k] + e[k]\). Output \(y\) is driven by past input \(u\), past output, and equation error \(e\); parameters enter **linearly** in \(A,B\), so **least squares (LS)** gives the solution in one step, no iteration.
- **ARMAX:** Adds a moving-average term in the noise, \(C(q^{-1}) e[k]\), i.e. \(A(q^{-1}) y[k] = B(q^{-1}) u[k] + C(q^{-1}) e[k]\). Parameters in \(A,B,C\) are **nonlinear** (because of \(C\)), so **iterative optimisation** (e.g. Gauss–Newton, prediction-error minimisation, or MLE) is required.
- **Output-only (AR / ARMA):** No input term: \(A(q^{-1}) y[k] = e[k]\) (AR) or \(A(q^{-1}) y[k] = C(q^{-1}) e[k]\) (ARMA). Used when only response is measured (e.g. ambient excitation); estimate \(A\) (and \(C\)) from response autocorrelation or directly by PEM, then get modal frequency and damping from the poles.

**Relation to state space** — Difference-equation and discrete state-space models are equivalent (convertible); PEM can also be applied directly to state-space parameters. Here the focus is on input–output forms (ARX/ARMAX) for simplicity and the fact that ARX has a closed-form LS solution.

---

## Algorithm in brief

**Prediction and prediction error** — For given data and parameters, the model yields a **one-step-ahead prediction** \(\hat{y}[k|k-1]\); the **prediction error** is \(\varepsilon[k] = y[k] - \hat{y}[k|k-1]\). PEM chooses parameters to minimise a criterion (usually sum of squared errors or negative log-likelihood).

**ARX: least squares** — ARX can be written as a linear regression in the parameter vector; the prediction error is **linear** in the parameters. Minimising \(\sum_k \varepsilon[k]^2\) leads to the **normal equations**; solve the linear system (or use **recursive least squares, RLS**, for online/streaming). No iteration, no initial guess; numerically stable and cheap.

**ARMAX / ARMA: iterative optimisation** — The prediction \(\hat{y}[k|k-1]\) depends on past \(\varepsilon\), which in turn depends on parameters, so the problem is **coupled** and **nonlinear** in parameters. Use **iterative methods** (e.g. Gauss–Newton, Levenberg–Marquardt, or gradient-based MLE) until convergence. Initialise from ARX or a high-order ARX; ensure convergence and numerical stability (e.g. keep \(A,C\) stable).

**Modal parameter extraction** — From the estimated \(A(q^{-1})\) (and \(B,C\)), find the **roots** of \(A(z^{-1})=0\) (discrete poles). Map poles to continuous time (e.g. \(s = \ln(z)/T_s\)) to get **continuous poles** \(s_r = -\zeta_r \omega_r \pm j\omega_r\sqrt{1-\zeta_r^2}\), hence **undamped natural frequency** \(\omega_r\) and **damping ratio** \(\zeta_r\). For multiple outputs, **mode shapes** come from the MIMO model or the equivalent state-space realisation (observability matrix or residue directions).

---

## Procedure (outline)

1. **Data and structure:** Prepare input–output sequences \(\{u[k],y[k]\}\) (or output-only \(\{y[k]\}\)); choose model type (ARX / ARMAX / AR / ARMA) and **orders** (\(n_a,n_b,n_c\), etc.).
2. **Order selection:** Use AIC, BIC, cross-validation, or identify at several orders and plot a **stabilisation diagram** (poles vs order); pick an order where true modes stabilise and spurious poles are few.
3. **Estimation:** ARX: solve normal equations (LS or RLS); ARMAX/ARMA: minimise prediction error or MLE iteratively, initialised e.g. by ARX or zero.
4. **Validation:** Residual whiteness, fit quality, or comparison with held-out data; impose stability on \(A,C\) in the optimisation if needed.
5. **Modal extraction:** Roots of \(A(z^{-1})\) → discrete poles → map to continuous poles → \(\omega_r,\zeta_r\); for multiple channels, get mode shapes from MIMO or state-space form.

---

## When to use and limitations

**Use when** — Input–output data are available and a **parametric model** is needed for control or system ID; modal **frequency and damping** (and mode shapes for MIMO) are required; data length and SNR are adequate; a difference-equation form and a chosen order are acceptable. For output-only data, AR/ARMA can be used for OMA under ambient excitation.

**Limitations** — **Order sensitivity:** too low underfits, too high introduces spurious modes; use stabilisation diagram or criteria. **Excitation:** input–output identification needs **persistent excitation (PE)**; otherwise parameters are not identifiable or ill-conditioned. **ARMAX/ARMA:** iteration, initialisation, convergence, and numerical stability (stable \(A,C\)) require care. **Mode shapes:** SISO gives only frequency and damping; need multiple channels or MIMO for mode shapes. **Nonlinearity and time variation:** model is LTI; strong nonlinearity or fast time variation need other methods.

---

## Engineering practice: practical notes

| Aspect | Notes |
|--------|--------|
| Excitation and identifiability | Input must be persistently exciting (PE), covering the band of interest; output-only assumes excitation approximately white or broadband. |
| Order and stabilisation diagram | Try several orders; use AIC/BIC or stabilisation diagram; stable poles that persist across orders are true modes; drifting or scattered poles are often spurious. |
| Initialisation and iteration (ARMAX/ARMA) | Initialise from high-order ARX or zero; monitor convergence and residuals; regularise or constrain stability if needed. |
| Numerical stability | Normal equations can be ill-conditioned; use QR/SVD for LS; in iteration, keep roots of \(A,C\) inside the unit circle. |
| Sampling and delay | Sampling rate must satisfy the band of interest (e.g. Nyquist); when using \(B(q^{-1})\), set input–output delay \(n_k\) appropriately. |
| Multi-channel and mode shapes | Multi-output ARX/ARMAX or equivalent state-space yield mode shapes; need enough channels and excitation to identify each mode. |

---

## Edge and online

**Suitability** — **Moderately well suited.** The model is compact (difference equation) and data size is controllable; **low-order ARX** needs only LS, no iteration, so compute and memory are modest and suitable for MCUs. ARMAX/ARMA and high-order ARX require iteration and larger matrices, so edge deployment needs limits (max order, max iterations) or ARX-only approximation.

**Potential** — **Recursive least squares (RLS)** can update ARX parameters online as new samples arrive, suited to streaming and edge. Low-order ARX can be implemented in fixed-point or integer arithmetic; identified \(\omega_r,\zeta_r\) can be compared with a baseline for online “frequency/damping shift?” screening, with suspect segments uploaded for heavier methods (e.g. SSI). ARX can also provide initial values for ARMAX or state-space PEM (coarse on edge, refine in cloud).

**Challenges** — **Order selection** on the edge is hard to do with full stabilisation diagrams; often rely on prior or fixed low order. **MLE/iteration** (ARMAX, ARMA) on MCUs needs tight control of iteration count and numerical stability. **Persistent excitation** under ambient excitation may be weak, reducing quality. **Multi-channel and mode shapes** increase data and compute. A practical strategy is a **tiered pipeline**: low-order ARX at the edge for fast screening, with suspect data uploaded for higher-order or ARMAX/SSI.
