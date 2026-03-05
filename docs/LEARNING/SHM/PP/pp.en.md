# Peak Picking (PP)

This page expands **2.1 PP** from the [SHM roadmap](../shm.en.md): how to quickly estimate modal frequencies and mode shapes using peak picking.

---

## Basic idea

- **Frequency:** Compute the FFT or power spectrum of each channel; the frequencies of the dominant peaks correspond to the **damped natural frequencies** of the modes (for light damping, approximately the undamped natural frequency \(\omega_n\)), and are used as the modal frequency estimates.

- **Mode shape:** At each peak frequency \(\omega_r\), form the ratio of each channel’s **amplitude** to the reference; that ratio approximates the relative mode shape. Physically the quantity used must be “amplitude” (same as the vibration amplitude at that frequency component), not raw power—**the ratio of FFT magnitudes and the ratio of PSD values correspond to different physical quantities, so mode-shape estimation must treat them differently** (below).

Under a single-mode-dominant approximation, the response at that frequency is dominated by that mode, so the mode-shape ratio equals the amplitude ratio:

\[
\frac{\phi_j}{\phi_{\mathrm{ref}}} \approx \frac{A_j(\omega_r)}{A_{\mathrm{ref}}(\omega_r)},
\]

where \(A_j\) is the amplitude at sensor \(j\) at \(\omega_r\) (given by \(|X_j|\) or \(\sqrt{S_j}\) as below).

**FFT vs power spectrum (PSD)**

- **FFT (Fourier spectrum)**  
  The ratio of \(\lvert X_j(\omega_r) \rvert\) across channels equals the **ratio of amplitudes** at that frequency (the modulus of \(X\) is proportional to amplitude, with the same scaling for all channels for a given window). So take \(A_j = \lvert X_j(\omega_r) \rvert\) to estimate the mode-shape ratio:

$$
\frac{\phi_j}{\phi_{\mathrm{ref}}} \approx \frac{\lvert X_j(\omega_r) \rvert}{\lvert X_{\mathrm{ref}}(\omega_r) \rvert}.
$$

- **Power spectrum (PSD)**  
  \(S\) has dimension [signal²/frequency] (e.g. (m/s²)²/Hz); \(S \propto \lvert X \rvert^2\), so **amplitude** \(\propto \sqrt{S}\)—only the ratio of \(\sqrt{S_j}\) across channels equals the amplitude ratio. With \(A_j = \sqrt{S_j(\omega_r)}\), use the **ratio of square roots**:

$$
\frac{\phi_j}{\phi_{\mathrm{ref}}} \approx \frac{\sqrt{S_j(\omega_r)}}{\sqrt{S_{\mathrm{ref}}(\omega_r)}}.
$$

Using \(S_j/S_{\mathrm{ref}}\) instead gives \((\phi_j/\phi_{\mathrm{ref}})^2\)—wrong dimension and wrong physical meaning.

---

## Procedure (outline)

1. FFT or power spectrum (e.g. Welch) for each channel.
2. Identify dominant peaks in the frequency range of interest; record their frequencies \(\omega_r\).
3. Choose a reference channel; at each \(\omega_r\), read the **amplitude** at every channel: use \(|X_j(\omega_r)|\) when using FFT, or \(\sqrt{S_j(\omega_r)}\) when using PSD (do not use \(S_j\) directly—wrong dimension).
4. Form the ratio of each channel’s amplitude to the reference; normalize the vector (e.g. by max or unit norm).
5. (Optional) If complex spectra are available, use the complex ratio to retain relative phase.

---

## When it works and limitations

- **Use when:** Modes are well separated, damping is low, and SNR is good; suited to quick checks or first estimates.
- **Limitations:** No damping estimate; peak overlap or close modes are hard to separate; sensitive to noise and reference choice; mode shape is an amplitude-ratio approximation only.

---

## Engineering practice: shortcomings and mitigations

In practice, PP’s weaknesses and the measures that can mitigate them are as follows; understanding them helps choose when PP is “good enough” versus when to add complexity.

**Main shortcomings** — Spectrum variance and noise (single-segment FFT has high variance, risk of false/missed detections); frequency resolution vs window length (\(\Delta f \approx 1/T\): too short blurs close modes, too long increases latency); peak-picking strategy (simple threshold mistakes harmonics/noise for peaks, discrete bins add scalloping); reference and mode-shape sign (amplitude ratio gives only \(|\phi_j/\phi_{\mathrm{ref}}|\), and reference near a node makes the ratio unstable); close modes (PP cannot disentangle when two frequencies are within \(\sim 1/T\) or damping bandwidth); no damping (PP does not yield damping for assessment).

**Practical mitigations**

| Shortcoming | Mitigation | Notes |
|-------------|------------|--------|
| Spectrum variance, noise | Segment averaging (e.g. Welch), band-pass filter | Averaging stabilises peaks; band-pass limits band and suppresses out-of-band noise. |
| Peak reliability | Peak prominence, local fit (e.g. parabola) | Prominence rejects broad bumps; local fit gives sub-bin frequency. |
| Multiple peaks | Multi-window / multi-segment voting, prior band | Persistent peaks are more credible; search only in known band when available. |
| Reference and sign | Reference away from nodes, use complex ratio | Choose sensor with large amplitude; complex spectrum retains phase/sign. |
| Close modes | Use PP for screening only, upload suspect data | “Significant change yes/no”; use FDD/SSI for dense bands offline or in cloud. |
| Window length | Resolution–latency trade-off | \(\Delta f\) smaller than minimum mode spacing; shorten window when real-time matters. |

**Summary** — PP’s strength is simplicity and low compute, suited to edge and online screening; its weaknesses come from “no model, just peaks”. Averaging, filtering, robust peak picking, and reference choice improve usability; close modes and damping are better handled downstream.

---

## Edge and online

**Why it fits** — Compute and memory: only FFT (\(O(N\log N)\)) and peak search (\(O(N)\)), no matrix factorisations or iterations, memory linear in window length; suitable for MCUs and low-power SoCs. Streaming and latency: sliding or block FFT with fixed window gives predictable latency, “compute as you sample, alert on device”. Potential: when modes are well separated and excitation stable, PP at the edge can screen “has frequency or mode shape shifted noticeably?” and upload suspect segments for FDD/SSI, reducing bandwidth and cloud load as the first stage of a tiered SHM pipeline.

**Challenges and trade-offs** — Noise and robustness: edge SNR is often lower; trade off window length (resolution and noise rejection) vs latency; light improvements (smoothing, multi-window voting) help but add compute. Close modes and mode-shape quality: a single peak can mix two modes; PP is better for “significant change yes/no” than fine decomposition. No damping: if damage shows mainly in damping, combine with other features or uplink. Reference: choose a robust reference (away from nodes) from prior or offline analysis, or accept limited reliability for some modes.
