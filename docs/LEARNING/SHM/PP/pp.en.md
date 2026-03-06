# Peak Picking (PP)

This page expands **2.1 PP** from the [SHM roadmap](../shm.en.md): how to quickly estimate modal frequencies and mode shapes by picking peaks on a frequency-domain function (response spectrum or FRF).

---

## Tutorial video

<iframe width="800" height="450" src="https://www.youtube-nocookie.com/embed/5wpFgU_oAss" frameborder="0" allowfullscreen></iframe>

---

## Concept

**Core idea** — On a frequency-domain representation of the response (or of the FRF), the **frequencies of the dominant peaks** are the **damped natural frequencies** (for light damping, ≈ undamped \(\omega_n\)). At each peak frequency \(\omega_r\), the **ratio of amplitudes across channels** approximates the **relative mode shape** \(\phi_j/\phi_{\mathrm{ref}}\). So: peaks → frequencies; amplitude ratio at each peak → mode shape.

**Relation to FRF** — Peak picking is the same idea in two settings. In **FRF-based modal analysis** (e.g. impact or shaker tests), \(|H(\omega)|\) has peaks at the natural frequencies, and the FRF values across DOFs at each peak give the mode shape. In **output-only** (ambient or unknown excitation), the response FFT or power spectrum also peaks at those frequencies, and the amplitude ratio across channels approximates the mode shape. So **PP = “pick peaks on a frequency-domain function”**—either the FRF \(H(\omega)\) (input–output) or the response spectrum (output-only). With FRF you have a direct transfer function; with output-only data the response spectrum inherits the system’s poles (resonances).

**Formula** — Under single-mode dominance at \(\omega_r\), the mode-shape ratio equals the amplitude ratio:

\[\frac{\phi_j}{\phi_{\mathrm{ref}}} \approx \frac{A_j(\omega_r)}{A_{\mathrm{ref}}(\omega_r)},\]

where \(A_j\) is the amplitude at sensor \(j\) at \(\omega_r\). The quantity used must be **amplitude** (same dimension as the vibration at that frequency), not power—so \(A_j\) is given by \(|X_j|\) or \(\sqrt{S_j}\) as below, not by \(S_j\) directly.

---

## Implementation

**FFT vs power spectrum (PSD)** — Use the same physical quantity for all channels.

- **FFT (Fourier spectrum):** \(|X_j(\omega_r)|\) is proportional to amplitude; the ratio across channels is the amplitude ratio. So

\[\frac{\phi_j}{\phi_{\mathrm{ref}}} \approx \frac{|X_j(\omega_r)|}{|X_{\mathrm{ref}}(\omega_r)|}.\]

- **PSD:** \(S\) has dimension [signal²/frequency], and \(S \propto |X|^2\), so amplitude \(\propto \sqrt{S}\). Use the **ratio of square roots**:

\[\frac{\phi_j}{\phi_{\mathrm{ref}}} \approx \frac{\sqrt{S_j(\omega_r)}}{\sqrt{S_{\mathrm{ref}}(\omega_r)}}.\]

  Using \(S_j/S_{\mathrm{ref}}\) gives \((\phi_j/\phi_{\mathrm{ref}})^2\)—wrong dimension and meaning.

**Phase and sign (output-only)** — Without known input, **absolute phase** is not available. If only **magnitude** is used (PSD or \(|X_j|\)), the mode-shape ratio is estimated only up to **absolute value** \(|\phi_j/\phi_{\mathrm{ref}}|\): **direction (sign)** is unknown (in phase vs out of phase with the reference). So amplitude-only PP yields a mode shape with arbitrary global sign. If **complex spectra** are available, the **relative phase** from the complex ratio \(X_j/X_{\mathrm{ref}}\) is preserved (the unknown input phase cancels), so sign can be recovered; the procedure’s “use the complex ratio to retain relative phase” is for this.

---

## Procedure

1. Compute FFT or power spectrum (e.g. Welch) for each channel.
2. Identify dominant peaks in the frequency range of interest; record their frequencies \(\omega_r\).
3. Choose a reference channel. At each \(\omega_r\), read **amplitude** at every channel: \(|X_j(\omega_r)|\) when using FFT, or \(\sqrt{S_j(\omega_r)}\) when using PSD (do not use \(S_j\) directly).
4. Form the ratio of each channel’s amplitude to the reference; normalize the vector (e.g. by max or unit norm).
5. (Optional) If complex spectra are available, use the complex ratio to retain relative phase (and thus sign).

---

## When to use and limitations

- **Use when:** Modes are well separated, damping is low, and SNR is good; suited to quick checks or first estimates.
- **Limitations:** No damping estimate; peak overlap or close modes are hard to separate; sensitive to noise and reference choice; with magnitude only, mode shape is only \(|\phi_j/\phi_{\mathrm{ref}}|\) (sign unknown).

---

## Practical shortcomings and mitigations

| Shortcoming | Mitigation | Notes |
|-------------|------------|--------|
| Spectrum variance, noise | Segment averaging (e.g. Welch), band-pass filter | Averaging stabilises peaks; band-pass suppresses out-of-band noise. |
| Peak reliability | Peak prominence, local fit (e.g. parabola) | Prominence rejects broad bumps; local fit gives sub-bin frequency. |
| Multiple peaks | Multi-window / multi-segment voting, prior band | Persistent peaks are more credible; search only in known band when available. |
| Reference and sign | Reference away from nodes, use complex ratio | Choose sensor with large amplitude; complex spectrum retains phase/sign. |
| Close modes | Use PP for screening only, upload suspect data | “Significant change yes/no”; use FDD/SSI for dense bands offline or in cloud. |
| Window length | Resolution–latency trade-off | \(\Delta f\) smaller than minimum mode spacing; shorten window when real-time matters. |

**Summary** — PP’s strength is simplicity and low compute, suited to edge and online screening; its weaknesses come from “no model, just peaks”. Averaging, filtering, robust peak picking, and reference choice improve usability; close modes and damping are better handled downstream.

---

## Edge and online computing

**Why it fits** — Only FFT (\(O(N\log N)\)) and peak search (\(O(N)\)), no matrix factorisations; memory linear in window length; suitable for MCUs and low-power SoCs. Sliding or block FFT gives predictable latency (“compute as you sample, alert on device”). When modes are well separated and excitation stable, PP at the edge can screen “has frequency or mode shape shifted?” and upload suspect segments for FDD/SSI, as the first stage of a tiered SHM pipeline.

**Core peak-picking algorithm and edge implementation** — Peak picking is local-maximum search on the magnitude (or power) spectrum: for each frequency bin, compare with neighbours (e.g. one or two bins on each side); a bin larger than both sides is a candidate peak. Use **peak prominence** or a simple threshold to reject spurious peaks, harmonics, and noise; for **sub-bin frequency** accuracy, fit a parabola in the neighbourhood of the peak and take the vertex. On the edge: use **segmented FFT** (e.g. Welch with 50% overlap) or **sliding-window FFT**; for each segment or window compute the FFT and run the above peak search on the magnitude spectrum—no need to keep phase for frequency and amplitude-ratio estimates. Peak detection is only comparisons and a few add/multiply operations, so it is amenable to fixed-point or integer math; with a fixed window length, latency and memory are deterministic, so an MCU can "acquire one window, process one window, report one window" and form a tiered "edge screening + upload of suspect segments" pipeline with the cloud or host.

**Challenges** — Edge SNR is often lower; trade off window length (resolution and noise) vs latency. A single peak can mix two close modes; PP is better for “significant change yes/no” than fine decomposition. No damping from PP; if damage shows mainly in damping, combine with other features or uplink. Choose a robust reference (away from nodes) from prior or offline analysis.
