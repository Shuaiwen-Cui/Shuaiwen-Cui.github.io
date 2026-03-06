# __STRUCTURAL HEALTH MONITORING__

## I. SHM and basic models

### 1.1 Brief introduction to SHM  
Structural Health Monitoring (SHM) aims at safety assessment, life prediction, and maintenance optimization of in-service structures. Sensors are deployed on structures to continuously collect vibration and response data, which are then combined with models and algorithms to identify and monitor the structural condition and potential damage.

### 1.2 Basic models in structural dynamics  
For an \(n\)-degree-of-freedom structure, the matrix form is

\[
\mathbf{M}\ddot{\mathbf{x}}(t) + \mathbf{C}\dot{\mathbf{x}}(t) + \mathbf{K}\mathbf{x}(t) = \mathbf{f}(t) \tag{1}
\]

**Notation:**

- \(\mathbf{x}(t),\dot{\mathbf{x}}(t),\ddot{\mathbf{x}}(t) \in \mathbb{R}^n\): displacement, velocity, acceleration;
- \(\mathbf{M},\mathbf{C},\mathbf{K} \in \mathbb{R}^{n \times n}\): mass, damping, stiffness matrices;
- \(\mathbf{f}(t) \in \mathbb{R}^n\): external force;
- \(n\): number of degrees of freedom.

Modal frequencies, damping ratios, and mode shapes are determined by these matrices and are key quantities for system and modal identification.

### 1.3 State-space models  
A linear time-invariant system can be written in state-space form. In continuous time,

\[
\dot{\mathbf{x}}(t) = \mathbf{A}_c \mathbf{x}(t) + \mathbf{B}_c \mathbf{u}(t) \tag{2}
\]

\[
\mathbf{y}(t) = \mathbf{C}_c \mathbf{x}(t) + \mathbf{D}_c \mathbf{u}(t) \tag{3}
\]

**Continuous-time (2)(3) notation:**

- \(\mathbf{x}(t) \in \mathbb{R}^n\): state vector;
- \(\mathbf{u}(t) \in \mathbb{R}^m\): input vector;
- \(\mathbf{y}(t) \in \mathbb{R}^p\): output vector;
- \(\mathbf{A}_c,\mathbf{B}_c,\mathbf{C}_c,\mathbf{D}_c\): continuous-time state, input, output, and feedthrough matrices;
- \(n,m,p\): dimensions of state, input, and output.

The discrete-time form \((k \in \mathbb{Z})\) is

\[
\mathbf{x}[k] = \mathbf{A}_d \mathbf{x}[k-1] + \mathbf{B}_d \mathbf{u}[k] \tag{4}
\]

\[
\mathbf{y}[k] = \mathbf{C}_d \mathbf{x}[k] + \mathbf{D}_d \mathbf{u}[k] \tag{5}
\]

**Discrete-time (4)(5) notation:**

- \(\mathbf{x}[k],\mathbf{u}[k],\mathbf{y}[k]\): state, input, and output at step \(k\);
- \(\mathbf{A}_d,\mathbf{B}_d,\mathbf{C}_d,\mathbf{D}_d\): discrete-time system matrices (from discretization);
- \(k\): discrete time index.

![State-space model block diagram](./SSM.png)  
*Figure 1: Block-diagram representation of a state-space model.*

In SHM, displacement and velocity are collected into a structural state \(\mathbf{z}(t) = [\mathbf{x}(t)^\top,\dot{\mathbf{x}}(t)^\top]^\top\); in discrete time the parameterized form is

\[
\mathbf{z}[k] = \mathbf{A}(\theta)\,\mathbf{z}[k-1] + \mathbf{B}(\theta)\,\mathbf{f}[k] \tag{6}
\]

\[
\mathbf{y}[k] = \mathbf{H}\mathbf{z}[k] + \mathbf{v}[k] \tag{7}
\]

**SHM parameterized form (6)(7) notation:**

- \(\mathbf{z}[k]\): structural state at step \(k\) (stack of displacement and velocity);
- \(\theta = \{\mathbf{M},\mathbf{C},\mathbf{K}\}\): structural parameters (mass, damping, stiffness matrices);
- \(\mathbf{A}(\theta),\mathbf{B}(\theta)\): state and input matrices depending on \(\theta\);
- \(\mathbf{f}[k]\): external force at step \(k\);
- \(\mathbf{H}\): observation matrix (sensor type and placement);
- \(\mathbf{v}[k]\): measurement noise.

This compact model is the basis for system and modal identification in SHM.

The following figures illustrate the block-diagram representation of this state-space model, and its specialization for SHM with parameter vector \(\theta\) and measurement operator \(\mathbf{H}\):

![State-space model for SHM](./SSM-SHM.png)  
*Figure 2: State-space model for SHM with parameters and measurement operator.*

---

## II. Common system/modal identification methods

Below are core, classical methods worth prioritising; each is given a one-sentence summary of the main idea and a one-sentence note on typical use.

### 2.1 PP

- **Name:** Peak Picking (PP).  
- **Main idea:** Read dominant peaks in the response power spectrum or Fourier spectrum for natural frequencies and use multi-channel amplitude ratios to estimate mode shapes.  
- **Typical use:** Quick, preliminary modal estimates when modes are well separated, damping is low, and SNR is good.  
- **Edge computing:** Well suited. Pros: FFT and peak search only; very low compute and memory; easy to run in real time. Cons: Low accuracy and robustness; sensitive to noise and close modes.

### 2.2 ARX / ARMAX (PEM)

- **Name:** ARX / ARMAX with prediction-error methods (PEM).  
- **Main idea:** Describe input–output (or output-only) dynamics with difference-equation models (ARX/ARMAX, etc.) and estimate parameters by minimising one-step-ahead prediction error; modal parameters follow from the identified model.  
- **Typical use:** Parametric system identification when input–output data are available; AR/ARMA variants for output-only cases.  
- **Edge computing:** Moderately suited. Pros: Compact model, controllable data size; modest cost for low order. Cons: Sensitive to order choice; MLE iteration can be heavy; numerical stability depends on implementation.

### 2.3 FRF

- **Name:** FRF curve fitting (Frequency Response Function curve fitting); common implementations include the Rational Fraction Polynomial (RFP) form.  
- **Main idea:** Fit measured FRFs in the frequency domain with rational functions (or orthogonal polynomials); extract modal frequency, damping, and mode shape from poles and residues.  
- **Typical use:** Laboratory modal testing with known excitation (impact, shaker) when accurate modal parameters are required.  
- **Edge computing:** Partially suited. Pros: Fitting can be limited to a frequency band; reduced-order fit lowers cost. Cons: Needs known excitation and FRF estimation; rational fit and root finding are non-trivial; full-band multi-channel is resource-heavy.

### 2.4 Transfer function fitting

- **Name:** Transfer function fitting (TF fitting).  
- **Main idea:** Fit a rational transfer function (e.g. \(H(s)=N(s)/D(s)\) or polynomial ratio in frequency) to measured FRF or input–output frequency-domain data; extract modal frequency, damping, and mode shape from the poles and residues of the fitted model.  
- **Typical use:** When FRF or frequency-domain data are available (e.g. after impact or shaker test); complements FRF estimation by providing a parametric modal model; used in classical modal analysis.  
- **Edge computing:** Partially suited. Pros: Can restrict to a frequency band; reduced-order fit lowers cost. Cons: Needs FRF or frequency data; rational fitting and root finding often require iteration or nonlinear optimisation; multi-channel full-band is resource-heavy.

### 2.5 FDD / EFDD

- **Name:** Frequency Domain Decomposition / Enhanced FDD.  
- **Main idea:** Decompose the output power spectral density matrix via SVD; dominant singular vectors approximate mode shapes and peak singular value curves give frequencies; EFDD refines damping by single-mode fitting in frequency bands.  
- **Typical use:** Output-only operational modal analysis (e.g. bridges, buildings under ambient excitation).  
- **Edge computing:** Moderately suited. Pros: Output-only; SVD per frequency line can be done in batches; manageable when channels and lines are limited. Cons: Buffer needed for PSD; EFDD damping fit adds an extra step.

### 2.6 Random Decrement (RDT / RDM)

- **Name:** Random Decrement Technique / Random Decrement Method (RDT / RDM).  
- **Main idea:** From long output-only response under ambient excitation, extract many short time segments triggered by a condition (e.g. level crossing) and average them to cancel the random part of excitation, yielding a **free-decay signature** (impulse-like response) that can be used to estimate modal frequency and damping (and, with multiple channels, mode shapes).  
- **Typical use:** Output-only operational modal analysis when you want a clean free-decay signal from ambient data; commonly used as a preprocessing step for **NExT/ERA** or for simple damping estimation.  
- **Edge computing:** Moderately suited. Pros: Mostly buffering + triggering + averaging; no heavy optimization. Cons: Needs sufficient data for averaging; trigger choice affects quality; multi-channel processing increases memory and bandwidth.

### 2.7 ERA / NExT-ERA

- **Name:** Eigensystem Realization Algorithm / Natural Excitation Technique – ERA.  
- **Main idea:** ERA forms Markov parameters and a Hankel matrix from impulse or free-decay response, then uses SVD and minimal realisation to obtain a state-space model and modal parameters; NExT-ERA under ambient excitation builds equivalent impulse responses from output auto- and cross-correlations, then applies ERA to them.  
- **Typical use:** ERA when impulse or free-decay data are available; NExT-ERA for output-only operational modal analysis (e.g. bridges, buildings under ambient excitation).  
- **Edge computing:** Moderately suited. Pros: Output-only; Hankel size can be limited; single SVD per run, no iteration. Cons: Correlation estimation needs buffer; order selection and multi-order trials increase cost.

### 2.8 SSI

- **Name:** Stochastic subspace identification (SSI: SSI-COV / SSI-DATA).  
- **Main idea:** Build block Hankel/Toeplitz matrices from output (or input–output) data, obtain observability/controllability subspaces via SVD, then recover a discrete state-space model and modal parameters from its system matrix.  
- **Typical use:** Output-only or input–output; robust modal identification with long, multi-channel data; widely used in operational modal analysis.  
- **Edge computing:** Poorly suited. Pros: Best identification quality; industry standard. Cons: Large Hankel and SVD; multi-order stabilization diagram; high compute and memory; heavy for low-power edge; may need reduced setup or cloud–edge split to offload.

### 2.9 BAYOMA

- **Name:** Bayesian Operational Modal Analysis (BAYOMA).  
- **Main idea:** In output-only setting, treat modal parameters as random variables; use Bayesian inference (e.g. MCMC, variational or Laplace approximation) to estimate their posterior from response data, yielding parameter estimates and uncertainty (e.g. credible intervals).  
- **Typical use:** Output-only OMA when uncertainty quantification or confidence intervals for modal parameters are needed; complements point-estimate methods such as SSI and FDD.  
- **Edge computing:** Poorly suited. Pros: Provides posterior and uncertainty for risk and decision-making. Cons: Bayesian computation (MCMC etc.) is heavy; high compute and memory; typically offline or cloud.


---

## III. Damage assessment methods

Methods below are ordered by dimension: detection first, then localisation, then quantification. Each entry is tagged with the dimension(s) it addresses (Detection / Localisation / Quantification).

### 3.1 Modal parameter change vs baseline

- **Main idea:** Compare current estimates of modal frequency, damping, or mode shape with a healthy baseline; form a scalar or vector index (e.g. frequency shift, damping change) and apply a threshold or simple statistical test to decide if damage is present.
- **Typical use:** When a healthy baseline exists and a quick yes/no or “significant change” decision is needed.
- **Dimensions:** Detection; can support rough quantification (e.g. magnitude of change).  
- **Edge computing:** Well suited. Pros: Comparison and threshold only; very low compute and memory. Cons: Depends on baseline quality; simple indices may miss damage. 
- **Online/offline:** Suited to online (real-time comparison).

### 3.2 Response- or feature-based anomaly detection

- **Main idea:** Use statistical features of multi-channel response (e.g. PCA residual, covariance change, Mahalanobis distance) or a novelty index; compare to the baseline distribution and flag anomaly when a threshold is exceeded, without explicit modal parameters.
- **Typical use:** Many channels, data-driven setup, when accurate modal identification is not required.
- **Dimensions:** Detection.  
- **Edge computing:** Moderately suited. Pros: Baseline/model can be trained offline; online step is projection and threshold only. Cons: Covariance or PCA cost grows with channel count. 
- **Online/offline:** Baseline offline, detection online; or fully offline.

### 3.3 MAC and mode-shape–based indices

- **Main idea:** Use Modal Assurance Criterion (MAC), mode shape difference, or mode slope to measure current vs baseline; regions or modes with large spatial change can indicate damage location.
- **Typical use:** When multiple mode shapes are available and damage detection plus rough localisation is needed.
- **Dimensions:** Detection, localisation.  
- **Edge computing:** Moderately suited. Pros: MAC etc. are inner products and norms; light compute. Cons: Needs current mode shapes; if online, depends on online modal ID. 
- **Online/offline:** Can be online (if modes are online) or offline.

### 3.4 Flexibility / stiffness matrix method

- **Main idea:** Flexibility and stiffness are inverses in the reduced DOF space; from identified modal parameters (frequency, mode shape, mass approximation) build the flexibility matrix or invert to get stiffness. Damage reduces stiffness and increases flexibility; compare current and baseline flexibility or stiffness (or their change) to find damage regions and to support quantification of stiffness loss.
- **Typical use:** When several modes and mass approximation are available; common for beams and frames; can be expressed as flexibility or stiffness, same in essence.
- **Dimensions:** Detection, localisation; can support quantification.  
- **Edge computing:** Partially suited. Pros: Clear physics; manageable when matrix size is small. Cons: Needs several modes and mass approximation; matrix build or inverse. 
- **Online/offline:** Often offline; can be online with light implementation.

### 3.5 Curvature and strain mode shape

- **Main idea:** Derive curvature (or strain mode) from displacement mode shape; damage causes local stiffness drop and a spike or discontinuity in curvature/strain, used to localise damage.
- **Typical use:** Dense sensors, beam- or plate-like structures, when localisation is the goal.
- **Dimensions:** Localisation.  
- **Edge computing:** Moderately suited. Pros: Curvature from mode difference or fit; modest compute; easy to distribute. Cons: Needs dense sensors and mode shape estimate. 
- **Online/offline:** Can be online (if modes online) or offline.

### 3.6 Model updating and stiffness inversion

- **Main idea:** Use a parameterised FE or reduced model (e.g. element stiffness reduction factors); fit or invert against current test data (modes, response) via optimisation or Bayesian inference to estimate stiffness loss and thus damage location and severity.
- **Typical use:** When a reliable structural model exists and quantitative location and severity are required.
- **Dimensions:** Detection, localisation, quantification.  
- **Edge computing:** Poorly suited. Pros: Quantitative. Cons: Optimisation/inversion is iterative and heavy; high compute and memory. 
- **Online/offline:** Suited to offline; online usually needs cloud or high-end node.

---

## IV. MODAL ANALYSIS TEACHING VIDEOS

Amirali Najafi Series

<iframe width="800" height="450" src="https://www.youtube-nocookie.com/embed/Og_VZ-PiBdM" frameborder="0" allowfullscreen></iframe>

NPTEL IIT Delhi Series

<iframe width="800" height="450" src="https://www.youtube-nocookie.com/embed/tzhSU9CTl08" frameborder="0" allowfullscreen></iframe>