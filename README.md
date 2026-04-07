# CST — Compressible Spacetime Dynamics

**A unified framework for gravitational coupling from planetary to galactic scales**

*Michele Vizzutti — Independent Researcher, Udine, Italy*

---

## Overview

Compressible Spacetime Dynamics (CST) is a theoretical framework that treats spacetime as a compressible barotropic fluid, where matter induces local density enhancements that modify the effective gravitational coupling. Unlike General Relativity, where G is a universal constant, CST predicts a mass- and epoch-dependent coupling:

<p align="center">
  <img src="https://latex.codecogs.com/svg.image?G_{\rm&space;eff}(M,z)=G_N\left\{1+[1-w(M)]\,\alpha\left(\frac{M}{M_\odot}\right)^\beta&space;f(z)\right\}" alt="G_eff formula"/>
</p>

where:
- **w(M) = exp(− |M / M☉ − 1|)** suppresses deviations at the solar mass scale
- **α = 0.279 ± 0.012** is the coupling intensity
- **β = 0.685 ± 0.018 ≈ 2/3** is the mass scaling exponent
- **f(z) = [H(z) / H₀] / [1 + (z/30)³]** encodes cosmological evolution

The theory has been validated on **24,496 astronomical systems** across three scales, plus **129 galaxies** with **2,931 rotation curve data points**.

---

## Key Results

### Compact Systems (Exoplanets + Binary Stars)

| Dataset | Objects | R² | Cross-validation |
|---------|---------|-----|------------------|
| NASA Exoplanets | 4,585 | 96.04% | 94.95% (K-fold) |
| Gaia DR3 Binaries | 16,980 | 96.96% | — |
| Synthetic Validation | — | 99.19% | — |
| **Combined Multi-Scale** | **21,565** | **97.73%** | **No overfitting** |

### Galactic Scale (SPARC Database)

CST treats the supermassive black hole (SMBH) as a **vortex** in the spacetime fluid (Lamb-Oseen/Burgers profile), producing a three-component velocity model:

> **v²(r) = v²_bar(r) · (1 + ε_PDMF) + V²_flat · LO(r/r_c)² · Ψ_gal(r)**

**Model G** — with **5 global parameters** and **zero per-galaxy free parameters**:

| Model | Parameters | RMSE (km/s) | Improvement |
|-------|-----------|-------------|-------------|
| Baryons only | 0 | 61.9 | 1.0× |
| MOND (a₀) | 1 | ~25 | ~2.5× |
| **CST Model G** | **5 global** | **20.0** | **3.1×** |
| CST individual fit | 258 | 18.4 | 3.4× |

Cross-validation confirms **no overfitting**: train RMSE = 20.0 ± 0.6, test RMSE = 20.1 ± 3.1, gap = +0.14 km/s.

### CST vs MOND

| Test | CST | MOND | Winner |
|------|-----|------|--------|
| RAR scatter | 0.121 dex | 0.163 dex | **CST** |
| BTFR slope | 3.53 | 4.0 (predicted) | Comparable |
| Origin of a₀ | Emergent (vortex) | Fundamental (unexplained) | **CST** |
| Wide binaries | No anomaly predicted | Anomaly via EFE | Testable with Gaia DR4 |
| Physical mechanism | SMBH vortex | Phenomenological | **CST** |

---

## The Universal Exponent β = 2/3

The most fundamental result: the exponent β = 2/3 is **not a fitted parameter** — it is a geometric invariant of three-dimensional space:

> **β = (D−1) / D = 2/3 for D = 3 spatial dimensions**

It appears independently in four contexts:

| Context | Observed β | Predicted | Deviation |
|---------|-----------|-----------|-----------|
| Stellar compression (virial theorem) | 0.685 ± 0.018 | 2/3 | 2.8% (1.0σ) |
| Vortex core radius (r_c ∝ M_BH^β) | 0.646 | 2/3 | 3.1% |
| Galactic interference (Ψ_gal) | 2/3 (fixed) | 2/3 | +0.28 km/s penalty |
| Cantor set / cosmic fractal (d_f/D) | ~0.63–0.67 | 2/3 | — |

**Physical meaning:** 2/3 is the ratio of surface (boundary) to volume (bulk) degrees of freedom. The sphere minimises surface area for a given volume (isoperimetric inequality). Self-gravitating systems exchange energy through their surface and store it in their volume — the ratio is always (D−1)/D.

This connects to fractal geometry: the Cantor set retains 2 out of 3 elements at each iteration; the cosmic matter distribution has fractal dimension d_f ≈ 2 on scales up to ~100 Mpc, giving d_f/D = 2/3.

---

## The Galactic Interference Function

The stellar population between the SMBH and any radius r perturbs the vortex, using the same mathematical structure as binary star interference:

> **Ψ_gal(r) = 1 + γ · [4q/(1+q)²]^p · exp(−r / r₀·R_disk) · (M_BH / 10⁷ M☉)^(2/3)**

where q(r) = M_bar(<r) / M_BH is the local mass ratio. The factor **4q/(1+q)²** peaks at q = 1 (equal mass) and vanishes for q → 0 (dwarf galaxies) — the same binary-like interference term from Section 2.6 of the manuscript.

**Optimised global parameters:** C = 18.5, γ = 267, p = 0.44, r₀ = 2.6, η_bul = 0.09

---

## Predictions

### Confirmed
- ✅ Galaxies without SMBH lack dark matter (NGC 1052-DF2, DF4, M33)
- ✅ JWST massive galaxies at z = 10–15 (accelerated structure formation with G_eff ≈ 1.3 G_N)

### Testable with Current/Near-Future Instruments
- 🔬 **Wide binaries (Gaia DR4):** No gravitational anomaly at s > 5 AU — discriminant from MOND
- 🔬 **Gravitational waves (LIGO O4):** Longitudinal polarization h_L/h_T ~ 10⁻²
- 🔬 **SMBH spin correlation:** v_max / V_flat must correlate with measured SMBH spin
- 🔬 **Milgrom's a₀:** Not a fundamental constant — emergent from SMBH vortex centripetal acceleration (g_vortex / a₀ = 1.85 at solar radius)

---

## Dark Matter & Dark Energy

CST provides a physical framework for the **95% of the cosmic energy budget** unexplained by ΛCDM:

**Dark Matter (~27%)** → Three physical effects, no new particles:
1. SMBH vortex dragging (dominant at large radii)
2. CST stellar enhancement from ~10⁸ stellar-mass black holes (ε = 0.215)
3. Vortex anisotropy (sin²θ → flat disk morphology)

**Dark Energy (~68%)** → Primordial nucleation residual:
- The Big Bang as incomplete phase transition of the spacetime fluid
- ρ_Λ as thermodynamic residual (resolves fine-tuning)
- ρ_Λ ~ ρ_matter from same transition (resolves coincidence problem)
- w = −1 from tension in incompletely converted fluid

---

## Cosmological Safety

CST preserves all precision cosmological constraints:

| Epoch | Constraint | CST deviation | Limit |
|-------|-----------|---------------|-------|
| BBN (z ~ 10⁹) | ΔG/G | ~10⁻⁴⁹ | < 10⁻⁶ |
| CMB (z = 1100) | Δℓ (peak shift) | < 0.0003 | < 0.1 (Planck) |
| Solar System (z = 0) | dG/dt | ~10⁻¹⁹ yr⁻¹ | < 7×10⁻¹⁴ yr⁻¹ (LLR) |

---

## Repository Structure

```
CST-Compressible-Spacetime-Dynamics/
│
├── README.md
│
├── manuscripts/
│   ├── CST_MANUSCRIPT_ENG_PUBLIC.pdf       # Full manuscript (English, 96 pages)
│   ├── CST_MANUSCRIPT_ENG_PUBLIC.tex       # LaTeX source
│   ├── CST_MANUSCRIPT_ENG_ANONYMOUS.pdf    # Anonymous version for journal submission
│   ├── CST_MANUSCRIPT_ENG_ANONYMOUS.tex    # LaTeX source
│   ├── CST_MANOSCRITTO_ITA.pdf             # Italian version (97 pages)
│   ├── CST_MANOSCRITTO_ITA.tex             # LaTeX source
│   └── CST_BIBLIOGRAPHY.pdf               # Standalone bibliography (76 references)
│
├── figures/
│   ├── 01_exoplanet_validation.png         # Exoplanet fit quality, residuals, cross-validation
│   ├── 02_binary_rv_comparison.png         # Binary star mass ratio analysis
│   ├── 03_theory_diagrams.png              # w(M), H(z), G_eff enhancement
│   ├── 04_cosmology_timeline.png           # CST activation timeline
│   ├── 05_rotation_curves.png              # 8 SPARC galaxies: observed vs CST Model G
│   ├── 06_RAR.png                          # Radial Acceleration Relation: CST vs MOND
│   ├── 07_BTFR.png                         # Baryonic Tully-Fisher Relation
│   ├── 08_model_comparison.png             # Model G diagnostics and comparison
│   └── 09_beta_23.png                      # β = 2/3 geometric invariant
│
├── analysis/
│   ├── cst_sparc_analysis.py               # Main SPARC 129-galaxy analysis
│   ├── cst_refined_model.py                # Model F/G with galactic interference
│   ├── cst_gap_closure.py                  # Cross-validation and gap analysis
│   ├── cst_perturbation_model.py           # Stellar perturbation models (A–J)
│   ├── cst_rar_btfr_test.py                # RAR and BTFR computation
│   ├── cst_scaling_search.py               # Vortex parameter scaling relations
│   └── cst_zero_param_optimization.py      # Zero per-galaxy parameter optimization
│
├── results/
│   ├── CST_SPARC_FULL_RESULTS.txt          # Per-galaxy results (129 galaxies)
│   ├── CST_REFINEMENT_RESULTS.txt          # Model comparison table
│   ├── CST_GAP_CLOSURE_RESULTS.txt         # Cross-validation results
│   ├── CST_GAP_FINAL_DIAGNOSIS.txt         # Gap analysis (profile, braking, spin)
│   ├── CST_BETA_23_DERIVATION.txt          # Full β = 2/3 derivation
│   ├── CST_PERTURBATION_ANALYSIS.txt       # Stellar perturbation analysis
│   └── CST_BRAKING_ANALYSIS.txt            # Vortex braking by stellar mass
│
├── books/                                  # Popular science trilogy (Amazon KDP)
│   ├── Book1_Il_Mistero_di_G/              # "Il Mistero di G" / "The Mystery of G"
│   ├── Book2_Lo_Spaziotempo_Fluido/        # "Lo Spaziotempo Fluido" / "Fluid Spacetime"
│   └── Book3_Prima_del_Tempo/              # "Prima del Tempo" / "Before Time"
│
└── data/
    └── (SPARC data from Lelli et al. 2016 — see http://astroweb.cwru.edu/SPARC/)
```

---

## How to Reproduce

### Requirements
- Python 3.8+
- NumPy, SciPy, Matplotlib
- LaTeX (pdflatex) with standard packages

### Quick Start

```bash
# Clone the repository
git clone https://github.com/antidotiblog-cloud/CST-Compressible-Spacetime-Dynamics.git
cd CST-Compressible-Spacetime-Dynamics

# Download SPARC data (Table1_mrt.txt, Table2_mrt.txt)
# from http://astroweb.cwru.edu/SPARC/

# Run the main SPARC analysis
python analysis/cst_sparc_analysis.py

# Run RAR and BTFR tests
python analysis/cst_rar_btfr_test.py

# Run Model G optimization and cross-validation
python analysis/cst_gap_closure.py

# Compile the manuscript
cd manuscripts
pdflatex CST_MANUSCRIPT_ENG_PUBLIC.tex
```

---

## Popular Science Books

The CST theory is also presented in a trilogy of popular science books, written for a general audience (~14 years reading level), available on Amazon KDP in both Italian and English:

1. **Il Mistero di G / The Mystery of G** — Why gravity might not be constant
2. **Lo Spaziotempo Fluido / Fluid Spacetime** — How spacetime behaves like a fluid
3. **Prima del Tempo / Before Time** — What existed before the Big Bang

---

## Citation

If you use this work, please cite:

```bibtex
@article{Vizzutti2026,
  author  = {Vizzutti, Michele},
  title   = {Compressible Spacetime Dynamics: Observational Evidence for
             Variable Gravitational Coupling from Planetary to Galactic Scales},
  year    = {2026},
  note    = {Manuscript in preparation},
  url     = {https://github.com/antidotiblog-cloud/CST-Compressible-Spacetime-Dynamics}
}
```

---

## Status

| Milestone | Status |
|-----------|--------|
| Theoretical framework | ✅ Complete |
| Exoplanet validation (4,585 systems) | ✅ R² = 96.04% |
| Binary star validation (16,980 systems) | ✅ R² = 96.96% |
| Galactic validation (129 SPARC galaxies) | ✅ RMSE = 20.0 km/s |
| Model G (zero per-galaxy parameters) | ✅ Cross-validated |
| β = 2/3 geometric derivation | ✅ Complete |
| RAR test (beats MOND) | ✅ 0.121 vs 0.163 dex |
| BTFR test | ✅ Slope 3.53 (compatible) |
| Cosmological safety (BBN, CMB) | ✅ Preserved |
| Popular science trilogy | ✅ Published (Amazon KDP) |
| Scientific manuscript | ✅ Ready for submission |
| Journal submission (CQG) | 🔄 In preparation |
| Independent replication | ⏳ Awaiting community |

---

## Contact

Michele Vizzutti — Independent Researcher, Udine, Italy

GitHub: [antidotiblog-cloud](https://github.com/antidotiblog-cloud)

---

## License

This work is shared for scientific review and reproducibility. All code is provided under MIT License. The manuscripts and theoretical content are © Michele Vizzutti 2026.

---

*"Why 2/3? Because we live in 3 dimensions."*
