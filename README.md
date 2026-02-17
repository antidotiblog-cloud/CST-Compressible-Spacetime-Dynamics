# Compressible Spacetime Dynamics (CST)
## Observational Evidence for Mass-Dependent Gravitational Coupling

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

**Author:** Michele Vizzutti — Independent Research, Udine, Italy  
**Contact:** antidoti.blog@gmail.com  
**Version:** 2.1 — February 2026

---

## Overview

This repository contains the complete research package for **Compressible Spacetime Dynamics (CST)** theory, presenting observational evidence that the gravitational constant G is not fundamental but emerges from spacetime-matter coupling with mass-dependent and cosmologically-varying strength.

### Key Result

Validation across **21,565 astronomical systems** (4,585 exoplanets + 16,980 Gaia DR3 binary stars) demonstrates:

```
G_eff(M,z) = G_N × {1 + [1-w(M)] × α × (M/M☉)^β × f(z)}

α = 0.279 ± 0.012   (coupling strength)
β = 0.685 ± 0.018   (mass scaling, vs theoretical prediction β = 2/3 at 2.7%)
a₀ = 0.50 ± 0.03 AU (resonance scale, ab initio prediction confirmed)
R² = 97.73%         (multi-scale combined, vs 45.2% for pure Kepler)
```

---

## Repository Structure

```
CST_PACKAGE_v2_1/
│
├── README.md                          ← This file
│
├── manuscripts/                       ← Scientific papers
│   ├── CST_MANOSCRITTO_COMPLETO_ITA.md   ← Full Italian manuscript (v2.1, ~28,000 words)
│   ├── MANUSCRIPT_CST_CORRECTED.tex      ← LaTeX English version (AASTeX format)
│   └── README_MANUSCRIPTS.md             ← Description of manuscript versions
│
├── figures/                           ← Publication-quality figures (300 DPI)
│   ├── 01_exoplanet_validation.png       ← Figure 1: Exoplanet CST validation
│   ├── 02_binary_rv_comparison.png       ← Figure 2: Binary star RV vs photometry
│   ├── 03_theory_diagrams.png            ← Figure 3: CST theoretical framework
│   ├── 04_cosmology_timeline.png         ← Figure 4: Cosmic evolution timeline
│   ├── binaries_gaia_fit_diagnostic.png  ← Figure 5: Gaia DR3 fit diagnostics
│   └── FIGURE_CAPTIONS.md                ← Complete captions (IT + EN)
│
├── code/                              ← Python analysis pipeline
│   ├── fit_master.py                     ← MAIN: unified CST fitting script
│   ├── download_nasa_exoplanets.py       ← Download NASA Exoplanet Archive data
│   ├── fit_sparc_galaxies.py             ← Galaxy rotation curves (α_cosmo)
│   ├── analyze_alpha_stratified.py       ← Test α(M) mass dependence
│   ├── verify_r2.py                      ← Verify R² values from source data
│   ├── requirements.txt                  ← Python dependencies
│   └── README_CODE.md                    ← Code documentation and usage guide
│
└── docs/                              ← Additional documentation
    ├── THEORY_SUMMARY.md                 ← Concise theory reference (IT + EN)
    ├── PARAMETERS_REFERENCE.md           ← All parameters with derivations
    └── REVISION_HISTORY.md               ← Version history and corrections
```

---

## Quick Start

```bash
# 1. Clone repository
git clone https://github.com/michelevizzutti/CST-compressible-spacetime.git
cd CST-compressible-spacetime

# 2. Install dependencies
pip install -r code/requirements.txt

# 3. Download exoplanet data
python code/download_nasa_exoplanets.py

# 4. Run main fit
python code/fit_master.py

# 5. Verify R² values
python code/verify_r2.py
```

---

## Core Formula

**For compact objects (planetary systems and binary stars):**
```
G_eff(M,z) = G_N × {1 + [1-w(M)] × α × (M/M☉)^β × f(z)}

where:
  w(M) = exp(-|M/M☉ - 1|)          weight function
  f(z) = √[Ω_m(1+z)³ + Ω_Λ]       cosmological transition
         ─────────────────────
             1 + (z/30)³
```

**For binary star interference:**
```
G_eff^(bin) = G_N × {1 + [1-w(M)] × α × Ψ(q,a,M) × f(z)}

Ψ(q,a,M) = 1 + γ₀ × M^η × [4q/(1+q)²] × exp(-a/a₀) × M^β
```

**Cosmological safety:** The transition function f(z) guarantees:
- At BBN (z~10⁹): ΔG/G < 10⁻¹⁴ → primordial abundances preserved ✅
- At CMB (z=1100): ΔG/G = 2.76×10⁻⁶ → Planck peaks unchanged ✅

---

## Key Predictions (Testable with Current Facilities)

| Prediction | Facility | Timeline | Signal |
|---|---|---|---|
| Exponential decay at a₀ = 0.50 AU | **Gaia DR4** | 2027 | >10σ |
| Longitudinal GW mode h_L/h_T ~ 0.01 | **LIGO O4** | 2025 | Stack ~200 events |
| Enhanced structure growth f·σ₈ | **Euclid** | 2028 | 5-10% |
| Massive galaxies z>10 | **JWST** | now | Already seen ✅ |
| Orbital period drift | **PLATO** | 2027 | μs/yr |

---

## Statistical Validation

| Test | Result | Threshold | Status |
|---|---|---|---|
| K-fold CV (10-fold) | 0.67% overfitting | < 2% | ✅ |
| Bootstrap CI (1000 iter) | excludes zero | zero excluded | ✅ |
| Metallicity killer test | α changes 3% | < 10% | ✅ |
| AIC vs Kepler | ΔAIC = -22.1 | < -10 "very strong" | ✅ |
| Combined significance | p < 10⁻²⁵⁰ | — | ✅ |

---

## Data Sources

- **Exoplanets:** NASA Exoplanet Archive (https://exoplanetarchive.ipac.caltech.edu/)  
  N = 4,585 confirmed systems with stellar age, mass, orbital parameters
- **Binary Stars:** ESA Gaia DR3 Non-Single Stars catalog  
  N = 16,980 systems with spectroscopic orbits
- **Galaxies (future):** SPARC database (Lelli et al. 2016)  
  N = 175 galaxies with rotation curves

---

## Citation

If you use this work, please cite:

```bibtex
@article{vizzutti2026cst,
  author  = {Vizzutti, Michele},
  title   = {Compressible Spacetime Dynamics: Observational Evidence for 
             Mass-Dependent and Cosmologically-Varying Gravitational Coupling},
  journal = {Monthly Notices of the Royal Astronomical Society},
  year    = {2026},
  note    = {Submitted. arXiv:2602.XXXXX}
}
```

---

## License

- **Code:** MIT License (see `code/LICENSE`)
- **Manuscripts:** Creative Commons CC BY 4.0
- **Figures:** Creative Commons CC BY 4.0

---

## Contact

**Michele Vizzutti**  
Independent Researcher, Udine, Italy  
Email: antidoti.blog@gmail.com

*This research was conducted independently without institutional affiliation. All data used are publicly available.*
