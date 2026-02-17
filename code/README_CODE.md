# Code Documentation — CST Analysis Pipeline

## Overview

This folder contains the complete Python analysis pipeline for the **Compressible Spacetime Dynamics (CST)** theory validation. All scripts are self-contained and documented.

**Python version:** 3.8+  
**Author:** Michele Vizzutti  
**License:** MIT

---

## Files

### `fit_master.py` ⭐ MAIN SCRIPT
**Purpose:** Unified CST fitting script for all astronomical datasets  
**Function:** Fits the G_eff(M,z) formula to exoplanets and binary stars

**What it does:**
1. Loads NASA exoplanet data (4,585 systems)
2. Loads Gaia DR3 binary data (16,980 systems)
3. Computes cosmological parameters (H(z)/H₀, transition function f(z))
4. Fits CST formula using OLS for exoplanets
5. Fits interference theory Ψ(q,a,M) using non-linear least squares for binaries
6. Runs statistical validation: bootstrap (1,000 iter), K-fold CV (10-fold), AIC/BIC
7. Generates all publication-quality figures (300 DPI)
8. Outputs complete results summary

**Usage:**
```bash
# Full analysis
python fit_master.py

# With figure generation
python fit_master.py --generate-figures

# Exoplanets only
python fit_master.py --dataset exoplanets

# Binaries only
python fit_master.py --dataset binaries

# Verbose mode
python fit_master.py --verbose
```

**Key outputs:**
- `results/fit_summary.txt` — complete statistical results
- `results/parameters.csv` — calibrated parameters with uncertainties
- `figures/` — all publication figures (PNG, 300 DPI)

**Core formula implemented:**
```python
def G_eff(M_solar, z, alpha=0.279, beta=0.685):
    """Effective gravitational constant for compact objects."""
    w = np.exp(-np.abs(M_solar - 1.0))          # weight function
    f = f_transition(z)                           # cosmological transition
    return 1.0 + (1.0 - w) * alpha * M_solar**beta * f

def f_transition(z, z_trans=30, n=3):
    """Transition function f(z) — suppresses G_eff at high z."""
    H_ratio = np.sqrt(0.315*(1+z)**3 + 0.685)    # H(z)/H₀
    S = 1.0 / (1.0 + (z/z_trans)**n)              # structure suppression
    return H_ratio * S
```

---

### `download_nasa_exoplanets.py`
**Purpose:** Download complete exoplanet dataset from NASA Exoplanet Archive  
**Function:** Retrieves confirmed exoplanets with stellar and orbital parameters

**What it does:**
1. Queries NASA TAP API for confirmed planets
2. Filters: requires stellar mass, age, orbital semi-major axis
3. Applies quality cuts (age < 13.8 Gyr, mass > 0, semi-major axis > 0)
4. Computes derived quantities: formation redshift, H(z)/H₀
5. Saves clean dataset to CSV

**Usage:**
```bash
python download_nasa_exoplanets.py
# Output: data/exoplanets_NASA_archive.csv (N ≈ 4,585 systems)
```

**API endpoint:** https://exoplanetarchive.ipac.caltech.edu/TAP/sync

**Columns retrieved:**
- `pl_name`: Planet name
- `st_mass`: Stellar mass [M☉]
- `st_age`: Stellar age [Gyr]
- `pl_orbsmax`: Orbital semi-major axis [AU]
- `pl_orbper`: Orbital period [days]
- `st_met`: Stellar metallicity [Fe/H]
- `st_logg`: Surface gravity [log(cm/s²)]
- `st_rad`: Stellar radius [R☉]

---

### `analyze_alpha_stratified.py`
**Purpose:** Test whether coupling strength α depends on system mass  
**Function:** Stratified analysis to determine if α is constant or α = α(M)

**What it does:**
1. Divides binary dataset into mass bins (0.5-1.0, 1.0-1.5, 1.5-2.0, 2.0+ M☉)
2. Fits α independently for each mass bin
3. Tests for systematic trend: is α(M) = constant?
4. Compares to theory: compact regime α ≈ 0.279 vs extended α_cosmo ≈ 0.07

**Key finding:** α is CONSTANT at 0.279 across all mass bins (p_trend > 0.3).  
The amplification in binaries comes from the interference factor Ψ(q,a,M), NOT from varying α.

**Usage:**
```bash
python analyze_alpha_stratified.py
# Output: results/alpha_stratified_analysis.pdf
```

---

### `verify_r2.py`
**Purpose:** Verify R² values from source data  
**Function:** Critical verification of statistical results reported in manuscript

**What it does:**
1. Loads the Gaia DR3 fit results
2. Recomputes R² from first principles
3. Verifies: R²(Gaia binaries) = 96.96% (NOT 99.38% as in earlier versions)
4. Documents the correction and its origin
5. Cross-checks all reported R² values

**Usage:**
```bash
python verify_r2.py
# Prints verification report with correct vs previously reported values
```

**Critical finding:** The R² = 99.38% reported in the original submission was the value from the *synthetic* dataset validation, not from real Gaia DR3 data. The correct value is 96.96%.

---

### `fit_sparc_galaxies.py`
**Purpose:** Test CST on galaxy rotation curves to calibrate α_cosmo  
**Function:** Extended-structure regime analysis using SPARC database

**What it does:**
1. Downloads/loads SPARC galaxy rotation curve data (175 galaxies)
2. Computes G_eff,extended(z) = G_N × [1 + α_cosmo × f(z)]
3. Fits α_cosmo to minimize residuals on rotation curves
4. Tests hypothesis: α_cosmo ≈ 0.05–0.10 << α = 0.279

**Note:** This script requires SPARC data download from http://astroweb.cwru.edu/SPARC/

**Usage:**
```bash
# First download SPARC data
python fit_sparc_galaxies.py --download-data
# Then run analysis
python fit_sparc_galaxies.py
```

**Expected result:** α_cosmo ≈ 0.05–0.10, confirming extended-structure regime has weaker coupling than compact objects.

---

## Installation

```bash
# Install all dependencies
pip install -r requirements.txt

# Or install manually
pip install numpy>=1.21 pandas>=1.3 scipy>=1.7 matplotlib>=3.4 astropy>=4.3 astroquery>=0.4
```

---

## Directory Structure (after running)

```
code/
├── fit_master.py
├── download_nasa_exoplanets.py
├── analyze_alpha_stratified.py
├── verify_r2.py
├── fit_sparc_galaxies.py
├── requirements.txt
├── README_CODE.md
│
├── data/                          (created on first run)
│   ├── exoplanets_NASA_archive.csv
│   └── binaries_gaia_dr3.csv
│
└── results/                       (created on first run)
    ├── fit_summary.txt
    ├── parameters.csv
    └── validation_report.txt
```

---

## Core Functions Reference

```python
# === COSMOLOGICAL FUNCTIONS ===

def age_to_redshift(age_Gyr, age_universe=13.8):
    """Convert stellar age to formation redshift (matter-dominated approx)."""
    t_form = age_universe - age_Gyr
    return max(0, (age_universe / t_form)**(2/3) - 1)

def hubble_ratio(z, Omega_m=0.315):
    """H(z)/H₀ for flat ΛCDM cosmology."""
    return np.sqrt(Omega_m * (1+z)**3 + (1-Omega_m))

def f_transition(z, z_trans=30, n=3):
    """Transition function f(z). Key property: f(z→∞)→0 (BBN safe)."""
    H_ratio = np.sqrt(0.315*(1+z)**3 + 0.685)
    S = 1.0 / (1.0 + (z/z_trans)**n)
    return H_ratio * S

# === G_eff FUNCTIONS ===

def w_mass(M_solar):
    """Weight function w(M). Key: w(M☉)=1 → G_eff=G_N at solar mass."""
    return np.exp(-np.abs(M_solar - 1.0))

def G_eff_planetary(M_solar, z, alpha=0.279, beta=0.685):
    """G_eff for single-star (planetary) systems.
    CORRECT formula includes [1-w(M)] term!
    """
    w = w_mass(M_solar)
    f = f_transition(z)
    return 1.0 + (1.0 - w) * alpha * (M_solar**beta) * f

def Psi_interference(q, a_AU, M_tot, gamma0=8.0, a0=0.50, beta=0.685):
    """Interference amplification factor for binary stars.
    Ab initio predictions: gamma0=8.0, a0=0.50 AU
    """
    f_q = 4*q / (1+q)**2
    f_a = np.exp(-a_AU / a0)
    f_M = M_tot**beta
    return 1.0 + gamma0 * f_q * f_a * f_M

def G_eff_binary(M_tot, z, q, a_AU, alpha=0.279, **psi_params):
    """G_eff for binary star systems with interference."""
    w = w_mass(M_tot)
    psi = Psi_interference(q, a_AU, M_tot, **psi_params)
    f = f_transition(z)
    return 1.0 + (1.0 - w) * alpha * psi * f

# === OBSERVABLES ===

def xi_ratio(M_solar, z, alpha=0.279, beta=0.685):
    """Observable velocity ratio xi = v_obs/v_Kep = sqrt(G_eff/G_N)."""
    g_ratio = G_eff_planetary(M_solar, z, alpha, beta)
    return np.sqrt(g_ratio)
```

---

## Tests

```bash
# Run unit tests
python -m pytest tests/ -v

# Quick sanity check
python -c "
import numpy as np
# Test: G_eff(M☉, z=0) should be exactly 1.0
from fit_master import G_eff_planetary
g = G_eff_planetary(1.0, 0.0)
assert abs(g - 1.0) < 1e-10, f'FAILED: G_eff(M☉, 0) = {g}'
print(f'G_eff(M☉, z=0) = {g:.10f} ✅ (expected 1.0)')

# Test: G_eff at BBN should be ~1.0
from fit_master import f_transition
f_bbn = f_transition(1e9)
print(f'f(z=10^9) = {f_bbn:.2e} ✅ (should be ~0)')
"
```

---

## Reproducibility

All results in the manuscript can be reproduced by:
1. Running `download_nasa_exoplanets.py` to get fresh data
2. Running `fit_master.py` to reproduce all fits
3. Running `verify_r2.py` to confirm R² values

Computational environment tested on:
- Python 3.9.7, NumPy 1.21.2, SciPy 1.7.1, Matplotlib 3.4.3
- Ubuntu 20.04 / macOS 12.0 / Windows 10

---

*Last updated: 17 February 2026*
