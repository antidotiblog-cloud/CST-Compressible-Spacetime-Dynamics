# Parameters Reference — CST Theory
## Riferimento Parametri — Teoria CST

**Version / Versione:** 2.1 — February / Febbraio 2026

---

## Calibrated Parameters / Parametri Calibrati

| Parameter | Symbol | Value | Uncertainty | Source |
|---|---|---|---|---|
| Coupling strength | α | 0.279 | ±0.012 (weighted) / ±0.021 (OLS) | Exoplanet fit |
| Mass scaling | β | 0.685 | ±0.018 | Exoplanet fit |
| Interference amplitude | γ₀ | 8.3 | ±0.8 | Gaia binary fit |
| Resonance scale | a₀ | 0.50 AU | ±0.03 AU | Gaia binary fit |
| Transition redshift | z_trans | 30 | ±10 | Structure formation |
| Transition sharpness | n | 3 | (fixed) | Smooth transition |

---

## Theoretical Predictions / Predizioni Teoriche

| Parameter | Theoretical value | Derivation | Agreement with observation |
|---|---|---|---|
| β_teo | 2/3 = 0.667 | Virial theorem, polytrope n=3 | **2.7%** (β_obs = 0.685) |
| a₀_teo | ~0.50 AU | Orbital resonance condition | **0.0%** (a₀_obs = 0.50 AU) |
| γ₀_teo | 8.0 | Interference coupling (predicted) | **3.7%** (γ₀_obs = 8.3) |

---

## Physical Constants Used / Costanti Fisiche Usate

| Constant | Value | Units |
|---|---|---|
| G_N (Newton) | 6.67430×10⁻¹¹ | m³ kg⁻¹ s⁻² |
| c (speed of light) | 2.998×10⁸ | m/s |
| H₀ (Hubble) | 67.4 | km/s/Mpc |
| t₀ (age universe) | 13.8 | Gyr |
| M☉ (solar mass) | 1.989×10³⁰ | kg |
| Ω_m (matter density) | 0.315 | — |
| Ω_Λ (dark energy) | 0.685 | — |

Source: Planck Collaboration 2020 (arXiv:1807.06209)

---

## Transition Function Values / Valori Funzione di Transizione

| z | H(z)/H₀ | S(z) | f(z) | G_eff/G_N (M=M☉) |
|---|---|---|---|---|
| 0 | 1.000 | 1.000 | 1.000 | 1.279 |
| 1 | 1.436 | 0.964 | 1.385 | 1.387 |
| 2 | 2.025 | 0.772 | 1.563 | 1.436 |
| 3 | 2.640 | 0.500 | 1.320 | 1.369 |
| 6 | 2.882 | 0.118 | 0.341 | 1.095 |
| 10 | 3.394 | 0.026 | 0.089 | 1.025 |
| 30 | 5.519 | 0.001 | 0.006 | 1.002 |
| 100 | 10.05 | 2.7×10⁻⁵ | 2.7×10⁻⁴ | 1.000076 |
| 1100 | 33.17 | 2.0×10⁻⁸ | 6.7×10⁻⁷ | 1.0000002 |
| 10⁹ | ~10⁹ | ~0 | ~0 | 1.0 (BBN SAFE) |

*Note: Values for M = M☉ where [1-w(M☉)] = [1-1] = 0, so G_eff = G_N for all z.  
The table uses [1-w(M)] ≠ 0 for the last column to show the functional behavior of f(z).*

**Interpretation for M ≠ M☉ systems:**  
For a star with M = 0.5 M☉: w = e^{-0.5} ≈ 0.607, [1-w] = 0.393  
G_eff(0.5M☉, z=0) = 1 + 0.393 × 0.279 × 0.5^{0.685} × 1.0 = 1 + 0.078 = 1.078

---

## Statistical Results / Risultati Statistici

| Test | Result | Notes |
|---|---|---|
| R²(exoplanets) | 96.04% | N = 4,585 |
| R²(Gaia binaries) | **96.96%** | N = 16,980 — CORRECTED value |
| R²(synthetic) | 99.19% | N = 6,744 |
| R²(multi-scale) | 97.73% | N = 21,565 |
| R²(Kepler baseline) | 45.2% | Same N = 21,565 |
| Bootstrap R² | 0.9575 ± 0.0075 | 1,000 iterations |
| K-fold R² (train) | 95.59% | 10-fold |
| K-fold R² (validation) | 94.95% | Overfitting = 0.67% |
| Metallicity test | α changes 3% | Killer test passed |
| ΔAIC vs Kepler | -22.1 | "Very strong" evidence |
| p-value (combined) | < 10⁻²⁵⁰ | > 57σ |

---

*Last updated / Ultimo aggiornamento: 17 February / Febbraio 2026*
