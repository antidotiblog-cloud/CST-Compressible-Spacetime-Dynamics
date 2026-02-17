# CST Theory Summary — Reference Card
## Riassunto della Teoria CST — Scheda di Riferimento

**Version / Versione:** 2.1 — February / Febbraio 2026

---

## Core Formula / Formula Principale

```
G_eff(M,z) = G_N × {1 + [1-w(M)] × α × (M/M☉)^β × f(z)}
```

| Component / Componente | Formula | Value / Valore |
|---|---|---|
| Weight function / Funzione peso | `w(M) = exp(-\|M/M☉ - 1\|)` | w(M☉) = 1, w(0)→0 |
| Coupling / Accoppiamento | `α` | 0.279 ± 0.012 |
| Mass scaling / Scaling massa | `β` | 0.685 ± 0.018 |
| Transition / Transizione | `f(z) = H(z)/H₀ / [1+(z/30)³]` | f(0)=1, f(∞)→0 |

**Critical:** Without `[1-w(M)]`, G_eff(M☉) = 1.279 G_N (WRONG). With it: G_eff(M☉) = G_N ✅

---

## Binary Interference / Interferenza Binarie

```
G_eff^(bin) = G_N × {1 + [1-w(M)] × α × Ψ(q,a,M) × f(z)}
Ψ(q,a,M) = 1 + γ₀ × M^η × [4q/(1+q)²] × exp(-a/a₀) × M^β
```

| Parameter / Parametro | Ab initio prediction / Predizione | Observed / Osservato |
|---|---|---|
| γ₀ (intensity / intensità) | 8.0 | 8.3 ± 0.8 |
| a₀ (resonance / risonanza) | ~0.50 AU | 0.50 ± 0.03 AU |
| β (mass scaling) | 2/3 = 0.667 | 0.685 ± 0.018 |

---

## Results Summary / Riepilogo Risultati

| Dataset | N | R² | Note |
|---|---|---|---|
| NASA exoplanets | 4,585 | **96.04%** | OLS + bootstrap |
| Gaia DR3 binaries | 16,980 | **96.96%** | ← CORRECT value |
| Synthetic validation | 6,744 | **99.19%** | Theory self-test |
| **Multi-scale combined** | **21,565** | **97.73%** | Main result |
| Kepler only (baseline) | 21,565 | 45.2% | Comparison |

---

## Cosmological Safety / Sicurezza Cosmologica

| Epoch / Epoca | z | ΔG/G | Status |
|---|---|---|---|
| Today / Oggi | 0 | 15–28% | Self-consistent ✅ |
| CMB | 1100 | 2.76×10⁻⁶ | SAFE ✅ |
| BBN | ~10⁹ | < 10⁻¹⁴ | SAFE ✅ |

---

## Key Physical Interpretation / Interpretazione Fisica Chiave

**English:** Spacetime behaves as a compressible fluid (equation of state P_ST = c_s² ρ_ST). Matter induces local compression waves. In binary systems, two masses create overlapping waves with constructive interference peaked at a₀ = 0.50 AU. Systems formed during rapid cosmic expansion (high H(z)/H₀) crystallize stronger gravitational coupling through a "lock-in" mechanism.

**Italiano:** Lo spaziotempo si comporta come un fluido compressibile (equazione di stato P_ST = c_s² ρ_ST). La materia induce onde di compressione locali. Nei sistemi binari, due masse creano onde sovrapposte con interferenza costruttiva con picco a a₀ = 0,50 AU. I sistemi formati durante rapida espansione cosmica (H(z)/H₀ alto) cristallizzano accoppiamento gravitazionale più forte attraverso un meccanismo di "lock-in".

---

## Why M☉ Is Special / Perché M☉ È Speciale

Stellar p-mode frequencies at M = M☉ resonate with cosmological mode frequencies at H₀:

```
ω_star(M☉) ≈ ω_cosmo = H₀
→ Geometric resonance → w(M☉) = 1 → G_eff(M☉) = G_N
```

This is NOT an arbitrary calibration — it follows from the resonance condition.

---

## Testable Predictions / Predizioni Verificabili

| Prediction / Predizione | Facility / Strumento | Timeline | Signal / Segnale |
|---|---|---|---|
| Exponential decay ∝ exp(-a/0.5 AU) | Gaia DR4 | 2027 | >10σ |
| Longitudinal GW mode h_L/h_T ~ 0.01 | LIGO O4 | 2025 | Stack 200 events |
| Enhanced structure growth | Euclid | 2028 | 5-10% |
| Massive galaxies z>10 | JWST | now | Partly observed ✅ |

---

*Last updated / Ultimo aggiornamento: 17 February / Febbraio 2026*
