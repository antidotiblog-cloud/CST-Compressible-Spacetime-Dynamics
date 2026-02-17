# Manuscripts — Compressible Spacetime Dynamics

This folder contains the complete scientific manuscripts for the CST theory.

---

## Files

### 1. `CST_MANOSCRITTO_COMPLETO_ITA.md` — **PRIMARY MANUSCRIPT** ⭐
**Language:** Italian  
**Version:** 2.1 (revised 17 February 2026)  
**Length:** ~28,000 words, 3,120 lines  
**Status:** READY FOR SUBMISSION

**Contents:**
- Abstract (~800 words)
- Section 1: Introduction — gravitational constant problem, astrophysical anomalies, previous approaches (~5,300 words)
- Section 2: Theoretical Framework — barotropic spacetime, G_eff formula, β=2/3 derivation, binary interference theory, a₀=0.50 AU ab initio derivation (~5,200 words)
- Section 3: Cosmological Safety — BBN, CMB compatibility, JWST structure tension (~5,000 words)
- Section 4: Data and Methodology — NASA exoplanets, Gaia DR3, synthetic dataset, statistical methods (~4,200 words)
- Section 5: Results — R²=96.04% (exoplanets), R²=96.96% (Gaia), R²=97.73% (combined), killer tests (~3,800 words)
- Section 6: Discussion — ab initio agreements, dark matter implications, gravitational waves, cyclic cosmology (~3,500 words)
- Section 7: Observational Predictions — 10 falsifiable predictions with Gaia DR4, LIGO O4, Euclid, JWST, SKA, PLATO, Vera Rubin, CMB-S4, Einstein Telescope (~3,400 words)
- Section 8: Conclusions (~2,500 words)
- Appendix A: Mathematical derivations (barotropic equations, a₀ derivation without circularity, w(M) differentiability)
- Appendix B: Stellar resonance and cosmological modes (why M☉ is special)
- Appendix C: Numerical tables (parameters, f(z) values, comparison table)
- Bibliography: 52 references (MNRAS format, A-Z order)

**Corrections applied (v2.0 → v2.1, 17/02/2026):**
1. Fixed critical formula error in Abstract: added missing `[1-w(M)]` term
2. Removed internal revision history references
3. Fixed number format (Italian: period as thousands separator)
4. Fixed author initials: M.N. → M.V.
5. Added missing bibliography entry: Quinn et al. (2013), PRL 111, 101102
6. Fixed spelling: "Exoplaneti" → "Esopianeti"
7. Clarified α uncertainty: ±0.021 (single OLS) vs ±0.012 (weighted multi-dataset)
8. Notation consistency: β_theory → β_teo
9. Removed "(Corretta)" from version header

---

### 2. `MANUSCRIPT_CST_CORRECTED.tex` — English LaTeX Version
**Language:** English  
**Format:** AASTeX 6.3.1 (compatible with ApJ, AJ)  
**Version:** Corrected (February 2026)  
**Status:** Earlier version — R² values and formula updated vs original rejection version

**Contents:**
- Abstract with corrected R² values (96.96% for Gaia, not 99.38%)
- Introduction with astrophysical motivation
- Theoretical framework derivation
- Observational validation (exoplanets + binaries)
- Predictions and conclusions
- Full bibliography in BibTeX format

**Note:** This LaTeX file corresponds to an earlier stage of the manuscript. The Italian `.md` file contains the most complete and up-to-date version of the theory. For submission, the Italian manuscript should be translated/adapted to this LaTeX template.

**Important:** The LaTeX abstract formula `G_eff = G_N[1 + α(M/M☉)^β × H(z)/H₀]` is an earlier form. The correct complete formula includes the weight function term `[1-w(M)]` and the transition function `f(z)`, as detailed in the Italian manuscript.

---

## Formula Reference (Correct Version)

The scientifically correct formula, as in the Italian manuscript v2.1:

```
G_eff(M,z) = G_N × {1 + [1-w(M)] × α × (M/M☉)^β × f(z)}

w(M) = exp(-|M/M☉ - 1|)                          weight function
f(z) = √[Ω_m(1+z)³ + Ω_Λ] / [1 + (z/30)³]      transition function
α = 0.279 ± 0.012                                  coupling strength
β = 0.685 ± 0.018                                  mass scaling
```

**WITHOUT** `[1-w(M)]`, the formula incorrectly predicts G_eff(M☉) = 1.279 G_N  
**WITH** `[1-w(M)]`, correctly: G_eff(M☉) = G_N (because w(M☉)=1 → [1-w]=0) ✅

---

## Suggested Journals for Submission

1. **MNRAS** (Monthly Notices of the Royal Astronomical Society) — recommended for multi-scale astronomical data work
2. **ApJ** (The Astrophysical Journal) — top-tier, strong for N=21,565 system dataset
3. **CQG** (Classical and Quantum Gravity) — for theoretical gravity aspects

---

*Last updated: 17 February 2026*
