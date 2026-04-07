# CST — Sessione di Lavoro 1 Aprile 2026

## Riepilogo dei Risultati

### 1. Test RAR (Radial Acceleration Relation)
- CST scatter: **0.121 dex** — MIGLIORE di MOND (0.163 dex)
- Senza parametro a₀: la RAR emerge dal vortice naturalmente

### 2. Test BTFR (Baryonic Tully-Fisher Relation)
- Pendenza CST: 3.53 vs osservata 3.43 — compatibile

### 3. Origine fisica di a₀ (Milgrom)
- g_vortex(R_sun) = 2.22 × 10⁻¹⁰ m/s²
- g_vortex / a₀ = 1.85
- **a₀ non è una costante fondamentale — è l'accelerazione centripeta del vortice**

### 4. Wide Binaries
- CST predice NESSUNA deviazione a s > 5 AU (discriminante da MOND)

### 5. Scoperta: v_max = V_flat (r = 0.94)
- Il vortice È la curva piatta — zero parametri per v_max

### 6. Scoperta: r_c ∝ M_BH^0.646 ≈ M_BH^(2/3)
- L'esponente 2/3 appare per la TERZA volta nella teoria

### 7. Modello a zero parametri per galassia

| Modello | Param glob | RMSE | Miglior. |
|---------|-----------|------|----------|
| Solo barioni | 0 | 61.9 | 1.0× |
| MOND (a₀) | 1 | ~25 | ~2.5× |
| CST semplice (r_c = C×Rd) | 1 | 25.6 | 2.4× |
| CST perturbazione locale (A) | 2 | 25.3 | 2.5× |
| CST soglia peso (C) | 3 | 25.2 | 2.5× |
| CST interferenza completa (B) | 4 | 22.2 | 2.8× |
| CST fit individuale | 258 | 18.4 | 3.4× |

### 8. Formula galattica CST con perturbazione stellare

```
v²_tot(r) = v²_bar(r) × (1 + ε_PDMF) + v²_vortice(r)

v_vortice(r) = V_flat × LO(r / r_c_eff(r))
r_c_eff(r) = C × R_disk × [1 + η × w(r)]
w(r) = 1 - exp(-(v_bar(r) / v_th)²)

Parametri fissi:
  ε_PDMF = 0.215 (dalla popolazione stellare)
  
Parametri globali (3):
  C = 4.42 (core del vortice in unità di R_disk)
  η = 0.81 (intensità perturbazione stellare)  
  v_th = 72 km/s (soglia di attivazione)
```

### 9. Intuizione di Michele confermata
- L'errore del modello correla con la densità centrale (r = 0.64)
- La soglia 72 km/s separa nane (perturbazione OFF) da massive (ON)
- Stessa struttura dell'interferenza binaria

### 10. Prossimi passi — Raffinamento
- Il gap 22→18 km/s (modello B vs fit individuale) va colmato
- Il modello B suggerisce che la perturbazione è radiale, non legata a v_bar locale
- Esplorare combinazione radiale + densità locale
- Esplorare effetto dello spin del SMBH sullo scatter residuo

---

## Risultati del Raffinamento (sessione pomeriggio)

### Modello F: Interferenza Galattica Completa
**Il modello di Michele** — stessa struttura matematica delle binarie

```
v²_tot(r) = v²_bar(r)(1+ε_PDMF) + v²_vortice(r) × Ψ_gal(r)

Ψ_gal(r) = 1 + γ_gal × [4q/(1+q)²]^p × exp(-r/r₀Rd) × (M_BH/M_ref)^β_gal

q(r) = M_bar(<r) / M_BH

v_max = V_flat
r_c = 19.7 × R_disk

Parametri globali (5):
  γ_gal = 110.2 (intensità interferenza)
  p = 0.24 (esponente massa relativa)
  r₀ = 2.47 (scala radiale in unità di R_disk)
  β_gal = 0.26 (dipendenza da M_BH)
  C = 19.7 (core del vortice)
```

### Risultato: RMSE = 19.8 km/s (3.13× miglioramento)

Con **5 parametri globali** (zero per galassia), CST raggiunge 
prestazioni quasi identiche al fit individuale (18.4 km/s con 258 parametri).

### Confronto finale

| Modello | Param | RMSE | Migl. |
|---------|-------|------|-------|
| Solo barioni | 0 | 61.9 | 1.0× |
| MOND (a₀) | 1 | ~25 | ~2.5× |
| CST 1-param | 1 | 25.6 | 2.4× |
| CST 3-param | 3 | 24.0 | 2.6× |
| **CST interf-F** | **5** | **19.8** | **3.1×** |
| CST individuale | 258 | 18.4 | 3.4× |

### Significato
- Il gap 22→19.8 km/s è stato chiuso con il modello di interferenza
- La stessa struttura matematica (4q/(1+q)², exp(-r/r₀)) funziona
  sia per le binarie stellari che per le galassie
- L'unità concettuale della teoria è preservata: un solo meccanismo
  (interferenza nel fluido spaziotempo) opera a tutte le scale
