import numpy as np
from scipy.optimize import minimize

PDMF_enh = 0.215; ML_disk = 0.5; ML_bul = 0.7

def v_LO(r, vm, rc):
    x = r / rc
    if x < 1e-6: return 0.0
    return vm * (1 - np.exp(-x**2)) / x / 0.6382
v_LO_v = np.vectorize(v_LO)

# Parse data (same as before)
with open('/mnt/user-data/uploads/Table1_mrt.txt','r') as f: lines1 = f.readlines()
last_sep = max(i for i,l in enumerate(lines1) if l.strip().startswith('---'))
t1 = {}
for i in range(last_sep+1, len(lines1)):
    p = lines1[i].split()
    if len(p) < 18: continue
    try: t1[p[0]] = {'T':int(p[1]),'D':float(p[2]),'Rd':float(p[11]),'Vf':float(p[15]),
                      'Q':int(p[17]),'L':float(p[7]),'MHI':float(p[13]),'RHI':float(p[14])}
    except: continue

t2 = {}
with open('/mnt/user-data/uploads/Table2_mrt.txt','r') as f:
    for line in f:
        p = line.split()
        if len(p) < 9: continue
        try:
            n=p[0]; float(p[1])
            if n not in t2: t2[n]={'R':[],'Vo':[],'eV':[],'Vg':[],'Vd':[],'Vb':[]}
            t2[n]['R'].append(float(p[2])); t2[n]['Vo'].append(float(p[3]))
            t2[n]['eV'].append(float(p[4])); t2[n]['Vg'].append(float(p[5]))
            t2[n]['Vd'].append(float(p[6])); t2[n]['Vb'].append(float(p[7]))
        except: continue
for g in t2:
    for k in t2[g]: t2[g][k] = np.array(t2[g][k])

good = [g for g in t1 if g in t2 and t1[g]['Q']<=2 and len(t2[g]['R'])>=5 and t1[g]['Vf']>0]

# ============================================================
# APPROACH: FIT GLOBAL PARAMETERS, ZERO FREE PER GALAXY
# ============================================================
# Model: v_max = A * V_flat^alpha_v
#         r_c  = B * V_flat^alpha_r * R_disk^beta_r
# 
# Global parameters: A, alpha_v, B, alpha_r, beta_r (+ PDMF_enh)
# Per-galaxy parameters: ZERO

print("="*70)
print("OTTIMIZZAZIONE GLOBALE: ZERO PARAMETRI PER GALASSIA")
print("="*70)

def global_chi2(params):
    """Compute total chi2 across ALL galaxies with global params only"""
    A, alpha_v, B, alpha_r, beta_r = params
    
    if A < 0 or B < 0 or alpha_v < 0 or alpha_r < -2 or beta_r < -2:
        return 1e20
    
    total_chi2 = 0
    n_total = 0
    
    for gal in good:
        g1 = t1[gal]; g2 = t2[gal]
        Vf = g1['Vf']; Rd = max(g1['Rd'], 0.01)
        
        # Predict vortex parameters from global formula
        vm_pred = A * Vf**alpha_v
        rc_pred = B * Vf**alpha_r * Rd**beta_r
        rc_pred = max(rc_pred, 0.05)
        
        # Baryonic
        vg2 = np.sign(g2['Vg']) * g2['Vg']**2
        vbar2 = vg2 + ML_disk*g2['Vd']**2 + ML_bul*g2['Vb']**2
        vbar_cst = np.sqrt(np.maximum(vbar2, 0)) * np.sqrt(1+PDMF_enh)
        
        # Vortex
        vv = v_LO_v(g2['R'], vm_pred, rc_pred)
        vpred = np.sqrt(vbar_cst**2 + vv**2)
        
        # Chi2
        w = 1.0 / np.maximum(g2['eV'], 1)**2
        total_chi2 += np.sum(w * (vpred - g2['Vo'])**2)
        n_total += len(g2['R'])
    
    return total_chi2

# Grid search for initial guess
print("\nGrid search per punto iniziale...")
best = (1e30, None)
for A in [0.7, 0.8, 0.9, 1.0, 1.1]:
    for av in [0.9, 1.0, 1.1]:
        for B in [0.01, 0.02, 0.05, 0.1, 0.2]:
            for ar in [1.0, 1.3, 1.5, 1.8]:
                for br in [0.0, 0.3, 0.5, 0.8]:
                    c2 = global_chi2([A, av, B, ar, br])
                    if c2 < best[0]:
                        best = (c2, [A, av, B, ar, br])

print(f"Miglior punto iniziale: {best[1]}, chi2={best[0]:.0f}")

# Refine with Nelder-Mead
res = minimize(global_chi2, best[1], method='Nelder-Mead', 
               options={'maxiter': 10000, 'xatol': 1e-6, 'fatol': 1})
A_f, av_f, B_f, ar_f, br_f = res.x

print(f"\nPARAMETRI GLOBALI OTTIMIZZATI:")
print(f"  v_max = {A_f:.4f} × V_flat^{av_f:.3f}")
print(f"  r_c   = {B_f:.4f} × V_flat^{ar_f:.3f} × R_disk^{br_f:.3f}")

# Compute RMSE with optimized global params
tot_ss = 0; tot_ss_bar = 0; tot_n = 0
per_gal_results = {}

for gal in good:
    g1 = t1[gal]; g2 = t2[gal]
    Vf = g1['Vf']; Rd = max(g1['Rd'], 0.01)
    
    vm_pred = A_f * Vf**av_f
    rc_pred = max(B_f * Vf**ar_f * Rd**br_f, 0.05)
    
    vg2 = np.sign(g2['Vg']) * g2['Vg']**2
    vbar2 = vg2 + ML_disk*g2['Vd']**2 + ML_bul*g2['Vb']**2
    vbar = np.sqrt(np.maximum(vbar2, 0))
    vbar_cst = vbar * np.sqrt(1+PDMF_enh)
    vv = v_LO_v(g2['R'], vm_pred, rc_pred)
    vpred = np.sqrt(vbar_cst**2 + vv**2)
    
    rmse = np.sqrt(np.mean((g2['Vo'] - vpred)**2))
    rmse_bar = np.sqrt(np.mean((g2['Vo'] - vbar)**2))
    
    tot_ss += np.sum((g2['Vo'] - vpred)**2)
    tot_ss_bar += np.sum((g2['Vo'] - vbar)**2)
    tot_n += len(g2['R'])
    per_gal_results[gal] = {'rmse': rmse, 'rmse_bar': rmse_bar, 'Vf': Vf}

rmse_global = np.sqrt(tot_ss / tot_n)
rmse_bar_global = np.sqrt(tot_ss_bar / tot_n)

print(f"\n{'='*60}")
print(f"RISULTATI: MODELLO A ZERO PARAMETRI LIBERI PER GALASSIA")
print(f"{'='*60}")
print(f"  Parametri GLOBALI: {5} (A, α_v, B, α_r, β_r)")
print(f"  Parametri per galassia: 0")
print(f"  Punti dati: {tot_n}")
print(f"")
print(f"  RMSE globale CST (zero-param) = {rmse_global:.2f} km/s")
print(f"  RMSE globale solo barioni     = {rmse_bar_global:.2f} km/s")
print(f"  Miglioramento: {rmse_bar_global/rmse_global:.2f}x")
print(f"")
print(f"  Confronto:")
print(f"    Solo barioni:           {rmse_bar_global:.2f} km/s")
print(f"    CST zero-param (questo): {rmse_global:.2f} km/s ({rmse_bar_global/rmse_global:.2f}x)")
print(f"    CST 2-param/galassia:    18.44 km/s (3.36x)")

# Also try: v_max = V_flat (fixed, no free param) and optimize only r_c formula
print(f"\n\n{'='*60}")
print(f"VARIANTE: v_max = V_flat (esatto), solo r_c ottimizzato")
print(f"{'='*60}")

def global_chi2_simple(params):
    B, ar, br = params
    if B < 0: return 1e20
    total = 0
    for gal in good:
        g1 = t1[gal]; g2 = t2[gal]
        Vf = g1['Vf']; Rd = max(g1['Rd'], 0.01)
        vm_pred = Vf  # EXACT: v_max = V_flat
        rc_pred = max(B * Vf**ar * Rd**br, 0.05)
        vg2 = np.sign(g2['Vg']) * g2['Vg']**2
        vbar2 = vg2 + ML_disk*g2['Vd']**2 + ML_bul*g2['Vb']**2
        vbar_cst = np.sqrt(np.maximum(vbar2, 0)) * np.sqrt(1+PDMF_enh)
        vv = v_LO_v(g2['R'], vm_pred, rc_pred)
        vpred = np.sqrt(vbar_cst**2 + vv**2)
        w = 1.0 / np.maximum(g2['eV'], 1)**2
        total += np.sum(w * (vpred - g2['Vo'])**2)
    return total

best2 = (1e30, None)
for B in [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]:
    for ar in [0.5, 1.0, 1.3, 1.5, 2.0]:
        for br in [0.0, 0.3, 0.5, 0.7, 1.0]:
            c2 = global_chi2_simple([B, ar, br])
            if c2 < best2[0]: best2 = (c2, [B, ar, br])

res2 = minimize(global_chi2_simple, best2[1], method='Nelder-Mead',
                options={'maxiter': 10000})
B2, ar2, br2 = res2.x

tot_ss2 = 0; tot_n2 = 0
for gal in good:
    g1 = t1[gal]; g2 = t2[gal]
    Vf = g1['Vf']; Rd = max(g1['Rd'], 0.01)
    vm_pred = Vf
    rc_pred = max(B2 * Vf**ar2 * Rd**br2, 0.05)
    vg2 = np.sign(g2['Vg']) * g2['Vg']**2
    vbar2 = vg2 + ML_disk*g2['Vd']**2 + ML_bul*g2['Vb']**2
    vbar_cst = np.sqrt(np.maximum(vbar2, 0)) * np.sqrt(1+PDMF_enh)
    vv = v_LO_v(g2['R'], vm_pred, rc_pred)
    vpred = np.sqrt(vbar_cst**2 + vv**2)
    tot_ss2 += np.sum((g2['Vo'] - vpred)**2)
    tot_n2 += len(g2['R'])

rmse2 = np.sqrt(tot_ss2 / tot_n2)
print(f"  v_max = V_flat (fisso)")
print(f"  r_c   = {B2:.4f} × V_flat^{ar2:.3f} × R_disk^{br2:.3f}")
print(f"  Parametri globali: 3 (B, α_r, β_r)")
print(f"  RMSE = {rmse2:.2f} km/s")
print(f"  Miglioramento: {rmse_bar_global/rmse2:.2f}x")

# ULTIMATE: v_max = V_flat, r_c = C * R_disk (just 1 free param!)
print(f"\n\n{'='*60}")
print(f"VERSIONE MINIMALE: v_max = V_flat, r_c = C × R_disk")
print(f"(UN SOLO parametro globale!)")
print(f"{'='*60}")

def chi2_1param(C):
    C = C[0]
    if C < 0: return 1e20
    total = 0
    for gal in good:
        g1 = t1[gal]; g2 = t2[gal]
        Vf = g1['Vf']; Rd = max(g1['Rd'], 0.01)
        vm_pred = Vf
        rc_pred = max(C * Rd, 0.05)
        vg2 = np.sign(g2['Vg']) * g2['Vg']**2
        vbar2 = vg2 + ML_disk*g2['Vd']**2 + ML_bul*g2['Vb']**2
        vbar_cst = np.sqrt(np.maximum(vbar2, 0)) * np.sqrt(1+PDMF_enh)
        vv = v_LO_v(g2['R'], vm_pred, rc_pred)
        vpred = np.sqrt(vbar_cst**2 + vv**2)
        w = 1.0 / np.maximum(g2['eV'], 1)**2
        total += np.sum(w * (vpred - g2['Vo'])**2)
    return total

# Scan C
best_C = min([(chi2_1param([c]), c) for c in np.arange(0.5, 10, 0.1)])
res3 = minimize(chi2_1param, [best_C[1]], method='Nelder-Mead')
C_opt = res3.x[0]

tot_ss3 = 0; tot_n3 = 0
for gal in good:
    g1 = t1[gal]; g2 = t2[gal]
    Vf = g1['Vf']; Rd = max(g1['Rd'], 0.01)
    vm_pred = Vf
    rc_pred = max(C_opt * Rd, 0.05)
    vg2 = np.sign(g2['Vg']) * g2['Vg']**2
    vbar2 = vg2 + ML_disk*g2['Vd']**2 + ML_bul*g2['Vb']**2
    vbar_cst = np.sqrt(np.maximum(vbar2, 0)) * np.sqrt(1+PDMF_enh)
    vv = v_LO_v(g2['R'], vm_pred, rc_pred)
    vpred = np.sqrt(vbar_cst**2 + vv**2)
    tot_ss3 += np.sum((g2['Vo'] - vpred)**2)
    tot_n3 += len(g2['R'])

rmse3 = np.sqrt(tot_ss3 / tot_n3)
print(f"  v_max = V_flat")
print(f"  r_c   = {C_opt:.2f} × R_disk")
print(f"  Parametro globale UNICO: C = {C_opt:.2f}")
print(f"  RMSE = {rmse3:.2f} km/s")
print(f"  Miglioramento: {rmse_bar_global/rmse3:.2f}x")

print(f"\n\n{'='*70}")
print(f"RIEPILOGO COMPARATIVO")
print(f"{'='*70}")
print(f"  {'Modello':40s} {'Param':>6s} {'RMSE':>8s} {'Migliora':>8s}")
print(f"  {'─'*65}")
print(f"  {'Solo barioni':40s} {'0':>6s} {rmse_bar_global:8.2f} {'1.00x':>8s}")
print(f"  {'MOND (a₀ fisso)':40s} {'1':>6s} {'~25':>8s} {'~2.5x':>8s}")
print(f"  {'CST 1-param (v_max=Vf, r_c=C×Rd)':40s} {'1':>6s} {rmse3:8.2f} {f'{rmse_bar_global/rmse3:.2f}x':>8s}")
print(f"  {'CST 3-param (v_max=Vf, r_c=B×Vf^a×Rd^b)':40s} {'3':>6s} {rmse2:8.2f} {f'{rmse_bar_global/rmse2:.2f}x':>8s}")
print(f"  {'CST 5-param (v_max e r_c globali)':40s} {'5':>6s} {rmse_global:8.2f} {f'{rmse_bar_global/rmse_global:.2f}x':>8s}")
print(f"  {'CST 2-param/galassia (258 totali)':40s} {'258':>6s} {'18.44':>8s} {'3.36x':>8s}")
print(f"  {'NFW 2-param/galassia (258 totali)':40s} {'258':>6s} {'~15-20':>8s} {'~3-4x':>8s}")

