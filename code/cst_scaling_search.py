import numpy as np
from scipy.optimize import minimize

PDMF_enh = 0.215; ML_disk = 0.5; ML_bul = 0.7

def v_LO(r, vm, rc):
    x = r / rc
    if x < 1e-6: return 0.0
    return vm * (1 - np.exp(-x**2)) / x / 0.6382
v_LO_v = np.vectorize(v_LO)

# Parse data
with open('/mnt/user-data/uploads/Table1_mrt.txt','r') as f: lines1 = f.readlines()
last_sep = max(i for i,l in enumerate(lines1) if l.strip().startswith('---'))
t1 = {}
for i in range(last_sep+1, len(lines1)):
    p = lines1[i].split()
    if len(p) < 18: continue
    try: t1[p[0]] = {'T':int(p[1]),'D':float(p[2]),'Rd':float(p[11]),'Vf':float(p[15]),
                      'Q':int(p[17]),'L':float(p[7]),'Reff':float(p[9]),'MHI':float(p[13]),
                      'RHI':float(p[14]),'SBdisk':float(p[12]),'SBeff':float(p[10])}
    except: continue

t2 = {}
with open('/mnt/user-data/uploads/Table2_mrt.txt','r') as f:
    for line in f:
        p = line.split()
        if len(p) < 9: continue
        try:
            n=p[0]; float(p[1]); R=float(p[2]); Vo=float(p[3]); eV=float(p[4])
            Vg=float(p[5]); Vd=float(p[6]); Vb=float(p[7])
            if n not in t2: t2[n]={'R':[],'Vo':[],'eV':[],'Vg':[],'Vd':[],'Vb':[]}
            t2[n]['R'].append(R); t2[n]['Vo'].append(Vo); t2[n]['eV'].append(eV)
            t2[n]['Vg'].append(Vg); t2[n]['Vd'].append(Vd); t2[n]['Vb'].append(Vb)
        except: continue
for g in t2:
    for k in t2[g]: t2[g][k] = np.array(t2[g][k])

good = [g for g in t1 if g in t2 and t1[g]['Q']<=2 and len(t2[g]['R'])>=5 and t1[g]['Vf']>0]

# Fit vortex for each galaxy
results = {}
for gal in good:
    g1=t1[gal]; g2=t2[gal]; n=len(g2['R'])
    vg2=np.sign(g2['Vg'])*g2['Vg']**2
    vbar2=vg2 + ML_disk*g2['Vd']**2 + ML_bul*g2['Vb']**2
    vbar = np.sqrt(np.maximum(vbar2, 0))
    vbar_cst = vbar * np.sqrt(1+PDMF_enh)
    
    def chi2(p):
        vm,rc = p
        if vm<0 or rc<0.05: return 1e15
        vv = v_LO_v(g2['R'], vm, rc)
        vpred = np.sqrt(vbar_cst**2 + vv**2)
        return np.sum(((vpred-g2['Vo'])/np.maximum(g2['eV'],1))**2)
    
    rmax=max(g2['R']); vmax_g=max(g2['Vo'])
    best=(1e20,(vmax_g*0.5, rmax*0.5))
    for vm in np.arange(5, vmax_g*1.5, max(vmax_g/8,3)):
        for rc in np.arange(max(0.2,rmax*0.05), rmax*3, max(rmax*0.08,0.3)):
            c2=chi2([vm,rc])
            if c2<best[0]: best=(c2,(vm,rc))
    res=minimize(chi2, best[1], method='Nelder-Mead', options={'maxiter':2000})
    vm_f,rc_f = max(res.x[0],0), max(res.x[1],0.05)
    
    # Filter out crazy fits (runaway parameters)
    if vm_f > 1000 or rc_f > 500:
        continue
    
    results[gal] = {
        'vm': vm_f, 'rc': rc_f,
        'Vf': g1['Vf'], 'Rd': g1['Rd'], 'L': g1['L'], 
        'MHI': g1['MHI'], 'RHI': g1['RHI'], 'T': g1['T'],
        'Reff': g1['Reff'], 'SBdisk': g1['SBdisk'], 'SBeff': g1['SBeff'],
        'D': g1['D'], 'Rmax': max(g2['R']), 'npts': n
    }

print(f"Galassie con fit buono (v_max < 1000, r_c < 500): {len(results)}")

# ============================================================
# EXPLORE ALL POSSIBLE CORRELATIONS
# ============================================================
print(f"\n{'='*70}")
print(f"CORRELAZIONI TRA PARAMETRI DEL VORTICE E PROPRIETÀ GALATTICHE")
print(f"{'='*70}")

# Gather arrays
gals = sorted(results.keys())
vm = np.array([results[g]['vm'] for g in gals])
rc = np.array([results[g]['rc'] for g in gals])
Vf = np.array([results[g]['Vf'] for g in gals])
Rd = np.array([results[g]['Rd'] for g in gals])
L = np.array([results[g]['L'] for g in gals])  # 10^9 Lsun
MHI = np.array([results[g]['MHI'] for g in gals])  # 10^9 Msun
RHI = np.array([results[g]['RHI'] for g in gals])
Reff = np.array([results[g]['Reff'] for g in gals])
T = np.array([results[g]['T'] for g in gals])
Rmax = np.array([results[g]['Rmax'] for g in gals])

# Derived quantities
MBH = np.array([10**(7.22 + 1.65*np.log10(max(v,10)/200)) for v in Vf])
Mbar = (ML_disk * L + MHI) * 1e9  # total baryonic mass in Msun (approximate)

# Correlations for v_max
print(f"\n--- CORRELAZIONI per v_max ---")
print(f"{'Proprietà':>15s}  {'r (Pearson)':>12s}  {'slope':>8s}  {'intercept':>10s}")
print(f"{'─'*50}")

props = {
    'V_flat': Vf,
    'log(V_flat)': np.log10(Vf),
    'R_disk': Rd,
    'log(L)': np.log10(np.maximum(L, 1e-4)),
    'log(M_bar)': np.log10(np.maximum(Mbar, 1e4)),
    'log(M_BH)': np.log10(MBH),
    'R_HI': RHI,
    'R_eff': Reff,
    'R_max': Rmax,
    'Type T': T,
}

for name, prop in props.items():
    mask = np.isfinite(prop) & np.isfinite(np.log10(np.maximum(vm, 0.1)))
    if np.sum(mask) < 10: continue
    
    lvm = np.log10(vm[mask])
    lprop = prop[mask] if 'log' in name else np.log10(np.maximum(prop[mask], 0.01))
    
    r = np.corrcoef(lprop, lvm)[0,1]
    c = np.polyfit(lprop, lvm, 1)
    print(f"{name:>15s}  {r:12.4f}  {c[0]:8.3f}  {c[1]:10.3f}")

print(f"\n--- CORRELAZIONI per r_c ---")
print(f"{'Proprietà':>15s}  {'r (Pearson)':>12s}  {'slope':>8s}  {'intercept':>10s}")
print(f"{'─'*50}")

for name, prop in props.items():
    mask = np.isfinite(prop) & np.isfinite(np.log10(np.maximum(rc, 0.01)))
    if np.sum(mask) < 10: continue
    
    lrc = np.log10(rc[mask])
    lprop = prop[mask] if 'log' in name else np.log10(np.maximum(prop[mask], 0.01))
    
    r_corr = np.corrcoef(lprop, lrc)[0,1]
    c = np.polyfit(lprop, lrc, 1)
    print(f"{name:>15s}  {r_corr:12.4f}  {c[0]:8.3f}  {c[1]:10.3f}")

# ============================================================
# TRY SIMPLE PHYSICAL MODEL: v_max = A * V_flat^B, r_c = C * Rd^D
# ============================================================
print(f"\n\n{'='*70}")
print(f"MODELLO PREDITTIVO: v_max e r_c da proprietà osservabili")
print(f"{'='*70}")

# Best correlations should be with V_flat (for v_max) and R_disk/R_eff (for r_c)
# Fit: log(v_max) = a * log(V_flat) + b
mask_good = (vm > 1) & (Vf > 10) & (Rd > 0.01) & (rc > 0.1)

log_vm = np.log10(vm[mask_good])
log_vf = np.log10(Vf[mask_good])
log_rc = np.log10(rc[mask_good])
log_rd = np.log10(Rd[mask_good])
log_rhi = np.log10(np.maximum(RHI[mask_good], 0.1))
log_rmax = np.log10(Rmax[mask_good])
log_reff = np.log10(np.maximum(Reff[mask_good], 0.1))

# v_max from V_flat
c_vm = np.polyfit(log_vf, log_vm, 1)
vm_pred = 10**(c_vm[0] * log_vf + c_vm[1])
resid_vm = log_vm - (c_vm[0] * log_vf + c_vm[1])
print(f"\n  v_max = {10**c_vm[1]:.4f} × V_flat^{c_vm[0]:.3f}")
print(f"  Scatter: {np.std(resid_vm):.4f} dex")
print(f"  r = {np.corrcoef(log_vf, log_vm)[0,1]:.4f}")

# r_c from R_disk
c_rc1 = np.polyfit(log_rd, log_rc, 1)
resid_rc1 = log_rc - (c_rc1[0] * log_rd + c_rc1[1])
print(f"\n  r_c = {10**c_rc1[1]:.4f} × R_disk^{c_rc1[0]:.3f}")
print(f"  Scatter: {np.std(resid_rc1):.4f} dex")
print(f"  r = {np.corrcoef(log_rd, log_rc)[0,1]:.4f}")

# r_c from R_HI  
c_rc2 = np.polyfit(log_rhi, log_rc, 1)
resid_rc2 = log_rc - (c_rc2[0] * log_rhi + c_rc2[1])
print(f"\n  r_c = {10**c_rc2[1]:.4f} × R_HI^{c_rc2[0]:.3f}")
print(f"  Scatter: {np.std(resid_rc2):.4f} dex")
print(f"  r = {np.corrcoef(log_rhi, log_rc)[0,1]:.4f}")

# r_c from Rmax (last measured radius)
c_rc3 = np.polyfit(log_rmax, log_rc, 1)
resid_rc3 = log_rc - (c_rc3[0] * log_rmax + c_rc3[1])
print(f"\n  r_c = {10**c_rc3[1]:.4f} × R_max^{c_rc3[0]:.3f}")
print(f"  Scatter: {np.std(resid_rc3):.4f} dex")
print(f"  r = {np.corrcoef(log_rmax, log_rc)[0,1]:.4f}")

# r_c from V_flat (simpler?)
c_rc4 = np.polyfit(log_vf, log_rc, 1)
resid_rc4 = log_rc - (c_rc4[0] * log_vf + c_rc4[1])
print(f"\n  r_c = {10**c_rc4[1]:.4f} × V_flat^{c_rc4[0]:.3f}")
print(f"  Scatter: {np.std(resid_rc4):.4f} dex")
print(f"  r = {np.corrcoef(log_vf, log_rc)[0,1]:.4f}")

# ============================================================
# MULTIVARIATE: r_c from V_flat + R_disk
# ============================================================
print(f"\n--- MODELLO MULTIVARIATO per r_c ---")
X = np.column_stack([log_vf, log_rd])
c_multi = np.linalg.lstsq(np.column_stack([X, np.ones(len(X))]), log_rc, rcond=None)[0]
rc_pred_multi = c_multi[0]*log_vf + c_multi[1]*log_rd + c_multi[2]
resid_multi = log_rc - rc_pred_multi
print(f"  r_c = {10**c_multi[2]:.4f} × V_flat^{c_multi[0]:.3f} × R_disk^{c_multi[1]:.3f}")
print(f"  Scatter: {np.std(resid_multi):.4f} dex")

# ============================================================
# ZERO FREE PARAMETERS MODEL: test on all galaxies
# ============================================================
print(f"\n\n{'='*70}")
print(f"TEST: MODELLO A ZERO PARAMETRI LIBERI PER GALASSIA")
print(f"{'='*70}")

# Use best single-variable formulas
# v_max from V_flat, r_c from R_disk (or V_flat)
best_vm_formula = c_vm  # log(vm) = slope * log(Vf) + intercept
# Try both r_c formulas
for rc_name, rc_formula in [('R_disk', c_rc1), ('V_flat', c_rc4), ('R_HI', c_rc2)]:
    tot_ss_cst = 0; tot_ss_bar = 0; tot_n = 0
    
    for gal in good:
        if gal not in results: continue
        g1=t1[gal]; g2=t2[gal]; n=len(g2['R'])
        
        # Predict v_max and r_c from observables ONLY
        vm_pred = 10**(best_vm_formula[0] * np.log10(g1['Vf']) + best_vm_formula[1])
        
        if rc_name == 'R_disk':
            if g1['Rd'] < 0.01: continue
            rc_pred = 10**(rc_formula[0] * np.log10(g1['Rd']) + rc_formula[1])
        elif rc_name == 'V_flat':
            rc_pred = 10**(rc_formula[0] * np.log10(g1['Vf']) + rc_formula[1])
        elif rc_name == 'R_HI':
            if g1['RHI'] < 0.01: continue
            rc_pred = 10**(rc_formula[0] * np.log10(g1['RHI']) + rc_formula[1])
        
        # CST prediction with PREDICTED parameters
        vg2=np.sign(g2['Vg'])*g2['Vg']**2
        vbar2=vg2 + ML_disk*g2['Vd']**2 + ML_bul*g2['Vb']**2
        vbar = np.sqrt(np.maximum(vbar2, 0))
        vbar_cst = vbar * np.sqrt(1+PDMF_enh)
        vv = v_LO_v(g2['R'], vm_pred, rc_pred)
        vpred = np.sqrt(vbar_cst**2 + vv**2)
        
        tot_ss_cst += np.sum((g2['Vo'] - vpred)**2)
        tot_ss_bar += np.sum((g2['Vo'] - vbar)**2)
        tot_n += n
    
    rmse_cst = np.sqrt(tot_ss_cst / tot_n)
    rmse_bar = np.sqrt(tot_ss_bar / tot_n)
    print(f"\n  r_c da {rc_name}:")
    print(f"    RMSE zero-param CST = {rmse_cst:.2f} km/s")
    print(f"    RMSE barioni soli   = {rmse_bar:.2f} km/s")
    print(f"    Miglioramento: {rmse_bar/rmse_cst:.2f}x")

