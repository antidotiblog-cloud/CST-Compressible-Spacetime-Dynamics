import numpy as np
from scipy.optimize import minimize, differential_evolution

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
                      'Q':int(p[17]),'L':float(p[7]),'MHI':float(p[13]),
                      'SBdisk':float(p[12]),'SBeff':float(p[10]),'Reff':float(p[9])}
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

# Precompute baryonic arrays for speed
precomp = {}
for gal in good:
    g2 = t2[gal]
    vg2 = np.sign(g2['Vg'])*g2['Vg']**2
    vbar2 = vg2 + ML_disk*g2['Vd']**2 + ML_bul*g2['Vb']**2
    vbar = np.sqrt(np.maximum(vbar2, 0))
    vbar_cst = vbar * np.sqrt(1+PDMF_enh)
    # Cumulative baryonic mass proxy: sum of v²_bar(r)*r up to each point
    # M_bar(<r) ∝ v²_bar(r) × r / G
    M_cum = np.cumsum(vbar**2 * g2['R'])
    M_cum_norm = M_cum / max(M_cum[-1], 1)
    precomp[gal] = {'vbar': vbar, 'vbar_cst': vbar_cst, 'vbar2': vbar2,
                    'M_cum_norm': M_cum_norm, 'vbul': g2['Vb'], 'vdisk': g2['Vd']}

print("="*70)
print("RAFFINAMENTO DEL MODELLO CST GALATTICO")
print("Obiettivo: chiudere il gap 22 → 18 km/s")
print("="*70)

# ============================================================
# ANALISI DIAGNOSTICA: dove perde il modello B?
# ============================================================
print(f"\n{'─'*70}")
print(f"DIAGNOSTICA: dove il modello B (RMSE=22.2) sbaglia?")
print(f"{'─'*70}")

# Model B params from previous run
C_B = 10.31; eta_B = 7.05; rfrac_B = 1.67

per_gal_B = []
for gal in good:
    g1=t1[gal]; g2=t2[gal]; pc=precomp[gal]
    Vf=g1['Vf']; Rd=max(g1['Rd'],0.01)
    rc=C_B*Rd; rscale=rfrac_B*Rd
    
    vpred = np.zeros(len(g2['R']))
    for i in range(len(g2['R'])):
        vv0 = v_LO(g2['R'][i], Vf, rc)
        Psi = 1 + eta_B * np.exp(-g2['R'][i]/max(rscale,0.01))
        vv = vv0 * np.sqrt(max(Psi, 0.01))
        vpred[i] = np.sqrt(pc['vbar_cst'][i]**2 + vv**2)
    
    rmse = np.sqrt(np.mean((g2['Vo']-vpred)**2))
    rmse_bar = np.sqrt(np.mean((g2['Vo']-pc['vbar'][i])**2))
    
    # Where does it fail? inner vs outer
    n = len(g2['R'])
    mid = n//2
    rmse_inner = np.sqrt(np.mean((g2['Vo'][:mid]-vpred[:mid])**2)) if mid > 1 else 0
    rmse_outer = np.sqrt(np.mean((g2['Vo'][mid:]-vpred[mid:])**2)) if mid > 1 else 0
    
    # Average residual sign (overpredict or underpredict?)
    mean_resid = np.mean(g2['Vo'] - vpred)
    
    per_gal_B.append({'name':gal, 'rmse':rmse, 'Vf':Vf, 'Rd':Rd, 'T':g1['T'],
                      'rmse_in':rmse_inner, 'rmse_out':rmse_outer,
                      'mean_resid':mean_resid, 'n':n,
                      'has_bulge': max(g2['Vb']) > 10})

# Sort by error
per_gal_B.sort(key=lambda x: -x['rmse'])

print(f"\n  Le 15 galassie più problematiche:")
print(f"  {'Galassia':12s} {'Vf':>5s} {'Rd':>5s} {'T':>2s} {'Bulge':>5s} │ {'RMSE':>5s} {'Inner':>5s} {'Outer':>5s} {'<resid>':>7s}")
for g in per_gal_B[:15]:
    tn = {0:'S0',1:'Sa',2:'Sab',3:'Sb',4:'Sbc',5:'Sc',6:'Scd',7:'Sd',8:'Sdm',9:'Sm',10:'Im'}
    print(f"  {g['name']:12s} {g['Vf']:5.0f} {g['Rd']:5.2f} {tn.get(g['T'],'?'):>2s} {'Y' if g['has_bulge'] else 'N':>5s} │ {g['rmse']:5.1f} {g['rmse_in']:5.1f} {g['rmse_out']:5.1f} {g['mean_resid']:+7.1f}")

# Count patterns
n_inner_worse = sum(1 for g in per_gal_B if g['rmse_in'] > g['rmse_out'] * 1.5)
n_outer_worse = sum(1 for g in per_gal_B if g['rmse_out'] > g['rmse_in'] * 1.5)
n_overpredict = sum(1 for g in per_gal_B if g['mean_resid'] < -5)
n_underpredict = sum(1 for g in per_gal_B if g['mean_resid'] > 5)
n_bulge = sum(1 for g in per_gal_B[:20] if g['has_bulge'])

print(f"\n  Pattern diagnostico:")
print(f"    Errore peggiore al centro: {n_inner_worse} galassie")
print(f"    Errore peggiore al bordo:  {n_outer_worse} galassie")
print(f"    Sovra-predice (v_CST > v_obs): {n_overpredict} galassie")
print(f"    Sotto-predice (v_CST < v_obs): {n_underpredict} galassie")
print(f"    Con bulge tra le top-20 errori: {n_bulge}/20")

# ============================================================
# MODELLO D: Vortice + perturbazione INTEGRATA (massa cumulata)
# ============================================================
print(f"\n\n{'='*70}")
print(f"MODELLO D: Perturbazione dalla massa barionica CUMULATA")
print(f"{'='*70}")
print("""
  Idea di Michele: non è la densità LOCALE che conta, è la MASSA
  TOTALE delle stelle fra il BH e il raggio r.
  
  Come nelle binarie: l'effetto dipende dalla massa relativa q = M2/M1.
  Qui: q(r) = M_bar(<r) / M_BH
  
  Quando q grande (tante stelle fra te e il BH) → perturbazione forte
  Quando q piccolo (poche stelle, nane) → perturbazione debole
""")

def model_D(params, return_details=False):
    """Cumulative mass perturbation model"""
    C, eta, q_scale = params
    if C < 0.5 or eta < -2 or q_scale < 0: return 1e20
    
    total_chi2 = 0
    details = []
    
    for gal in good:
        g1=t1[gal]; g2=t2[gal]; pc=precomp[gal]
        Vf=g1['Vf']; Rd=max(g1['Rd'],0.01)
        MBH = 10**(7.22+1.65*np.log10(max(Vf,10)/200))
        rc0 = C * Rd
        n = len(g2['R'])
        
        # Cumulative baryonic mass at each radius (proxy from velocity)
        # M_bar(<r) ~ v²_bar(r) × r / G
        G_N = 6.674e-11; Msun = 1.989e30; kpc_m = 3.086e19
        
        vpred = np.zeros(n)
        for i in range(n):
            # Mass ratio: baryonic mass inside r vs SMBH mass
            r_m = g2['R'][i] * kpc_m
            vbar_ms = pc['vbar'][i] * 1e3  # km/s to m/s
            M_bar_r = vbar_ms**2 * r_m / G_N / Msun  # in solar masses
            q = M_bar_r / max(MBH, 1)
            
            # Weight: analogous to binary interference
            # w(q) = 1 - exp(-q/q_scale) 
            # q small (dwarfs, few stars near BH) → w ≈ 0 → pure vortex
            # q large (massive, many stars) → w → 1 → full perturbation
            w_q = 1 - np.exp(-q / max(q_scale, 0.01))
            
            rc_eff = max(rc0 * (1 + eta * w_q), 0.05)
            vv = v_LO(g2['R'][i], Vf, rc_eff)
            vpred[i] = np.sqrt(pc['vbar_cst'][i]**2 + vv**2)
        
        total_chi2 += np.sum(((g2['Vo']-vpred)/np.maximum(g2['eV'],1))**2)
        
        if return_details:
            rmse = np.sqrt(np.mean((g2['Vo']-vpred)**2))
            rmse_bar = np.sqrt(np.mean((g2['Vo']-pc['vbar'])**2))
            details.append({'name':gal,'rmse':rmse,'rmse_bar':rmse_bar,'Vf':Vf,'T':g1['T']})
    
    if return_details: return total_chi2, details
    return total_chi2

print("Grid search modello D...")
best_D = (1e30, 5, 1, 100)
for C in [3,4,5,6,7,8,9,10]:
    for eta in [0.3, 0.5, 1, 2, 3, 5]:
        for qs in [10, 50, 100, 500, 1000, 5000]:
            c2 = model_D([C, eta, qs])
            if c2 < best_D[0]: best_D = (c2, C, eta, qs)

res_D = minimize(model_D, [best_D[1], best_D[2], best_D[3]], method='Nelder-Mead', options={'maxiter':10000})
C_D, eta_D, qs_D = res_D.x

chi2_D, details_D = model_D([C_D, eta_D, qs_D], return_details=True)
tot_ss_D = sum(d['rmse']**2 * len(t2[d['name']]['R']) for d in details_D)
tot_ss_bar = sum(d['rmse_bar']**2 * len(t2[d['name']]['R']) for d in details_D)
tot_n = sum(len(t2[d['name']]['R']) for d in details_D)
rmse_D = np.sqrt(tot_ss_D/tot_n)
rmse_bar = np.sqrt(tot_ss_bar/tot_n)

print(f"\n  r_c(r) = {C_D:.2f} × Rd × [1 + {eta_D:.2f} × w(q(r))]")
print(f"  w(q) = 1 - exp(-q/{qs_D:.0f})")
print(f"  q(r) = M_bar(<r) / M_BH")
print(f"  Parametri globali: 3 (C, η, q_scale)")
print(f"  RMSE = {rmse_D:.2f} km/s")
print(f"  Miglioramento: {rmse_bar/rmse_D:.2f}x")

# ============================================================
# MODELLO E: Combinazione radiale + massa cumulata
# ============================================================
print(f"\n\n{'='*70}")
print(f"MODELLO E: Doppia perturbazione (radiale + massa cumulata)")
print(f"{'='*70}")

def model_E(params, return_details=False):
    """Combined radial + cumulative mass perturbation"""
    C, eta_rad, r_frac, eta_mass, q_scale = params
    if C < 0.5 or r_frac < 0.1 or q_scale < 1: return 1e20
    
    total_chi2 = 0
    details = []
    G_N = 6.674e-11; Msun = 1.989e30; kpc_m = 3.086e19
    
    for gal in good:
        g1=t1[gal]; g2=t2[gal]; pc=precomp[gal]
        Vf=g1['Vf']; Rd=max(g1['Rd'],0.01)
        MBH = 10**(7.22+1.65*np.log10(max(Vf,10)/200))
        rc0 = C * Rd
        rscale = r_frac * Rd
        n = len(g2['R'])
        
        vpred = np.zeros(n)
        for i in range(n):
            # Radial factor: vortex amplified near center
            f_rad = 1 + eta_rad * np.exp(-g2['R'][i]/max(rscale,0.01))
            
            # Mass ratio factor: perturbation from accumulated stellar mass
            r_m = g2['R'][i] * kpc_m
            vbar_ms = pc['vbar'][i] * 1e3
            M_bar_r = vbar_ms**2 * r_m / G_N / Msun
            q = M_bar_r / max(MBH, 1)
            w_q = 1 - np.exp(-q / max(q_scale, 0.01))
            f_mass = 1 + eta_mass * w_q
            
            # Combined: vortex modified by both factors
            vv0 = v_LO(g2['R'][i], Vf, rc0 * f_mass)
            vv = vv0 * np.sqrt(max(f_rad, 0.01))
            vpred[i] = np.sqrt(pc['vbar_cst'][i]**2 + vv**2)
        
        total_chi2 += np.sum(((g2['Vo']-vpred)/np.maximum(g2['eV'],1))**2)
        
        if return_details:
            rmse = np.sqrt(np.mean((g2['Vo']-vpred)**2))
            rmse_bar = np.sqrt(np.mean((g2['Vo']-pc['vbar'])**2))
            details.append({'name':gal,'rmse':rmse,'rmse_bar':rmse_bar,'Vf':Vf,'T':g1['T']})
    
    if return_details: return total_chi2, details
    return total_chi2

print("Grid search modello E (5 parametri)...")
best_E = (1e30, [8, 5, 2, 1, 100])
for C in [5,7,9,11]:
    for er in [2,5,8]:
        for rf in [1,2,3]:
            for em in [0.5,1,2]:
                for qs in [50,200,1000]:
                    c2 = model_E([C,er,rf,em,qs])
                    if c2 < best_E[0]: best_E = (c2, [C,er,rf,em,qs])

res_E = minimize(model_E, best_E[1], method='Nelder-Mead', options={'maxiter':20000})
C_E, er_E, rf_E, em_E, qs_E = res_E.x

chi2_E, details_E = model_E(res_E.x, return_details=True)
tot_ss_E = sum(d['rmse']**2*len(t2[d['name']]['R']) for d in details_E)
tot_n_E = sum(len(t2[d['name']]['R']) for d in details_E)
rmse_E = np.sqrt(tot_ss_E/tot_n_E)

print(f"\n  Vortice amplificato: v_vort × √[1 + {er_E:.2f} × exp(-r/{rf_E:.2f}×Rd)]")
print(f"  Core perturbato: r_c = {C_E:.2f}×Rd × [1 + {em_E:.2f} × w(q)]")
print(f"  w(q) = 1 - exp(-q/{qs_E:.0f}),  q = M_bar(<r)/M_BH")
print(f"  Parametri globali: 5")
print(f"  RMSE = {rmse_E:.2f} km/s")
print(f"  Miglioramento: {rmse_bar/rmse_E:.2f}x")

# ============================================================
# MODELLO F: Il "Modello di Michele" — interferenza galattica completa
# ============================================================
print(f"\n\n{'='*70}")
print(f"MODELLO F: INTERFERENZA GALATTICA COMPLETA")
print(f"  (Il modello di Michele: come binarie, ma per il vortice)")
print(f"{'='*70}")
print("""
  Nelle binarie CST: Ψ = 1 + γ₀ M^η [4q/(1+q)²] exp(-a/a₀) M^β
  
  Per le galassie, il vortice è perturbato dalle stelle con la STESSA
  struttura matematica:
  
  v²_tot(r) = v²_bar(r)(1+ε) + v²_vortice(r) × Ψ_gal(r)
  
  Ψ_gal(r) = 1 + γ_gal × [q(r)/(1+q(r))²]^p × exp(-r/r₀) × (M_BH/M_ref)^β_gal
  
  dove:
  - q(r) = M_bar(<r)/M_BH = rapporto massa (come nelle binarie)
  - exp(-r/r₀) = attenuazione con la distanza (come exp(-a/a₀))
  - (M_BH/M_ref)^β_gal = dipendenza dalla massa del BH
  - Il fattore 4q/(1+q)² è massimo quando q=1 (come nelle binarie!)
""")

def model_F(params, return_details=False):
    """Full galactic interference model — binary-like"""
    C, gamma_g, p_q, r0_frac, beta_g = params
    if C < 0.5 or gamma_g < 0 or p_q < 0 or r0_frac < 0.1: return 1e20
    
    total_chi2 = 0
    details = []
    G_N = 6.674e-11; Msun = 1.989e30; kpc_m = 3.086e19
    M_ref = 1e7  # reference SMBH mass
    
    for gal in good:
        g1=t1[gal]; g2=t2[gal]; pc=precomp[gal]
        Vf=g1['Vf']; Rd=max(g1['Rd'],0.01)
        MBH = 10**(7.22+1.65*np.log10(max(Vf,10)/200))
        rc = C * Rd
        r0 = r0_frac * Rd
        n = len(g2['R'])
        
        mass_factor = (MBH / M_ref)**beta_g
        
        vpred = np.zeros(n)
        for i in range(n):
            vv0 = v_LO(g2['R'][i], Vf, rc)
            
            # Mass ratio at this radius
            r_m = g2['R'][i] * kpc_m
            vbar_ms = pc['vbar'][i] * 1e3
            M_bar_r = max(vbar_ms**2 * r_m / G_N / Msun, 1)
            q = M_bar_r / max(MBH, 1)
            
            # Binary-like interference factor
            # 4q/(1+q)² is max at q=1, → 0 for q→0 and q→∞
            binary_factor = (4*q / (1+q)**2)**p_q
            radial_factor = np.exp(-g2['R'][i] / max(r0, 0.01))
            
            Psi = 1 + gamma_g * binary_factor * radial_factor * mass_factor
            
            vv = vv0 * np.sqrt(max(Psi, 0.01))
            vpred[i] = np.sqrt(pc['vbar_cst'][i]**2 + vv**2)
        
        total_chi2 += np.sum(((g2['Vo']-vpred)/np.maximum(g2['eV'],1))**2)
        
        if return_details:
            rmse = np.sqrt(np.mean((g2['Vo']-vpred)**2))
            rmse_bar = np.sqrt(np.mean((g2['Vo']-pc['vbar'])**2))
            details.append({'name':gal,'rmse':rmse,'rmse_bar':rmse_bar,
                          'Vf':Vf,'T':g1['T'],'MBH':MBH})
    
    if return_details: return total_chi2, details
    return total_chi2

print("Grid search modello F (5 parametri)...")
best_F = (1e30, [8, 3, 1, 2, 0.3])
for C in [5,7,8,10,12]:
    for gg in [1,3,5,8,12]:
        for pq in [0.3, 0.5, 1, 2]:
            for rf in [1,2,3,5]:
                for bg in [-0.5, 0, 0.3, 0.67]:
                    c2 = model_F([C,gg,pq,rf,bg])
                    if c2 < best_F[0]: best_F = (c2, [C,gg,pq,rf,bg])

res_F = minimize(model_F, best_F[1], method='Nelder-Mead', options={'maxiter':30000})

chi2_F, details_F = model_F(res_F.x, return_details=True)
tot_ss_F = sum(d['rmse']**2*len(t2[d['name']]['R']) for d in details_F)
tot_n_F = sum(len(t2[d['name']]['R']) for d in details_F)
rmse_F = np.sqrt(tot_ss_F/tot_n_F)

C_F, gg_F, pq_F, rf_F, bg_F = res_F.x

print(f"\n  MODELLO F — Interferenza galattica completa:")
print(f"  r_c = {C_F:.2f} × R_disk")
print(f"  Ψ(r) = 1 + {gg_F:.2f} × [4q/(1+q)²]^{pq_F:.2f} × exp(-r/{rf_F:.2f}×Rd) × (M_BH/10⁷)^{bg_F:.3f}")
print(f"  q(r) = M_bar(<r) / M_BH")
print(f"  Parametri globali: 5")
print(f"  RMSE = {rmse_F:.2f} km/s")
print(f"  Miglioramento: {rmse_bar/rmse_F:.2f}x")

# Show β_gal
print(f"\n  β_gal = {bg_F:.3f}")
if abs(bg_F - 0.667) < 0.2:
    print(f"  → COMPATIBILE con 2/3! L'esponente universale appare ancora!")
elif abs(bg_F) < 0.15:
    print(f"  → VICINO A ZERO: l'effetto non dipende dalla massa del BH")
    print(f"    (coerente con β_net ≈ 0 per BH: massa e gravità superficiale si cancellano)")

# Per-type performance
print(f"\n  Prestazioni per tipo:")
tn = {0:'S0',1:'Sa',2:'Sab',3:'Sb',4:'Sbc',5:'Sc',6:'Scd',7:'Sd',8:'Sdm',9:'Sm',10:'Im',11:'BCD'}
ts = {}
for d in details_F:
    T = d['T']
    if T not in ts: ts[T] = {'r':[],'rb':[]}
    ts[T]['r'].append(d['rmse']); ts[T]['rb'].append(d['rmse_bar'])

print(f"  {'Tipo':>4s} {'N':>3s} {'RMSE':>6s} {'RMSE_bar':>8s} {'Migl':>6s}")
for T in sorted(ts.keys()):
    s=ts[T]; n=len(s['r'])
    print(f"  {tn.get(T,'?'):>4s} {n:3d} {np.median(s['r']):6.1f} {np.median(s['rb']):8.1f} {np.median(np.array(s['rb'])/np.maximum(np.array(s['r']),0.1)):5.1f}x")

# ============================================================
# TABELLA FINALE
# ============================================================
print(f"\n\n{'='*70}")
print(f"TABELLA COMPARATIVA COMPLETA")
print(f"{'='*70}")
print(f"""
  Modello                                Par  RMSE    Migl.  
  ────────────────────────────────────────────────────────────
  Solo barioni                            0   61.9    1.0x   
  MOND (a₀)                               1   ~25     ~2.5x  
  CST semplice (r_c=C×Rd)                 1   25.6    2.4x   
  CST soglia-C (r_c + w(v_bar))           3   25.2    2.5x   
  CST massa-D (r_c + w(q))                3   {rmse_D:.1f}    {rmse_bar/rmse_D:.1f}x   
  CST doppia-E (radiale + massa)           5   {rmse_E:.1f}    {rmse_bar/rmse_E:.1f}x   
  CST interf-F (modello binario-like)      5   {rmse_F:.1f}    {rmse_bar/rmse_F:.1f}x   
  CST fit individuale                    258   18.4    3.4x   
""")

