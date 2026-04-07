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

print("="*70)
print("MODELLO CST CON FUNZIONE PESO STELLARE")
print("(Analogia con l'interferenza delle binarie)")
print("="*70)
print("""
  BINARIE: Ψ = 1 + γ₀ M^η [4q/(1+q)²] exp(-a/a₀) M^β
  
  La perturbazione:
  - Si accende quando le stelle sono VICINE (exp(-a/a₀) → 1)
  - Si spegne quando sono LONTANE (exp(-a/a₀) → 0)
  - Dipende dalla MASSA relativa (4q/(1+q)²)
  
  GALASSIE: stessa logica. A ogni raggio r, le stelle locali perturbano
  il vortice. La perturbazione dipende da:
  - Quanto è densa la materia stellare a quel raggio (= v_bar(r))
  - Quanto conta rispetto al vortice (= rapporto v_bar/v_vortex)
  - Si accende nelle regioni dense (centro) e si spegne nelle regioni
    vuote (periferia o galassie nane)
  
  Modello: il raggio effettivo del core del vortice viene MODIFICATO
  localmente dalla densità stellare:
  
  r_c_eff(r) = r_c₀ × [1 + η × w(r)]
  
  dove w(r) è la funzione peso = v²_bar(r) / v²_vortex(r)
  
  - Nelle nane: v_bar piccolo ovunque → w(r) ≈ 0 → r_c_eff = r_c₀ (vortice puro)
  - Nelle massive, al centro: v_bar grande → w(r) >> 0 → r_c allargato
  - Nelle massive, al bordo: v_bar cala → w(r) → 0 → vortice puro
""")

# ============================================================
# MODELLO 1: r_c effettivo modificato localmente
# ============================================================
print(f"\n{'─'*70}")
print(f"MODELLO A: Perturbazione locale del core")
print(f"   r_c_eff(r) = C × Rd × [1 + η × (v_bar(r)/V_flat)²]")
print(f"{'─'*70}")

def model_A(params, return_details=False):
    """Local perturbation model: r_c varies with local baryonic density"""
    C, eta = params
    if C < 0.1 or eta < -1: return 1e20
    
    total_chi2 = 0; total_n = 0
    details = []
    
    for gal in good:
        g1 = t1[gal]; g2 = t2[gal]
        Vf = g1['Vf']; Rd = max(g1['Rd'], 0.01)
        n = len(g2['R'])
        rc0 = C * Rd
        
        vg2 = np.sign(g2['Vg']) * g2['Vg']**2
        vbar2 = vg2 + ML_disk*g2['Vd']**2 + ML_bul*g2['Vb']**2
        vbar = np.sqrt(np.maximum(vbar2, 0))
        vbar_cst = vbar * np.sqrt(1 + PDMF_enh)
        
        vpred = np.zeros(n)
        for i in range(n):
            # Weight: how much baryonic matter is here relative to V_flat
            w_i = (vbar[i] / max(Vf, 1))**2
            # Effective core radius at this point
            rc_eff = max(rc0 * (1 + eta * w_i), 0.05)
            vv = v_LO(g2['R'][i], Vf, rc_eff)
            vpred[i] = np.sqrt(vbar_cst[i]**2 + vv**2)
        
        ss = np.sum((g2['Vo'] - vpred)**2)
        total_chi2 += np.sum(((g2['Vo'] - vpred)/np.maximum(g2['eV'], 1))**2)
        total_n += n
        
        if return_details:
            rmse = np.sqrt(np.mean((g2['Vo'] - vpred)**2))
            rmse_bar = np.sqrt(np.mean((g2['Vo'] - vbar)**2))
            details.append({'name': gal, 'rmse': rmse, 'rmse_bar': rmse_bar, 
                          'Vf': Vf, 'T': g1['T']})
    
    if return_details:
        return total_chi2, details
    return total_chi2

# Grid search
print("\nGrid search...")
best_A = (1e30, 5, 1)
for C in np.arange(2, 15, 1):
    for eta in np.arange(-1, 5, 0.2):
        c2 = model_A([C, eta])
        if c2 < best_A[0]: best_A = (c2, C, eta)

res_A = minimize(model_A, [best_A[1], best_A[2]], method='Nelder-Mead', options={'maxiter':10000})
C_A, eta_A = res_A.x

chi2_A, details_A = model_A([C_A, eta_A], return_details=True)
rmse_arr = [d['rmse'] for d in details_A]
rmse_bar_arr = [d['rmse_bar'] for d in details_A]
tot_ss = sum(d['rmse']**2 * len(t2[d['name']]['R']) for d in details_A)
tot_ss_bar = sum(d['rmse_bar']**2 * len(t2[d['name']]['R']) for d in details_A)
tot_n = sum(len(t2[d['name']]['R']) for d in details_A)
rmse_A = np.sqrt(tot_ss / tot_n)
rmse_bar = np.sqrt(tot_ss_bar / tot_n)

print(f"\n  r_c_eff(r) = {C_A:.2f} × R_disk × [1 + {eta_A:.3f} × (v_bar(r)/V_flat)²]")
print(f"  Parametri globali: 2 (C={C_A:.2f}, η={eta_A:.3f})")
print(f"  RMSE = {rmse_A:.2f} km/s")
print(f"  Miglioramento: {rmse_bar/rmse_A:.2f}x")

# ============================================================
# MODELLO B: Interferenza completa (come le binarie)
# ============================================================
print(f"\n{'─'*70}")
print(f"MODELLO B: Interferenza stellare completa")
print(f"   v²_tot = v²_bar×(1+ε_PDMF) + v²_vortex × Ψ_gal(r)")
print(f"   Ψ_gal(r) = 1 + η × (v_bar/V_flat)^α × exp(-r/r_scale)")
print(f"{'─'*70}")

def model_B(params, return_details=False):
    """Full interference model: vortex amplitude modified by baryonic density"""
    C, eta, alpha_w, r_frac = params
    if C < 0.1 or eta < -2 or alpha_w < 0 or r_frac < 0.1: return 1e20
    
    total_chi2 = 0
    details = []
    
    for gal in good:
        g1 = t1[gal]; g2 = t2[gal]
        Vf = g1['Vf']; Rd = max(g1['Rd'], 0.01)
        n = len(g2['R'])
        rc = C * Rd
        r_scale = r_frac * Rd  # scale where perturbation activates
        
        vg2 = np.sign(g2['Vg']) * g2['Vg']**2
        vbar2 = vg2 + ML_disk*g2['Vd']**2 + ML_bul*g2['Vb']**2
        vbar = np.sqrt(np.maximum(vbar2, 0))
        vbar_cst = vbar * np.sqrt(1 + PDMF_enh)
        
        vpred = np.zeros(n)
        for i in range(n):
            vv0 = v_LO(g2['R'][i], Vf, rc)
            # Interference weight: strong where baryons are dense AND close to center
            w_bar = (vbar[i] / max(Vf, 1))**alpha_w
            w_radial = np.exp(-g2['R'][i] / max(r_scale, 0.01))
            Psi = 1 + eta * w_bar * w_radial
            vv = vv0 * np.sqrt(max(Psi, 0.01))  # modify vortex amplitude
            vpred[i] = np.sqrt(vbar_cst[i]**2 + vv**2)
        
        total_chi2 += np.sum(((g2['Vo'] - vpred)/np.maximum(g2['eV'], 1))**2)
        
        if return_details:
            rmse = np.sqrt(np.mean((g2['Vo'] - vpred)**2))
            rmse_bar = np.sqrt(np.mean((g2['Vo'] - vbar)**2))
            details.append({'name': gal, 'rmse': rmse, 'rmse_bar': rmse_bar,
                          'Vf': Vf, 'T': g1['T']})
    
    if return_details:
        return total_chi2, details
    return total_chi2

print("\nGrid search...")
best_B = (1e30, 7, 0.5, 2, 3)
for C in [4,5,6,7,8,9,10]:
    for eta in [0.2, 0.5, 1, 2, 3]:
        for alpha_w in [1, 2, 3]:
            for r_frac in [1, 2, 3, 5]:
                c2 = model_B([C, eta, alpha_w, r_frac])
                if c2 < best_B[0]: best_B = (c2, C, eta, alpha_w, r_frac)

res_B = minimize(model_B, [best_B[1], best_B[2], best_B[3], best_B[4]], 
                 method='Nelder-Mead', options={'maxiter':20000})
C_B, eta_B, alpha_B, rfrac_B = res_B.x

chi2_B, details_B = model_B([C_B, eta_B, alpha_B, rfrac_B], return_details=True)
tot_ss_B = sum(d['rmse']**2 * len(t2[d['name']]['R']) for d in details_B)
tot_n_B = sum(len(t2[d['name']]['R']) for d in details_B)
rmse_B = np.sqrt(tot_ss_B / tot_n_B)

print(f"\n  v²_tot = v²_bar(1+ε) + v²_vortex × Ψ(r)")
print(f"  Ψ(r) = 1 + {eta_B:.3f} × (v_bar/V_flat)^{alpha_B:.2f} × exp(-r/{rfrac_B:.2f}×Rd)")
print(f"  r_c = {C_B:.2f} × R_disk")
print(f"  Parametri globali: 4 (C, η, α_w, r_frac)")
print(f"  RMSE = {rmse_B:.2f} km/s")
print(f"  Miglioramento: {rmse_bar/rmse_B:.2f}x")

# ============================================================
# MODELLO C: Peso binario-like con soglia di attivazione
# ============================================================
print(f"\n{'─'*70}")
print(f"MODELLO C: Funzione peso con soglia (come w(M) per le binarie)")
print(f"   Il peso si 'accende' quando la densità barionica supera una soglia")
print(f"   e si 'spegne' sotto — esattamente come w(M) = exp(-|M/M☉-1|)")
print(f"{'─'*70}")

def model_C(params, return_details=False):
    """Weight function model with activation threshold"""
    C, eta, v_threshold = params
    if C < 0.1 or eta < -2 or v_threshold < 1: return 1e20
    
    total_chi2 = 0
    details = []
    
    for gal in good:
        g1 = t1[gal]; g2 = t2[gal]
        Vf = g1['Vf']; Rd = max(g1['Rd'], 0.01)
        n = len(g2['R'])
        rc0 = C * Rd
        
        vg2 = np.sign(g2['Vg']) * g2['Vg']**2
        vbar2 = vg2 + ML_disk*g2['Vd']**2 + ML_bul*g2['Vb']**2
        vbar = np.sqrt(np.maximum(vbar2, 0))
        vbar_cst = vbar * np.sqrt(1 + PDMF_enh)
        
        vpred = np.zeros(n)
        for i in range(n):
            # Weight function analogous to w(M) = exp(-|M/M☉ - 1|)
            # Here: w = exp(-|v_bar/v_threshold - 1|)
            # When v_bar ≈ v_threshold: w ≈ 1 (perturbation ON)
            # When v_bar << v_threshold: w → exp(-1 + v_bar/v_th) → small (OFF in dwarfs)
            # When v_bar >> v_threshold: w → small (OFF in very dense regions)
            # Actually simpler: w = 1 - exp(-(v_bar/v_threshold)²)
            # This is 0 when v_bar << v_threshold, and ~1 when v_bar >> v_threshold
            w_i = 1 - np.exp(-(vbar[i]/v_threshold)**2)
            
            rc_eff = max(rc0 * (1 + eta * w_i), 0.05)
            vv = v_LO(g2['R'][i], Vf, rc_eff)
            vpred[i] = np.sqrt(vbar_cst[i]**2 + vv**2)
        
        total_chi2 += np.sum(((g2['Vo'] - vpred)/np.maximum(g2['eV'], 1))**2)
        
        if return_details:
            rmse = np.sqrt(np.mean((g2['Vo'] - vpred)**2))
            rmse_bar = np.sqrt(np.mean((g2['Vo'] - vbar)**2))
            details.append({'name': gal, 'rmse': rmse, 'rmse_bar': rmse_bar,
                          'Vf': Vf, 'T': g1['T']})
    
    if return_details:
        return total_chi2, details
    return total_chi2

print("\nGrid search...")
best_C = (1e30, 7, 1, 50)
for C in [3,4,5,6,7,8,9,10,12]:
    for eta in [0.5, 1, 2, 3, 5, 8]:
        for vth in [20, 40, 60, 80, 100, 130, 160, 200]:
            c2 = model_C([C, eta, vth])
            if c2 < best_C[0]: best_C = (c2, C, eta, vth)

res_C = minimize(model_C, [best_C[1], best_C[2], best_C[3]], 
                 method='Nelder-Mead', options={'maxiter':10000})
C_C, eta_C, vth_C = res_C.x

chi2_C, details_C = model_C([C_C, eta_C, vth_C], return_details=True)
tot_ss_C = sum(d['rmse']**2 * len(t2[d['name']]['R']) for d in details_C)
tot_n_C = sum(len(t2[d['name']]['R']) for d in details_C)
rmse_C = np.sqrt(tot_ss_C / tot_n_C)

print(f"\n  r_c_eff(r) = {C_C:.2f} × Rd × [1 + {eta_C:.2f} × w(r)]")
print(f"  w(r) = 1 - exp(-(v_bar(r)/{vth_C:.0f})²)")
print(f"  Soglia v_threshold = {vth_C:.0f} km/s")
print(f"  Parametri globali: 3 (C, η, v_threshold)")
print(f"  RMSE = {rmse_C:.2f} km/s")
print(f"  Miglioramento: {rmse_bar/rmse_C:.2f}x")

print(f"""
  INTERPRETAZIONE FISICA:
  
  La soglia v_threshold = {vth_C:.0f} km/s è la velocità barionica
  sopra la quale le stelle perturbano significativamente il vortice.
  
  - Nane (v_bar ~ 20-50 km/s): w ≈ {1-np.exp(-(30/vth_C)**2):.3f} → perturbazione quasi OFF
  - Spirali medie (v_bar ~ 100 km/s): w ≈ {1-np.exp(-(100/vth_C)**2):.3f} → perturbazione parziale  
  - Spirali massive (v_bar ~ 200 km/s): w ≈ {1-np.exp(-(200/vth_C)**2):.3f} → perturbazione piena
""")

# ============================================================
# ANALISI PER TIPO MORFOLOGICO
# ============================================================
print(f"\n{'='*70}")
print(f"PRESTAZIONI PER TIPO MORFOLOGICO")
print(f"{'='*70}")

type_names = {0:'S0',1:'Sa',2:'Sab',3:'Sb',4:'Sbc',5:'Sc',6:'Scd',7:'Sd',8:'Sdm',9:'Sm',10:'Im',11:'BCD'}

# Use best model (C)
ts = {}
for d in details_C:
    T = d['T']
    if T not in ts: ts[T] = {'rmse': [], 'rmse_bar': [], 'imp': []}
    ts[T]['rmse'].append(d['rmse'])
    ts[T]['rmse_bar'].append(d['rmse_bar'])
    ts[T]['imp'].append(d['rmse_bar'] / max(d['rmse'], 0.1))

print(f"\n  {'Tipo':>4s} {'Nome':>4s} {'N':>3s} │ {'RMSE_CST':>8s} {'RMSE_bar':>8s} │ {'Migl':>6s}")
print(f"  {'─'*48}")
for T in sorted(ts.keys()):
    s = ts[T]
    print(f"  {T:4d} {type_names.get(T,'?'):>4s} {len(s['rmse']):3d} │ {np.median(s['rmse']):8.1f} {np.median(s['rmse_bar']):8.1f} │ {np.median(s['imp']):5.1f}x")

# ============================================================
# TABELLA COMPARATIVA FINALE
# ============================================================
print(f"\n\n{'='*70}")
print(f"TABELLA COMPARATIVA FINALE — TUTTI I MODELLI")
print(f"{'='*70}")
print(f"""
  Modello                                   Par   RMSE    Migl.  Note
  ──────────────────────────────────────────────────────────────────────
  Solo barioni                               0    61.9    1.0x   baseline
  MOND (a₀)                                  1    ~25     ~2.5x  1 param fenomenologico
  CST semplice: rc=C×Rd                      1    25.6    2.4x   1 param universale
  CST locale-A: rc_eff = C×Rd×[1+η×(vb/Vf)²] 2   {rmse_A:.1f}    {rmse_bar/rmse_A:.1f}x   perturbazione proporzionale
  CST soglia-C: rc_eff + w(v_bar)            3    {rmse_C:.1f}    {rmse_bar/rmse_C:.1f}x   con soglia di attivazione
  CST interf-B: Ψ(r) completa               4    {rmse_B:.1f}    {rmse_bar/rmse_B:.1f}x   interferenza completa
  CST fit individuale                      258    18.4    3.4x   benchmark (overfitting)
  
  
  FORMULA FISICA COMPLETA DEL MODELLO CST GALATTICO:
  ═══════════════════════════════════════════════════
  
  v²_tot(r) = v²_bar(r) × (1 + ε_PDMF) + v²_vortice(r)
  
  dove:
    v_vortice(r) = V_flat × LO(r / r_c_eff(r))
    r_c_eff(r) = C × R_disk × [1 + η × w(r)]
    w(r) = 1 - exp(-(v_bar(r) / v_th)²)
    
    ε_PDMF = 0.215  (calcolato dalla popolazione stellare)
    C = {C_C:.2f}        (dimensione base del vortice in unità di R_disk)
    η = {eta_C:.2f}        (intensità della perturbazione stellare)
    v_th = {vth_C:.0f} km/s  (soglia di attivazione della perturbazione)
    
  SIGNIFICATO FISICO:
  - Il vortice SMBH ha un core di raggio ~ {C_C:.0f}×R_disk
  - Dove la materia barionica è densa (v_bar > {vth_C:.0f} km/s), le stelle
    perturbano il vortice e ne allargano il core di un fattore 1+η
  - Dove la materia è rarefatta (nane, periferia), il vortice è "puro"
  
  Questa è la STESSA struttura dell'interferenza binaria:
    Binarie: Ψ = 1 + γ₀ × (fattore massa) × exp(-a/a₀)
    Galassie: w = 1 - exp(-(v_bar/v_th)²)
  
  In entrambi i casi: una funzione che si accende quando la materia
  è vicina e intensa, e si spegne quando è lontana o rarefatta.
""")

