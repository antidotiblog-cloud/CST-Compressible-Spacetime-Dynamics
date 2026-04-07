import numpy as np
from scipy.optimize import minimize

PDMF_enh = 0.215; ML_disk = 0.5; ML_bul = 0.7
G_N = 6.674e-11; Msun = 1.989e30; kpc_m = 3.086e19

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

# Precompute
precomp = {}
for gal in good:
    g2 = t2[gal]
    vg2 = np.sign(g2['Vg'])*g2['Vg']**2
    vbar2 = vg2 + ML_disk*g2['Vd']**2 + ML_bul*g2['Vb']**2
    vbar = np.sqrt(np.maximum(vbar2, 0))
    vbar_cst = vbar * np.sqrt(1+PDMF_enh)
    precomp[gal] = {'vbar': vbar, 'vbar_cst': vbar_cst}

# ============================================================
# STEP 1: Diagnostica — dove perde il modello F?
# ============================================================
print("="*70)
print("STEP 1: DIAGNOSTICA DETTAGLIATA DEL MODELLO F")
print("="*70)

# Model F params
C_F=19.69; gg_F=110.18; pq_F=0.24; rf_F=2.47; bg_F=0.258; M_ref=1e7

per_gal = []
all_resid_r = []  # (r/Rd, residual, gal_type)

for gal in good:
    g1=t1[gal]; g2=t2[gal]; pc=precomp[gal]
    Vf=g1['Vf']; Rd=max(g1['Rd'],0.01)
    MBH = 10**(7.22+1.65*np.log10(max(Vf,10)/200))
    rc=C_F*Rd; r0=rf_F*Rd
    mass_factor = (MBH/M_ref)**bg_F
    n=len(g2['R'])
    
    vpred = np.zeros(n)
    for i in range(n):
        vv0 = v_LO(g2['R'][i], Vf, rc)
        r_m = g2['R'][i]*kpc_m
        vbar_ms = pc['vbar'][i]*1e3
        M_bar_r = max(vbar_ms**2*r_m/G_N/Msun, 1)
        q = M_bar_r/max(MBH,1)
        binary_factor = (4*q/(1+q)**2)**pq_F
        radial_factor = np.exp(-g2['R'][i]/max(r0,0.01))
        Psi = 1 + gg_F*binary_factor*radial_factor*mass_factor
        vv = vv0*np.sqrt(max(Psi,0.01))
        vpred[i] = np.sqrt(pc['vbar_cst'][i]**2 + vv**2)
        
        # Store residual vs normalized radius
        all_resid_r.append((g2['R'][i]/Rd, g2['Vo'][i]-vpred[i], g1['T'], Vf))
    
    rmse = np.sqrt(np.mean((g2['Vo']-vpred)**2))
    rmse_bar = np.sqrt(np.mean((g2['Vo']-pc['vbar'])**2))
    per_gal.append({'name':gal,'rmse':rmse,'rmse_bar':rmse_bar,'Vf':Vf,'T':g1['T'],
                    'Rd':Rd,'MBH':MBH,'n':n})

# Residual pattern vs r/Rd
all_resid_r = np.array(all_resid_r, dtype=[('rRd','f8'),('resid','f8'),('T','i4'),('Vf','f8')])

print(f"\nResidui medi per intervallo di r/Rd:")
print(f"  {'r/Rd':>10s}  {'N':>5s}  {'<resid>':>8s}  {'σ':>6s}  {'Bias':>10s}")
for lo, hi in [(0,1),(1,2),(2,3),(3,5),(5,8),(8,15),(15,50)]:
    mask = (all_resid_r['rRd']>=lo) & (all_resid_r['rRd']<hi)
    if np.sum(mask) > 5:
        r = all_resid_r['resid'][mask]
        bias = "OVER" if np.mean(r) < -3 else "UNDER" if np.mean(r) > 3 else "OK"
        print(f"  [{lo:2d},{hi:2d})     {np.sum(mask):5d}  {np.mean(r):+8.2f}  {np.std(r):6.2f}  {bias:>10s}")

# By galaxy mass
print(f"\nResidui medi per range di V_flat:")
for vlo,vhi in [(0,80),(80,130),(130,200),(200,350)]:
    mask = (all_resid_r['Vf']>=vlo) & (all_resid_r['Vf']<vhi)
    if np.sum(mask)>5:
        r = all_resid_r['resid'][mask]
        print(f"  Vf [{vlo:3d},{vhi:3d}): N={np.sum(mask):4d}, <resid>={np.mean(r):+6.2f}, σ={np.std(r):5.2f}")

# ============================================================
# STEP 2: Provare β_gal = 2/3 (fissato, fisicamente motivato)
# ============================================================
print(f"\n\n{'='*70}")
print(f"STEP 2: β_gal FISSATO A 2/3 (predizione teorica)")
print(f"{'='*70}")

def model_F_fixed_beta(params):
    C, gamma_g, p_q, r0_frac = params
    beta_g = 2.0/3.0  # FISSO!
    if C<0.5 or gamma_g<0 or p_q<0 or r0_frac<0.1: return 1e20
    total=0
    for gal in good:
        g1=t1[gal]; g2=t2[gal]; pc=precomp[gal]
        Vf=g1['Vf']; Rd=max(g1['Rd'],0.01)
        MBH=10**(7.22+1.65*np.log10(max(Vf,10)/200))
        rc=C*Rd; r0=r0_frac*Rd; mf=(MBH/M_ref)**beta_g
        n=len(g2['R'])
        vpred=np.zeros(n)
        for i in range(n):
            vv0=v_LO(g2['R'][i],Vf,rc)
            r_m=g2['R'][i]*kpc_m; vbar_ms=pc['vbar'][i]*1e3
            M_bar_r=max(vbar_ms**2*r_m/G_N/Msun,1)
            q=M_bar_r/max(MBH,1)
            Psi=1+gamma_g*(4*q/(1+q)**2)**p_q*np.exp(-g2['R'][i]/max(r0,0.01))*mf
            vv=vv0*np.sqrt(max(Psi,0.01))
            vpred[i]=np.sqrt(pc['vbar_cst'][i]**2+vv**2)
        total+=np.sum(((g2['Vo']-vpred)/np.maximum(g2['eV'],1))**2)
    return total

best = (1e30,[15,80,0.3,2])
for C in [10,15,20,25]:
    for gg in [30,60,100,200]:
        for pq in [0.2,0.5,1]:
            for rf in [1,2,3,5]:
                c2=model_F_fixed_beta([C,gg,pq,rf])
                if c2<best[0]: best=(c2,[C,gg,pq,rf])
res=minimize(model_F_fixed_beta, best[1], method='Nelder-Mead', options={'maxiter':20000})
C2,gg2,pq2,rf2 = res.x

# RMSE
tot_ss=0; tot_n=0; tot_ss_bar=0
for gal in good:
    g1=t1[gal]; g2=t2[gal]; pc=precomp[gal]
    Vf=g1['Vf']; Rd=max(g1['Rd'],0.01)
    MBH=10**(7.22+1.65*np.log10(max(Vf,10)/200))
    rc=C2*Rd; r0=rf2*Rd; mf=(MBH/M_ref)**(2/3)
    n=len(g2['R'])
    vpred=np.zeros(n)
    for i in range(n):
        vv0=v_LO(g2['R'][i],Vf,rc)
        r_m=g2['R'][i]*kpc_m; vbar_ms=pc['vbar'][i]*1e3
        M_bar_r=max(vbar_ms**2*r_m/G_N/Msun,1)
        q=M_bar_r/max(MBH,1)
        Psi=1+gg2*(4*q/(1+q)**2)**pq2*np.exp(-g2['R'][i]/max(r0,0.01))*mf
        vv=vv0*np.sqrt(max(Psi,0.01))
        vpred[i]=np.sqrt(pc['vbar_cst'][i]**2+vv**2)
    tot_ss+=np.sum((g2['Vo']-vpred)**2)
    tot_ss_bar+=np.sum((g2['Vo']-pc['vbar'])**2)
    tot_n+=len(g2['R'])
rmse_fb=np.sqrt(tot_ss/tot_n); rmse_bar=np.sqrt(tot_ss_bar/tot_n)

print(f"\n  Con β_gal = 2/3 FISSO (4 parametri liberi):")
print(f"  Ψ = 1 + {gg2:.1f} × [4q/(1+q)²]^{pq2:.2f} × exp(-r/{rf2:.2f}×Rd) × (M_BH/10⁷)^(2/3)")
print(f"  r_c = {C2:.1f} × Rd")
print(f"  RMSE = {rmse_fb:.2f} km/s (miglioramento: {rmse_bar/rmse_fb:.2f}x)")
print(f"  vs β_gal libero (0.258): RMSE = 19.78 km/s")
print(f"  Differenza: {rmse_fb-19.78:+.2f} km/s")

# ============================================================
# STEP 3: Il bulge come perturbazione SEPARATA
# ============================================================
print(f"\n\n{'='*70}")
print(f"STEP 3: BULGE come perturbazione separata dal disco")
print(f"{'='*70}")
print("""
  La diagnostica mostra che i peggiori errori sono nelle galassie 
  con bulge massiccio (Sb, Sa). Il bulge è una popolazione stellare
  DIVERSA dal disco: concentrata, sferoidale, spessa.
  
  Idea: il bulge perturba il vortice in modo diverso dal disco.
  Il disco è piatto come il vortice (coplanare → perturbazione minima).
  Il bulge è sferico → interseca il piano del vortice → perturbazione forte.
""")

def model_G(params, return_details=False):
    """Separate bulge perturbation"""
    C, gamma_g, pq, rf, eta_bul = params
    beta_g = 2.0/3.0
    if C<0.5 or gamma_g<0 or pq<0 or rf<0.1: return 1e20
    
    total=0; details=[]
    for gal in good:
        g1=t1[gal]; g2=t2[gal]; pc=precomp[gal]
        Vf=g1['Vf']; Rd=max(g1['Rd'],0.01)
        MBH=10**(7.22+1.65*np.log10(max(Vf,10)/200))
        rc=C*Rd; r0=rf*Rd; mf=(MBH/M_ref)**beta_g
        n=len(g2['R'])
        
        vpred=np.zeros(n)
        for i in range(n):
            vv0=v_LO(g2['R'][i],Vf,rc)
            
            # Disk+gas mass ratio
            r_m=g2['R'][i]*kpc_m
            vbar_ms=pc['vbar'][i]*1e3
            M_bar_r=max(vbar_ms**2*r_m/G_N/Msun,1)
            q=M_bar_r/max(MBH,1)
            
            # Standard interference
            Psi=1+gamma_g*(4*q/(1+q)**2)**pq*np.exp(-g2['R'][i]/max(r0,0.01))*mf
            
            # EXTRA bulge perturbation: bulge velocity contribution
            vbul = ML_bul * g2['Vb'][i]**2  # v²_bulge contribution
            vbar2_total = pc['vbar'][i]**2
            f_bul = vbul / max(vbar2_total, 1)  # fraction from bulge
            
            # Bulge perturbs the core size (spherical vs planar)
            Psi_bul = 1 + eta_bul * f_bul
            rc_eff = rc * Psi_bul
            
            vv0_eff = v_LO(g2['R'][i], Vf, rc_eff)
            vv = vv0_eff * np.sqrt(max(Psi, 0.01))
            vpred[i] = np.sqrt(pc['vbar_cst'][i]**2 + vv**2)
        
        total+=np.sum(((g2['Vo']-vpred)/np.maximum(g2['eV'],1))**2)
        if return_details:
            rmse=np.sqrt(np.mean((g2['Vo']-vpred)**2))
            rmse_bar=np.sqrt(np.mean((g2['Vo']-pc['vbar'])**2))
            details.append({'name':gal,'rmse':rmse,'rmse_bar':rmse_bar,'Vf':Vf,'T':g1['T']})
    
    if return_details: return total, details
    return total

print("Grid search modello G (5 param: C, γ, p, r₀, η_bul)...")
best_G=(1e30,[15,80,0.3,2,1])
for C in [10,15,20]:
    for gg in [50,100,150]:
        for pq in [0.2,0.5]:
            for rf in [1.5,2.5,4]:
                for eb in [-2,-1,0,1,2,3]:
                    c2=model_G([C,gg,pq,rf,eb])
                    if c2<best_G[0]: best_G=(c2,[C,gg,pq,rf,eb])

res_G=minimize(model_G, best_G[1], method='Nelder-Mead', options={'maxiter':20000})

chi2_G, details_G = model_G(res_G.x, return_details=True)
tot_ss_G=sum(d['rmse']**2*len(t2[d['name']]['R']) for d in details_G)
tot_n_G=sum(len(t2[d['name']]['R']) for d in details_G)
rmse_G=np.sqrt(tot_ss_G/tot_n_G)
C_G,gg_G,pq_G,rf_G,eb_G = res_G.x

print(f"\n  Modello G — con perturbazione bulge separata:")
print(f"  r_c_eff = {C_G:.1f}×Rd × [1 + {eb_G:.2f} × f_bulge(r)]")
print(f"  Ψ = 1 + {gg_G:.1f} × [4q/(1+q)²]^{pq_G:.2f} × exp(-r/{rf_G:.1f}×Rd) × (M_BH/10⁷)^(2/3)")
print(f"  β_gal = 2/3 (fisso)")
print(f"  Parametri: 5 (C, γ, p, r₀, η_bul)")
print(f"  RMSE = {rmse_G:.2f} km/s")
print(f"  Miglioramento: {rmse_bar/rmse_G:.2f}x")

# Per-type performance
tn={0:'S0',1:'Sa',2:'Sab',3:'Sb',4:'Sbc',5:'Sc',6:'Scd',7:'Sd',8:'Sdm',9:'Sm',10:'Im',11:'BCD'}
ts={}
for d in details_G:
    T=d['T']
    if T not in ts: ts[T]={'r':[],'rb':[]}
    ts[T]['r'].append(d['rmse']); ts[T]['rb'].append(d['rmse_bar'])

print(f"\n  Per tipo morfologico:")
print(f"  {'Tipo':>4s} {'N':>3s} {'RMSE':>6s} {'RMSE_bar':>8s} {'Migl':>6s}")
for T in sorted(ts.keys()):
    s=ts[T]; n=len(s['r'])
    imp=np.median(np.array(s['rb'])/np.maximum(np.array(s['r']),0.1))
    print(f"  {tn.get(T,'?'):>4s} {n:3d} {np.median(s['r']):6.1f} {np.median(s['rb']):8.1f} {imp:5.1f}x")

# ============================================================
# STEP 4: Cross-validation (il test definitivo)
# ============================================================
print(f"\n\n{'='*70}")
print(f"STEP 4: CROSS-VALIDATION 5-FOLD")
print(f"{'='*70}")

# Use best model so far
def evaluate_model(params, gal_list):
    C,gg,pq,rf,eb = params
    beta_g=2/3
    tot_ss=0; tot_n=0
    for gal in gal_list:
        g1=t1[gal]; g2=t2[gal]; pc=precomp[gal]
        Vf=g1['Vf']; Rd=max(g1['Rd'],0.01)
        MBH=10**(7.22+1.65*np.log10(max(Vf,10)/200))
        rc=C*Rd; r0=rf*Rd; mf=(MBH/M_ref)**beta_g
        n=len(g2['R'])
        vpred=np.zeros(n)
        for i in range(n):
            vbul=ML_bul*g2['Vb'][i]**2; f_bul=vbul/max(pc['vbar'][i]**2,1)
            rc_eff=rc*(1+eb*f_bul)
            vv0=v_LO(g2['R'][i],Vf,rc_eff)
            r_m=g2['R'][i]*kpc_m; vbar_ms=pc['vbar'][i]*1e3
            M_bar_r=max(vbar_ms**2*r_m/G_N/Msun,1)
            q=M_bar_r/max(MBH,1)
            Psi=1+gg*(4*q/(1+q)**2)**pq*np.exp(-g2['R'][i]/max(r0,0.01))*mf
            vv=vv0*np.sqrt(max(Psi,0.01))
            vpred[i]=np.sqrt(pc['vbar_cst'][i]**2+vv**2)
        tot_ss+=np.sum((g2['Vo']-vpred)**2); tot_n+=n
    return np.sqrt(tot_ss/tot_n) if tot_n>0 else 999

def fit_model(gal_list):
    def chi2(params):
        C,gg,pq,rf,eb=params
        beta_g=2/3
        if C<0.5 or gg<0 or pq<0 or rf<0.1: return 1e20
        total=0
        for gal in gal_list:
            g1=t1[gal]; g2=t2[gal]; pc=precomp[gal]
            Vf=g1['Vf']; Rd=max(g1['Rd'],0.01)
            MBH=10**(7.22+1.65*np.log10(max(Vf,10)/200))
            rc=C*Rd; r0=rf*Rd; mf=(MBH/M_ref)**beta_g
            n=len(g2['R'])
            vpred=np.zeros(n)
            for i in range(n):
                vbul=ML_bul*g2['Vb'][i]**2; f_bul=vbul/max(pc['vbar'][i]**2,1)
                rc_eff=rc*(1+eb*f_bul)
                vv0=v_LO(g2['R'][i],Vf,rc_eff)
                r_m=g2['R'][i]*kpc_m; vbar_ms=pc['vbar'][i]*1e3
                M_bar_r=max(vbar_ms**2*r_m/G_N/Msun,1)
                q=M_bar_r/max(MBH,1)
                Psi=1+gg*(4*q/(1+q)**2)**pq*np.exp(-g2['R'][i]/max(r0,0.01))*mf
                vv=vv0*np.sqrt(max(Psi,0.01))
                vpred[i]=np.sqrt(pc['vbar_cst'][i]**2+vv**2)
            total+=np.sum(((g2['Vo']-vpred)/np.maximum(g2['eV'],1))**2)
        return total
    res=minimize(chi2, res_G.x, method='Nelder-Mead', options={'maxiter':10000})
    return res.x

np.random.seed(42)
gal_arr = np.array(good)
np.random.shuffle(gal_arr)
folds = np.array_split(gal_arr, 5)

print(f"\n  5-fold cross-validation:")
train_rmses = []; test_rmses = []
for fold_idx in range(5):
    test_gals = list(folds[fold_idx])
    train_gals = [g for i,f in enumerate(folds) if i!=fold_idx for g in f]
    
    # Fit on train
    params = fit_model(train_gals)
    
    # Evaluate on both
    rmse_train = evaluate_model(params, train_gals)
    rmse_test = evaluate_model(params, test_gals)
    train_rmses.append(rmse_train)
    test_rmses.append(rmse_test)
    
    print(f"    Fold {fold_idx+1}: train={rmse_train:.2f}, test={rmse_test:.2f} km/s  (N_test={len(test_gals)})")

print(f"\n  Media train: {np.mean(train_rmses):.2f} ± {np.std(train_rmses):.2f} km/s")
print(f"  Media test:  {np.mean(test_rmses):.2f} ± {np.std(test_rmses):.2f} km/s")
print(f"  Differenza:  {np.mean(test_rmses)-np.mean(train_rmses):+.2f} km/s")
print(f"  → {'NO overfitting' if np.mean(test_rmses)-np.mean(train_rmses) < 3 else '⚠ POSSIBILE overfitting'}")

# ============================================================
# FINAL TABLE
# ============================================================
print(f"\n\n{'='*70}")
print(f"TABELLA FINALE — MODELLO CST GALATTICO OTTIMIZZATO")
print(f"{'='*70}")
print(f"""
  ══════════════════════════════════════════════════════════════════
  FORMULA CST GALATTICA COMPLETA (Modello G)
  ══════════════════════════════════════════════════════════════════
  
  v²_tot(r) = v²_bar(r) × (1 + ε_PDMF) + v²_vortice(r) × Ψ_gal(r)
  
  v_vortice(r) = V_flat × LO(r / r_c_eff(r))
  
  r_c_eff(r) = C × R_disk × [1 + η_bul × f_bul(r)]
  f_bul(r) = v²_bul(r) / v²_bar(r)
  
  Ψ_gal(r) = 1 + γ × [4q/(1+q)²]^p × exp(-r/r₀×R_disk) × (M_BH/10⁷)^(2/3)
  q(r) = M_bar(<r) / M_BH
  
  PARAMETRI:
    Fissi (dalla teoria):
      ε_PDMF = 0.215 (popolazione stellare)
      β_gal  = 2/3   (esponente universale CST)
      v_max  = V_flat (il vortice è la curva piatta)
    
    Globali (5, fittati su 129 galassie):
      C      = {C_G:.1f}   (raggio core vortice in unità di Rd)
      γ      = {gg_G:.1f}  (intensità interferenza stellare)
      p      = {pq_G:.2f}  (esponente rapporto massa)
      r₀     = {rf_G:.1f}   (scala radiale interferenza in unità di Rd)
      η_bul  = {eb_G:.2f}  (perturbazione aggiuntiva del bulge)
    
    Per galassia: ZERO
  
  RISULTATO:
    RMSE = {rmse_G:.2f} km/s  (miglioramento {rmse_bar/rmse_G:.2f}x vs barioni)
    Cross-validation: train {np.mean(train_rmses):.1f} ± {np.std(train_rmses):.1f}, test {np.mean(test_rmses):.1f} ± {np.std(test_rmses):.1f}
    
  CONFRONTO:
    Solo barioni:       61.9 km/s (1.0x)
    MOND (1 param):     ~25  km/s (~2.5x)
    CST-G (5 param):    {rmse_G:.1f} km/s ({rmse_bar/rmse_G:.1f}x)
    CST individuale:    18.4 km/s (3.4x)  [258 param]
  ══════════════════════════════════════════════════════════════════
""")

