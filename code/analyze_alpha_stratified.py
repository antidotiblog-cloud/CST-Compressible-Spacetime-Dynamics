#!/usr/bin/env python3
"""
ANALYZE_ALPHA_STRATIFIED.PY - α(M) Mass-Dependent Analysis
============================================================

Tests whether coupling strength α varies with system mass.

Critical Question:
------------------
In original fits, α was used as GLOBAL constant (0.279) for ALL 21,565 systems.
But physically, does α depend on mass M?

Theory predicts:
- Compact objects (M ~ M_☉): α ~ 0.28 (observed)
- Extended structures (M >> M_☉): α_cosmo ~ 0.07 (theoretical)

This script:
1. Stratifies binaries by mass bins
2. Fits α separately for each bin
3. Tests for systematic mass dependence α(M)

Author: Michele Vizzutti
Date: February 10, 2026
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# Mass bins for stratification (solar masses)
MASS_BINS = [
    (0.5, 0.8),   # Low mass
    (0.8, 1.2),   # Solar-like
    (1.2, 1.8),   # Intermediate
    (1.8, 3.0)    # High mass
]

BETA_FIXED = 0.685  # Fix β from global fit
Z_TRANS = 30.0

# ============================================================================
# CORE FUNCTIONS (same as fit_master.py)
# ============================================================================

def age_to_redshift(age_Gyr, age_universe=13.8):
    t_formation = age_universe - age_Gyr
    t_formation = np.maximum(t_formation, 0.01)
    z = (age_universe / t_formation)**(2/3) - 1
    return np.maximum(z, 0)

def hubble_ratio(z, Omega_m=0.315, Omega_L=0.685):
    return np.sqrt(Omega_m * (1+z)**3 + Omega_L)

def transition_function(z, z_trans=Z_TRANS):
    H_ratio = hubble_ratio(z)
    S = 1.0 / (1.0 + (z/z_trans)**3)
    return H_ratio * S

def w_mass(M):
    return np.exp(-np.abs(M - 1))

def f_q(q):
    return 4*q / (1+q)**2

def f_a(a, a0=0.5):
    return np.exp(-a/a0)

def f_M(M_tot, beta=BETA_FIXED):
    return M_tot**beta

def Psi_interference(q, a, M_tot, gamma=8.0, a0=0.5):
    return 1 + gamma * f_q(q) * f_a(a, a0) * f_M(M_tot, BETA_FIXED)

def v_ratio_model(M_tot, z, q, a, alpha, gamma=8.0, a0=0.5):
    """
    Velocity ratio model with α as free parameter.
    β, γ, a₀ fixed from global fit.
    """
    w = w_mass(M_tot)
    f_z = transition_function(z)
    Psi = Psi_interference(q, a, M_tot, gamma, a0)
    alpha_eff = alpha * Psi
    G_ratio = 1 + (1 - w) * alpha_eff * f_z
    return np.sqrt(G_ratio)

# ============================================================================
# DATA LOADING
# ============================================================================

def load_all_binaries(directory='../02_DATI/binaries'):
    """Load and combine all binary CSV files."""
    import glob
    
    files = glob.glob(f'{directory}/*.csv')
    dfs = []
    
    for f in files:
        try:
            df = pd.read_csv(f)
            dfs.append(df)
        except Exception as e:
            print(f"⚠️  Skipping {f}: {e}")
    
    df_combined = pd.concat(dfs, ignore_index=True)
    
    # Ensure required columns
    required = ['M_total', 'q', 'SemiMajor_AU', 'age_Gyr', 'v_ratio']
    missing = [col for col in required if col not in df_combined.columns]
    
    if missing:
        print(f"❌ Missing columns: {missing}")
        # Attempt to derive
        if 'M1' in df_combined.columns and 'M2' in df_combined.columns:
            df_combined['M_total'] = df_combined['M1'] + df_combined['M2']
            df_combined['q'] = df_combined['M2'] / df_combined['M1']
        
        if 'age_Gyr' not in df_combined.columns and 'Age' in df_combined.columns:
            df_combined['age_Gyr'] = df_combined['Age']
    
    # Calculate cosmology
    df_combined['z_formation'] = age_to_redshift(df_combined['age_Gyr'])
    df_combined['f_z'] = transition_function(df_combined['z_formation'])
    
    return df_combined

# ============================================================================
# STRATIFIED FITTING
# ============================================================================

def fit_alpha_for_bin(df_bin, mass_range):
    """
    Fit α for a specific mass bin.
    
    Parameters:
    -----------
    df_bin : pandas.DataFrame
        Binary data
    mass_range : tuple
        (M_min, M_max) in solar masses
    
    Returns:
    --------
    alpha : float
        Fitted coupling strength
    alpha_err : float
        Standard error
    N : int
        Number of systems in bin
    R2 : float
        Fit quality
    """
    M_min, M_max = mass_range
    
    # Filter by mass
    mask = (df_bin['M_total'] >= M_min) & (df_bin['M_total'] < M_max)
    df_sub = df_bin[mask].copy()
    N = len(df_sub)
    
    if N < 10:
        print(f"  ⚠️  Only {N} systems in bin [{M_min:.1f}, {M_max:.1f}] M_☉ - SKIP")
        return None, None, N, None
    
    # Extract data
    M_tot = df_sub['M_total'].values
    z = df_sub['z_formation'].values
    q = df_sub['q'].values
    a = df_sub['SemiMajor_AU'].values
    v_obs = df_sub['v_ratio'].values
    
    # Fit function (only α free)
    def model(dummy, alpha):
        return v_ratio_model(M_tot, z, q, a, alpha, gamma=8.0, a0=0.5)
    
    x_dummy = np.arange(N)
    
    try:
        popt, pcov = curve_fit(
            model, x_dummy, v_obs,
            p0=[0.279],
            bounds=([0], [1]),
            maxfev=10000
        )
        
        alpha = popt[0]
        alpha_err = np.sqrt(pcov[0,0])
        
        # Calculate R²
        v_pred = model(x_dummy, alpha)
        ss_res = np.sum((v_obs - v_pred)**2)
        ss_tot = np.sum((v_obs - np.mean(v_obs))**2)
        R2 = 1 - (ss_res / ss_tot)
        
        return alpha, alpha_err, N, R2
        
    except Exception as e:
        print(f"  ❌ Fit failed: {e}")
        return None, None, N, None

# ============================================================================
# ANALYSIS
# ============================================================================

def analyze_alpha_stratified(df):
    """
    Main stratified analysis.
    
    Returns:
    --------
    results : pandas.DataFrame
        Table with columns: M_range, M_center, alpha, alpha_err, N, R2
    """
    print("\n" + "="*70)
    print("α(M) STRATIFIED ANALYSIS")
    print("="*70)
    print()
    
    results = []
    
    for M_min, M_max in MASS_BINS:
        M_center = (M_min + M_max) / 2
        print(f"Mass bin: [{M_min:.1f}, {M_max:.1f}] M_☉  (center: {M_center:.2f})")
        
        alpha, alpha_err, N, R2 = fit_alpha_for_bin(df, (M_min, M_max))
        
        if alpha is not None:
            print(f"  N = {N}")
            print(f"  α = {alpha:.3f} ± {alpha_err:.3f}")
            print(f"  R² = {R2:.4f}")
            print()
            
            results.append({
                'M_min': M_min,
                'M_max': M_max,
                'M_center': M_center,
                'alpha': alpha,
                'alpha_err': alpha_err,
                'N': N,
                'R2': R2
            })
    
    return pd.DataFrame(results)

# ============================================================================
# STATISTICAL TESTS
# ============================================================================

def test_alpha_dependence(results_df):
    """
    Test for systematic α(M) trend.
    
    Null hypothesis: α is constant (no mass dependence)
    Alternative: α = α₀ + β_M × log(M/M_☉)
    """
    print("\n" + "="*70)
    print("STATISTICAL TEST: α vs M")
    print("="*70)
    print()
    
    M_center = results_df['M_center'].values
    alpha = results_df['alpha'].values
    alpha_err = results_df['alpha_err'].values
    
    # Weighted linear fit
    weights = 1 / alpha_err**2
    
    # Log(M)
    log_M = np.log10(M_center)
    
    # Weighted mean (null hypothesis)
    alpha_mean = np.average(alpha, weights=weights)
    alpha_mean_err = 1 / np.sqrt(np.sum(weights))
    
    print(f"Weighted mean: α = {alpha_mean:.3f} ± {alpha_mean_err:.3f}")
    
    # Linear trend
    def linear_model(log_M, alpha0, beta_M):
        return alpha0 + beta_M * log_M
    
    try:
        from scipy.optimize import curve_fit
        popt, pcov = curve_fit(
            linear_model, log_M, alpha,
            sigma=alpha_err,
            p0=[0.279, 0],
            absolute_sigma=True
        )
        
        alpha0, beta_M = popt
        alpha0_err = np.sqrt(pcov[0,0])
        beta_M_err = np.sqrt(pcov[1,1])
        
        print(f"\nLinear fit: α = α₀ + β_M × log(M/M_☉)")
        print(f"  α₀ = {alpha0:.3f} ± {alpha0_err:.3f}")
        print(f"  β_M = {beta_M:.3f} ± {beta_M_err:.3f}")
        
        # Significance test
        t_statistic = beta_M / beta_M_err
        from scipy.stats import t as t_dist
        p_value = 2 * (1 - t_dist.cdf(abs(t_statistic), df=len(alpha)-2))
        
        print(f"\nSignificance test:")
        print(f"  t-statistic = {t_statistic:.2f}")
        print(f"  p-value = {p_value:.4f}")
        
        if p_value < 0.05:
            print(f"  ✅ SIGNIFICANT mass dependence (p < 0.05)")
        else:
            print(f"  ❌ NO significant mass dependence (p ≥ 0.05)")
        
        # Chi-squared goodness of fit
        alpha_pred = linear_model(log_M, alpha0, beta_M)
        chi2 = np.sum(((alpha - alpha_pred) / alpha_err)**2)
        dof = len(alpha) - 2
        chi2_reduced = chi2 / dof
        
        print(f"\nGoodness of fit:")
        print(f"  χ² = {chi2:.2f}")
        print(f"  dof = {dof}")
        print(f"  χ²/dof = {chi2_reduced:.2f}")
        
    except Exception as e:
        print(f"❌ Linear fit failed: {e}")
        beta_M = None
        p_value = None
    
    return alpha_mean, beta_M, p_value

# ============================================================================
# PLOTTING
# ============================================================================

def plot_alpha_vs_mass(results_df, output_file=None):
    """
    Plot α vs mass with error bars.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    M_center = results_df['M_center'].values
    alpha = results_df['alpha'].values
    alpha_err = results_df['alpha_err'].values
    
    # Data points
    ax.errorbar(M_center, alpha, yerr=alpha_err, 
                fmt='o', markersize=8, capsize=5, capthick=2,
                color='blue', label='Fitted α per bin')
    
    # Global mean
    alpha_global = 0.279
    ax.axhline(alpha_global, color='red', linestyle='--', lw=2,
               label=f'Global α = {alpha_global:.3f}')
    
    # Best fit line
    log_M = np.log10(M_center)
    
    try:
        from scipy.optimize import curve_fit
        def linear_model(log_M, alpha0, beta_M):
            return alpha0 + beta_M * log_M
        
        popt, _ = curve_fit(
            linear_model, log_M, alpha,
            sigma=alpha_err,
            absolute_sigma=True
        )
        
        alpha0, beta_M = popt
        
        M_plot = np.logspace(np.log10(M_center.min()), np.log10(M_center.max()), 100)
        alpha_plot = linear_model(np.log10(M_plot), alpha0, beta_M)
        
        ax.plot(M_plot, alpha_plot, 'g-', lw=2,
                label=f'Linear fit: α = {alpha0:.3f} + {beta_M:.3f}×log(M)')
    
    except:
        pass
    
    ax.set_xlabel('Total Mass (M_☉)', fontsize=12)
    ax.set_ylabel('Coupling Strength α', fontsize=12)
    ax.set_title('Mass Dependence of α', fontsize=14, weight='bold')
    ax.set_xscale('log')
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\n✅ Figure saved to {output_file}")
    else:
        plt.show()
    
    plt.close()

# ============================================================================
# MAIN
# ============================================================================

def main():
    print("="*70)
    print("α(M) STRATIFIED ANALYSIS - MASS DEPENDENCE TEST")
    print("="*70)
    print()
    
    # Load data
    print("1. Loading binary data...")
    df = load_all_binaries()
    print(f"   Loaded {len(df)} binary systems")
    print(f"   Mass range: {df['M_total'].min():.2f} - {df['M_total'].max():.2f} M_☉")
    print()
    
    # Stratified analysis
    results_df = analyze_alpha_stratified(df)
    
    # Statistical test
    alpha_mean, beta_M, p_value = test_alpha_dependence(results_df)
    
    # Plot
    plot_alpha_vs_mass(results_df, '../06_ANALISI/figures/alpha_vs_mass.png')
    
    # Save results
    results_df.to_csv('../04_RISULTATI/alpha_stratified.csv', index=False)
    print(f"\n✅ Results saved to ../04_RISULTATI/alpha_stratified.csv")
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print()
    print(f"Global α (all masses): {0.279:.3f} ± {0.012:.3f}")
    print(f"Weighted mean α (stratified): {alpha_mean:.3f}")
    
    if beta_M is not None and p_value is not None:
        if p_value < 0.05:
            print(f"\n✅ SIGNIFICANT mass dependence detected!")
            print(f"   β_M = {beta_M:.3f} (slope of α vs log M)")
            print(f"   p = {p_value:.4f}")
        else:
            print(f"\n❌ NO significant mass dependence")
            print(f"   β_M = {beta_M:.3f} ± error")
            print(f"   p = {p_value:.4f} (> 0.05)")
            print(f"\n   Conclusion: α appears CONSTANT across mass range")
    
    print()

if __name__ == "__main__":
    main()
