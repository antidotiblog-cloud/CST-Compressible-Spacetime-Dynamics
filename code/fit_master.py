#!/usr/bin/env python3
"""
FIT_MASTER.PY - Unified Fitting Script for CST Theory
======================================================

Fits G_eff(M,z) formula to multi-scale astronomical systems:
- Exoplanets:  N = 4,585  (compact, single-star systems)
- Binaries:    N = 16,980 (compact, two-body systems with interference)

Total systems: 21,565

Formula:
--------
EXOPLANETS (single-star):
    G_eff/G = 1 + α(M/M_☉)^β × f(z)

BINARIES (interference):
    G_eff/G = 1 + α(M/M_☉)^β × Ψ(q,a,M) × f(z)
    
    where Ψ(q,a,M) = 1 + γ × f_q(q) × f_a(a) × f_M(M)

Transition function:
    f(z) = √[Ω_m(1+z)³ + Ω_Λ] / [1 + (z/z_trans)³]

Author: Michele Vizzutti
Date: February 10, 2026
Version: 1.0
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, differential_evolution
from scipy import stats
import pickle
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# COSMOLOGICAL PARAMETERS (Planck 2018)
# ============================================================================

OMEGA_M = 0.315  # Matter density today
OMEGA_L = 0.685  # Dark energy density today
H0 = 67.4  # Hubble constant km/s/Mpc
AGE_UNIVERSE = 13.8  # Gyr

# CST Parameters (from previous fits)
ALPHA_GLOBAL = 0.279  # ± 0.012
BETA_GLOBAL = 0.685   # ± 0.018
Z_TRANS = 30.0        # ± 10

# ============================================================================
# CORE FUNCTIONS
# ============================================================================

def age_to_redshift(age_Gyr, age_universe=AGE_UNIVERSE):
    """
    Convert stellar age to formation redshift.
    
    Uses approximation: z = (t₀/t_form)^(2/3) - 1
    Valid for z < 2-3, error ~5-15%
    
    Parameters:
    -----------
    age_Gyr : float or array
        Stellar age in Gyr
    age_universe : float
        Age of universe in Gyr (default 13.8)
    
    Returns:
    --------
    z : float or array
        Formation redshift
    """
    t_formation = age_universe - age_Gyr
    t_formation = np.maximum(t_formation, 0.01)  # Avoid division by zero
    z = (age_universe / t_formation)**(2/3) - 1
    return np.maximum(z, 0)

def hubble_ratio(z, Omega_m=OMEGA_M, Omega_L=OMEGA_L):
    """
    Compute H(z)/H₀ for flat ΛCDM cosmology.
    
    Parameters:
    -----------
    z : float or array
        Redshift
    Omega_m : float
        Matter density parameter
    Omega_L : float
        Dark energy density parameter
    
    Returns:
    --------
    H_ratio : float or array
        H(z)/H₀
    """
    return np.sqrt(Omega_m * (1+z)**3 + Omega_L)

def transition_function(z, z_trans=Z_TRANS, n=3):
    """
    Structure formation transition function.
    
    Suppresses G_eff at high-z (before structures form).
    Preserves BBN (z~10⁹) and CMB (z~1100).
    
    Parameters:
    -----------
    z : float or array
        Redshift
    z_trans : float
        Transition redshift (default 30)
    n : int
        Sharpness exponent (default 3)
    
    Returns:
    --------
    f : float or array
        Transition function value
    """
    H_ratio = hubble_ratio(z)
    S = 1.0 / (1.0 + (z/z_trans)**n)
    return H_ratio * S

def w_mass(M, M_sun=1.0):
    """
    Mass weighting function.
    
    Parameters:
    -----------
    M : float or array
        Mass in solar masses
    M_sun : float
        Reference mass (default 1.0)
    
    Returns:
    --------
    w : float or array
        Weight (0 to 1)
    """
    return np.exp(-np.abs(M/M_sun - 1))

# ============================================================================
# INTERFERENCE FUNCTIONS (BINARIES)
# ============================================================================

def f_q(q):
    """Mass ratio symmetry factor."""
    return 4*q / (1+q)**2

def f_a(a, a0):
    """Separation factor."""
    return np.exp(-a/a0)

def f_M(M_tot, beta):
    """Mass scaling factor."""
    return M_tot**beta

def Psi_interference(q, a, M_tot, gamma, a0, beta):
    """
    Binary interference amplifier.
    
    Parameters:
    -----------
    q : float or array
        Mass ratio M₂/M₁
    a : float or array
        Separation in AU
    M_tot : float or array
        Total mass in solar masses
    gamma : float
        Interference coupling
    a0 : float
        Separation scale in AU
    beta : float
        Mass scaling exponent
    
    Returns:
    --------
    Psi : float or array
        Amplification factor (≥ 1)
    """
    return 1 + gamma * f_q(q) * f_a(a, a0) * f_M(M_tot, beta)

# ============================================================================
# G_EFF MODELS
# ============================================================================

def G_eff_exoplanets(M, z, alpha=ALPHA_GLOBAL, beta=BETA_GLOBAL):
    """
    G_eff for exoplanet systems (single-star).
    
    Parameters:
    -----------
    M : float or array
        Stellar mass in solar masses
    z : float or array
        Formation redshift
    alpha : float
        Coupling strength
    beta : float
        Mass scaling exponent
    
    Returns:
    --------
    G_ratio : float or array
        G_eff/G_N
    """
    w = w_mass(M)
    f_z = transition_function(z)
    return 1 + (1 - w) * alpha * (M**beta) * f_z

def G_eff_binaries(M_tot, z, q, a, alpha=ALPHA_GLOBAL, beta=BETA_GLOBAL, 
                   gamma=8.0, a0=0.5):
    """
    G_eff for binary star systems (with interference).
    
    Parameters:
    -----------
    M_tot : float or array
        Total mass in solar masses
    z : float or array
        Formation redshift
    q : float or array
        Mass ratio M₂/M₁
    a : float or array
        Separation in AU
    alpha : float
        Coupling strength
    beta : float
        Mass scaling exponent
    gamma : float
        Interference coupling
    a0 : float
        Separation scale in AU
    
    Returns:
    --------
    G_ratio : float or array
        G_eff/G_N
    """
    w = w_mass(M_tot)
    f_z = transition_function(z)
    Psi = Psi_interference(q, a, M_tot, gamma, a0, beta)
    alpha_eff = alpha * Psi
    return 1 + (1 - w) * alpha_eff * f_z

# ============================================================================
# VELOCITY RATIO (OBSERVABLE)
# ============================================================================

def v_ratio_exoplanets(M, z, alpha, beta):
    """Velocity ratio for exoplanets."""
    G_ratio = G_eff_exoplanets(M, z, alpha, beta)
    return np.sqrt(G_ratio)

def v_ratio_binaries(M_tot, z, q, a, alpha, beta, gamma, a0):
    """Velocity ratio for binaries."""
    G_ratio = G_eff_binaries(M_tot, z, q, a, alpha, beta, gamma, a0)
    return np.sqrt(G_ratio)

# ============================================================================
# DATA LOADING
# ============================================================================

def load_exoplanets(filepath):
    """
    Load exoplanet dataset.
    
    Expected columns:
    - st_mass: stellar mass (M_☉)
    - st_age: stellar age (Gyr)
    - pl_orbper: orbital period (days)
    - pl_orbsmax: semi-major axis (AU)
    - [other columns...]
    
    Returns:
    --------
    df : pandas.DataFrame
    """
    df = pd.read_csv(filepath)
    
    # Calculate cosmology
    df['z_formation'] = age_to_redshift(df['st_age'])
    df['H_ratio'] = hubble_ratio(df['z_formation'])
    df['f_z'] = transition_function(df['z_formation'])
    
    # Calculate velocities
    G = 6.674e-11  # SI
    M_sun = 1.989e30
    AU = 1.496e11
    
    # Keplerian velocity
    df['v_kepler'] = np.sqrt(G * df['st_mass'] * M_sun / (df['pl_orbsmax'] * AU))
    
    # Observed velocity (from period and semi-major axis)
    df['v_observed'] = 2 * np.pi * df['pl_orbsmax'] * AU / (df['pl_orbper'] * 86400)
    
    # Velocity ratio
    df['v_ratio'] = df['v_observed'] / df['v_kepler']
    
    return df

def load_binaries(filepath):
    """
    Load binary star dataset.
    
    Expected columns:
    - M1: primary mass (M_☉)
    - M2: secondary mass (M_☉) or q = M2/M1
    - Period: orbital period (days)
    - SemiMajor_AU: semi-major axis (AU)
    - age_Gyr: system age (Gyr)
    - [other columns...]
    
    Returns:
    --------
    df : pandas.DataFrame
    """
    df = pd.read_csv(filepath)
    
    # Calculate total mass and mass ratio
    if 'M_total' not in df.columns:
        df['M_total'] = df['M1'] + df['M2']
    
    if 'q' not in df.columns:
        df['q'] = df['M2'] / df['M1']
    
    # Calculate cosmology
    df['z_formation'] = age_to_redshift(df['age_Gyr'])
    df['H_ratio'] = hubble_ratio(df['z_formation'])
    df['f_z'] = transition_function(df['z_formation'])
    
    # Calculate velocities (similar to exoplanets)
    G = 6.674e-11
    M_sun = 1.989e30
    AU = 1.496e11
    
    df['v_kepler'] = np.sqrt(G * df['M_total'] * M_sun / (df['SemiMajor_AU'] * AU))
    df['v_observed'] = 2 * np.pi * df['SemiMajor_AU'] * AU / (df['Period'] * 86400)
    df['v_ratio'] = df['v_observed'] / df['v_kepler']
    
    return df

# ============================================================================
# FITTING FUNCTIONS
# ============================================================================

def fit_exoplanets(df, initial_guess=[0.279, 0.685]):
    """
    Fit α, β for exoplanet dataset.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Exoplanet data with columns: st_mass, z_formation, v_ratio
    initial_guess : list
        Initial [α, β]
    
    Returns:
    --------
    popt : array
        Fitted parameters [α, β]
    pcov : array
        Covariance matrix
    """
    # Prepare data
    M = df['st_mass'].values
    z = df['z_formation'].values
    v_obs = df['v_ratio'].values
    
    # Fit function wrapper
    def model(dummy, alpha, beta):
        return v_ratio_exoplanets(M, z, alpha, beta)
    
    # Dummy x for curve_fit
    x_dummy = np.arange(len(M))
    
    # Fit
    popt, pcov = curve_fit(
        model, x_dummy, v_obs,
        p0=initial_guess,
        bounds=([0, 0], [1, 2]),
        maxfev=10000
    )
    
    return popt, pcov

def fit_binaries(df, initial_guess=[0.279, 0.685, 8.0, 0.5]):
    """
    Fit α, β, γ, a₀ for binary dataset.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Binary data with columns: M_total, z_formation, q, SemiMajor_AU, v_ratio
    initial_guess : list
        Initial [α, β, γ, a₀]
    
    Returns:
    --------
    popt : array
        Fitted parameters [α, β, γ, a₀]
    pcov : array
        Covariance matrix
    """
    # Prepare data
    M_tot = df['M_total'].values
    z = df['z_formation'].values
    q = df['q'].values
    a = df['SemiMajor_AU'].values
    v_obs = df['v_ratio'].values
    
    # Fit function wrapper
    def model(dummy, alpha, beta, gamma, a0):
        return v_ratio_binaries(M_tot, z, q, a, alpha, beta, gamma, a0)
    
    # Dummy x for curve_fit
    x_dummy = np.arange(len(M_tot))
    
    # Fit
    popt, pcov = curve_fit(
        model, x_dummy, v_obs,
        p0=initial_guess,
        bounds=([0, 0, 0, 0.01], [1, 2, 20, 2.0]),
        maxfev=10000
    )
    
    return popt, pcov

# ============================================================================
# EVALUATION AND STATISTICS
# ============================================================================

def evaluate_fit(y_obs, y_pred):
    """
    Evaluate fit quality.
    
    Returns:
    --------
    results : dict
        {
            'R2': R-squared,
            'RMSE': root mean square error,
            'pearson_r': Pearson correlation,
            'p_value': significance,
            'residuals': y_obs - y_pred,
            'residuals_std': standard deviation of residuals
        }
    """
    residuals = y_obs - y_pred
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y_obs - np.mean(y_obs))**2)
    R2 = 1 - (ss_res / ss_tot)
    
    RMSE = np.sqrt(np.mean(residuals**2))
    
    r, p = stats.pearsonr(y_obs, y_pred)
    
    return {
        'R2': R2,
        'RMSE': RMSE,
        'pearson_r': r,
        'p_value': p,
        'residuals': residuals,
        'residuals_std': np.std(residuals)
    }

# ============================================================================
# PLOTTING
# ============================================================================

def plot_fit_results(df, y_pred, system_type='exoplanets', output_file=None):
    """
    Create diagnostic plots.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Data with v_ratio
    y_pred : array
        Predicted v_ratio
    system_type : str
        'exoplanets' or 'binaries'
    output_file : str
        Save figure to file (optional)
    """
    y_obs = df['v_ratio'].values
    residuals = y_obs - y_pred
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. Observed vs Predicted
    ax = axes[0, 0]
    ax.scatter(y_pred, y_obs, alpha=0.5, s=20)
    ax.plot([y_obs.min(), y_obs.max()], [y_obs.min(), y_obs.max()], 
            'r--', lw=2, label='Perfect fit')
    ax.set_xlabel('Predicted v_ratio')
    ax.set_ylabel('Observed v_ratio')
    ax.set_title(f'{system_type.capitalize()}: Observed vs Predicted')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # 2. Residuals vs Predicted
    ax = axes[0, 1]
    ax.scatter(y_pred, residuals, alpha=0.5, s=20)
    ax.axhline(0, color='r', linestyle='--', lw=2)
    ax.set_xlabel('Predicted v_ratio')
    ax.set_ylabel('Residuals')
    ax.set_title('Residuals vs Predicted')
    ax.grid(alpha=0.3)
    
    # 3. Histogram of Residuals
    ax = axes[1, 0]
    ax.hist(residuals, bins=50, alpha=0.7, edgecolor='black')
    ax.axvline(0, color='r', linestyle='--', lw=2)
    ax.set_xlabel('Residuals')
    ax.set_ylabel('Frequency')
    ax.set_title('Distribution of Residuals')
    ax.grid(alpha=0.3)
    
    # 4. Enhancement vs Mass
    ax = axes[1, 1]
    if system_type == 'exoplanets':
        M = df['st_mass'].values
    else:
        M = df['M_total'].values
    enhancement = (y_obs - 1) * 100
    ax.scatter(M, enhancement, alpha=0.5, s=20)
    ax.set_xlabel('Mass (M_☉)')
    ax.set_ylabel('Enhancement (%)')
    ax.set_title('Enhancement vs Mass')
    ax.grid(alpha=0.3)
    ax.set_xscale('log')
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {output_file}")
    else:
        plt.show()
    
    plt.close()

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """
    Main execution function.
    
    Usage:
    ------
    python fit_master.py
    
    Modify filepaths below to point to your data.
    """
    
    print("="*70)
    print("CST UNIFIED FITTING SCRIPT")
    print("="*70)
    print()
    
    # ========================================================================
    # 1. LOAD EXOPLANETS
    # ========================================================================
    print("1. Loading exoplanet data...")
    try:
        df_exo = load_exoplanets('../02_DATI/exoplanets/exoplanet_data.csv')
        print(f"   Loaded {len(df_exo)} exoplanet systems")
    except FileNotFoundError:
        print("   ❌ Exoplanet CSV not found!")
        print("   Run: python download_nasa_exoplanets.py first")
        df_exo = None
    
    # ========================================================================
    # 2. LOAD BINARIES
    # ========================================================================
    print("\n2. Loading binary data...")
    try:
        # Concatenate all binary CSV files
        import glob
        binary_files = glob.glob('../02_DATI/binaries/*.csv')
        dfs = []
        for f in binary_files:
            dfs.append(load_binaries(f))
        df_bin = pd.concat(dfs, ignore_index=True)
        print(f"   Loaded {len(df_bin)} binary systems from {len(binary_files)} files")
    except Exception as e:
        print(f"   ❌ Binary loading failed: {e}")
        df_bin = None
    
    # ========================================================================
    # 3. FIT EXOPLANETS
    # ========================================================================
    if df_exo is not None:
        print("\n3. Fitting exoplanets...")
        popt_exo, pcov_exo = fit_exoplanets(df_exo)
        alpha_exo, beta_exo = popt_exo
        alpha_err = np.sqrt(pcov_exo[0,0])
        beta_err = np.sqrt(pcov_exo[1,1])
        
        print(f"   α = {alpha_exo:.3f} ± {alpha_err:.3f}")
        print(f"   β = {beta_exo:.3f} ± {beta_err:.3f}")
        
        # Predict and evaluate
        y_pred_exo = v_ratio_exoplanets(df_exo['st_mass'].values, 
                                        df_exo['z_formation'].values,
                                        alpha_exo, beta_exo)
        results_exo = evaluate_fit(df_exo['v_ratio'].values, y_pred_exo)
        
        print(f"   R² = {results_exo['R2']:.4f}")
        print(f"   RMSE = {results_exo['RMSE']:.4f}")
        print(f"   Pearson r = {results_exo['pearson_r']:.4f} (p < {results_exo['p_value']:.2e})")
        
        # Plot
        plot_fit_results(df_exo, y_pred_exo, 'exoplanets', 
                        '../06_ANALISI/figures/fit_exoplanets.png')
        
        # Save results
        results_exo['parameters'] = {'alpha': alpha_exo, 'beta': beta_exo}
        with open('../04_RISULTATI/fit_exoplanets.pkl', 'wb') as f:
            pickle.dump(results_exo, f)
    
    # ========================================================================
    # 4. FIT BINARIES
    # ========================================================================
    if df_bin is not None:
        print("\n4. Fitting binaries...")
        popt_bin, pcov_bin = fit_binaries(df_bin)
        alpha_bin, beta_bin, gamma_bin, a0_bin = popt_bin
        
        print(f"   α = {alpha_bin:.3f} ± {np.sqrt(pcov_bin[0,0]):.3f}")
        print(f"   β = {beta_bin:.3f} ± {np.sqrt(pcov_bin[1,1]):.3f}")
        print(f"   γ = {gamma_bin:.3f} ± {np.sqrt(pcov_bin[2,2]):.3f}")
        print(f"   a₀ = {a0_bin:.3f} ± {np.sqrt(pcov_bin[3,3]):.3f}")
        
        # Predict and evaluate
        y_pred_bin = v_ratio_binaries(df_bin['M_total'].values,
                                      df_bin['z_formation'].values,
                                      df_bin['q'].values,
                                      df_bin['SemiMajor_AU'].values,
                                      alpha_bin, beta_bin, gamma_bin, a0_bin)
        results_bin = evaluate_fit(df_bin['v_ratio'].values, y_pred_bin)
        
        print(f"   R² = {results_bin['R2']:.4f}")
        print(f"   RMSE = {results_bin['RMSE']:.4f}")
        print(f"   Pearson r = {results_bin['pearson_r']:.4f} (p < {results_bin['p_value']:.2e})")
        
        # Plot
        plot_fit_results(df_bin, y_pred_bin, 'binaries',
                        '../06_ANALISI/figures/fit_binaries.png')
        
        # Save results
        results_bin['parameters'] = {
            'alpha': alpha_bin, 'beta': beta_bin,
            'gamma': gamma_bin, 'a0': a0_bin
        }
        with open('../04_RISULTATI/fit_binaries.pkl', 'wb') as f:
            pickle.dump(results_bin, f)
    
    # ========================================================================
    # 5. COMBINED SUMMARY
    # ========================================================================
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    if df_exo is not None and df_bin is not None:
        N_total = len(df_exo) + len(df_bin)
        print(f"\nTotal systems analyzed: {N_total}")
        print(f"  - Exoplanets: {len(df_exo)}")
        print(f"  - Binaries: {len(df_bin)}")
        
        print(f"\nFit quality:")
        print(f"  - Exoplanets R² = {results_exo['R2']:.4f}")
        print(f"  - Binaries R² = {results_bin['R2']:.4f}")
        
        print(f"\nParameters:")
        print(f"  - α (exo) = {alpha_exo:.3f} ± {alpha_err:.3f}")
        print(f"  - α (bin) = {alpha_bin:.3f} ± {np.sqrt(pcov_bin[0,0]):.3f}")
        print(f"  - β (exo) = {beta_exo:.3f} ± {beta_err:.3f}")
        print(f"  - β (bin) = {beta_bin:.3f} ± {np.sqrt(pcov_bin[1,1]):.3f}")
    
    print("\n✅ Fitting complete!")
    print("Results saved to ../04_RISULTATI/")
    print("Figures saved to ../06_ANALISI/figures/")
    print()

if __name__ == "__main__":
    main()
