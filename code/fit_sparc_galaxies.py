#!/usr/bin/env python3
"""
FIT_SPARC_GALAXIES.PY - α_cosmo Validation from Galaxy Rotation Curves
========================================================================

Tests whether α_cosmo ≈ 0.05-0.10 explains galaxy rotation curves.

SPARC Dataset:
--------------
- 175 galaxies with high-quality rotation curves
- Photometric data (3.6 μm Spitzer)
- Gas kinematics (HI 21cm line)
- Stellar mass profiles

Reference:
----------
Lelli et al. (2016), "SPARC: Mass Models for 175 Disk Galaxies"
AJ, 152, 157
http://astroweb.cwru.edu/SPARC/

Theory Test:
------------
Extended structures: G_eff = G_N × [1 + α_cosmo × f(z)]

α_cosmo << α_compact because:
1. Distributed mass (not single coherent M)
2. Geometric averaging over volume
3. Partial destructive interference

Prediction: α_cosmo ≈ 0.05-0.10  (α_compact/4)

Author: Michele Vizzutti
Date: February 10, 2026
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats
import requests
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# SPARC data URL (master table)
SPARC_URL = "http://astroweb.cwru.edu/SPARC/MasterTable.txt"

# Individual galaxy rotation curves
SPARC_RC_URL = "http://astroweb.cwru.edu/SPARC/{galaxy}_rotmod.dat"

OUTPUT_DIR = '../02_DATI/sparc/'

# CST parameters
ALPHA_COMPACT = 0.279  # From solar system fit
BETA = 0.685
Z_TRANS = 30.0

# ============================================================================
# DOWNLOAD FUNCTIONS
# ============================================================================

def download_sparc_master():
    """
    Download SPARC master table.
    
    Contains:
    - Galaxy names
    - Distances
    - Inclinations
    - Stellar masses
    - etc.
    
    Returns:
    --------
    df : pandas.DataFrame
    """
    print("Downloading SPARC master table...")
    
    try:
        # Download
        response = requests.get(SPARC_URL, timeout=30)
        response.raise_for_status()
        
        # Parse (fixed-width format)
        lines = response.text.strip().split('\n')
        
        # Skip header (first few lines)
        data_lines = [l for l in lines if not l.startswith('#') and len(l) > 10]
        
        # Parse columns (see SPARC documentation for format)
        galaxies = []
        for line in data_lines:
            parts = line.split()
            if len(parts) >= 12:
                galaxies.append({
                    'Galaxy': parts[0],
                    'Type': parts[1],
                    'Dist': float(parts[2]),  # Mpc
                    'e_Dist': float(parts[3]),
                    'Inc': float(parts[4]),  # deg
                    'e_Inc': float(parts[5]),
                    'L36': float(parts[6]),  # 10^9 L_sun
                    'e_L36': float(parts[7]),
                    'Vflat': float(parts[8]),  # km/s
                    'e_Vflat': float(parts[9]),
                    'Q': float(parts[10]),  # Quality flag
                    'Ref': parts[11] if len(parts) > 11 else ''
                })
        
        df = pd.DataFrame(galaxies)
        
        print(f"✅ Downloaded {len(df)} galaxies")
        print()
        
        return df
    
    except Exception as e:
        print(f"❌ Download failed: {e}")
        return None

def download_rotation_curve(galaxy_name):
    """
    Download individual rotation curve for galaxy.
    
    Parameters:
    -----------
    galaxy_name : str
        Galaxy identifier from master table
    
    Returns:
    --------
    df : pandas.DataFrame
        Columns: Rad (kpc), Vobs (km/s), errV, Vgas, Vdisk, Vbul
    """
    url = SPARC_RC_URL.format(galaxy=galaxy_name)
    
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        
        # Parse
        lines = response.text.strip().split('\n')
        data_lines = [l for l in lines if not l.startswith('#') and len(l) > 10]
        
        data = []
        for line in data_lines:
            parts = line.split()
            if len(parts) >= 6:
                data.append({
                    'Rad': float(parts[0]),  # kpc
                    'Vobs': float(parts[1]),  # km/s
                    'errV': float(parts[2]),
                    'Vgas': float(parts[3]),
                    'Vdisk': float(parts[4]),
                    'Vbul': float(parts[5]) if len(parts) > 5 else 0.0
                })
        
        return pd.DataFrame(data)
    
    except:
        return None

# ============================================================================
# COSMOLOGY FOR GALAXIES
# ============================================================================

def estimate_galaxy_age(galaxy_type):
    """
    Estimate typical age of galaxy by morphological type.
    
    Rough approximation:
    - Sa, Sb: old, z_form ~ 2-4
    - Sc, Sd: younger, z_form ~ 1-2
    - Irr, BCD: very young, z_form ~ 0.5-1
    
    Returns:
    --------
    age_Gyr : float
    """
    if galaxy_type in ['Sa', 'Sab', 'Sb']:
        return 10.0  # Old, z_form ~ 3
    elif galaxy_type in ['Sbc', 'Sc', 'Scd']:
        return 8.0  # Intermediate, z_form ~ 1.5
    elif galaxy_type in ['Sd', 'Sm', 'Im']:
        return 5.0  # Young, z_form ~ 0.8
    else:
        return 7.0  # Default

def calculate_cosmology_galaxy(df):
    """
    Calculate cosmology for each galaxy.
    """
    # Estimate ages
    df['age_Gyr'] = df['Type'].apply(estimate_galaxy_age)
    
    # z_formation
    age_universe = 13.8
    t_formation = age_universe - df['age_Gyr']
    t_formation = np.maximum(t_formation, 0.01)
    df['z_formation'] = (age_universe / t_formation)**(2/3) - 1
    df['z_formation'] = np.maximum(df['z_formation'], 0)
    
    # H(z)/H₀
    Omega_m = 0.315
    Omega_L = 0.685
    z = df['z_formation']
    df['H_ratio'] = np.sqrt(Omega_m * (1+z)**3 + Omega_L)
    
    # Transition function
    H_ratio = df['H_ratio']
    S = 1.0 / (1.0 + (z/Z_TRANS)**3)
    df['f_z'] = H_ratio * S

# ============================================================================
# FITTING MODEL
# ============================================================================

def v_model_CST(R_kpc, M_enc_Msun, alpha_cosmo, f_z):
    """
    Rotation velocity with G_eff for extended structures.
    
    V(R) = sqrt(G_eff × M_enc / R)
    
    where G_eff = G_N × [1 + α_cosmo × f(z)]
    
    Parameters:
    -----------
    R_kpc : float or array
        Radius in kpc
    M_enc_Msun : float or array
        Enclosed mass in M_☉
    alpha_cosmo : float
        Cosmological coupling (to fit)
    f_z : float
        Transition function value
    
    Returns:
    --------
    V : float or array
        Velocity in km/s
    """
    # Constants
    G_N = 4.302e-6  # kpc (km/s)² / M_☉
    
    # G_eff
    G_eff = G_N * (1 + alpha_cosmo * f_z)
    
    # Velocity
    V = np.sqrt(G_eff * M_enc_Msun / R_kpc)
    
    return V

def fit_galaxy_rotation_curve(rc_df, f_z):
    """
    Fit α_cosmo for a single galaxy's rotation curve.
    
    Assumes:
    - Total M_enc = M_baryon + M_DM
    - M_baryon from photometry (Vdisk, Vgas, Vbul)
    - M_DM = dark matter (fit as free parameter)
    
    Parameters:
    -----------
    rc_df : pandas.DataFrame
        Rotation curve data
    f_z : float
        Transition function for this galaxy
    
    Returns:
    --------
    alpha_cosmo : float
        Fitted coupling
    M_DM_profile : array
        Dark matter mass profile
    chi2 : float
        Goodness of fit
    """
    R = rc_df['Rad'].values  # kpc
    V_obs = rc_df['Vobs'].values  # km/s
    V_err = rc_df['errV'].values
    V_gas = rc_df['Vgas'].values
    V_disk = rc_df['Vdisk'].values
    V_bul = rc_df['Vbul'].values
    
    # Baryonic component (known from observations)
    V_baryon = np.sqrt(V_gas**2 + V_disk**2 + V_bul**2)
    
    # Enclosed baryonic mass (approximate from V_baryon)
    G_N = 4.302e-6  # kpc (km/s)² / M_☉
    M_baryon = V_baryon**2 * R / G_N
    
    # Fit function: Total velocity with DM + CST modification
    def model(R, alpha_cosmo, M_DM_norm):
        """
        V_tot = sqrt(V_baryon² + V_DM² × G_eff/G_N)
        
        where V_DM = sqrt(G_N × M_DM / R)
              G_eff = G_N × (1 + α_cosmo × f_z)
        """
        # DM profile (NFW-like, simplified)
        M_DM = M_DM_norm * np.log(1 + R / 10)  # Rough approximation
        
        # DM velocity with G_eff
        G_eff = G_N * (1 + alpha_cosmo * f_z)
        V_DM = np.sqrt(G_eff * M_DM / R)
        
        # Total
        V_tot = np.sqrt(V_baryon**2 + V_DM**2)
        
        return V_tot
    
    try:
        popt, pcov = curve_fit(
            model, R, V_obs,
            p0=[0.07, 1e10],  # Initial guess
            bounds=([0, 0], [0.5, 1e12]),
            sigma=V_err,
            absolute_sigma=True,
            maxfev=5000
        )
        
        alpha_cosmo, M_DM_norm = popt
        
        # Chi-squared
        V_model = model(R, alpha_cosmo, M_DM_norm)
        chi2 = np.sum(((V_obs - V_model) / V_err)**2)
        dof = len(R) - 2
        chi2_red = chi2 / dof
        
        return alpha_cosmo, M_DM_norm, chi2_red
    
    except:
        return None, None, None

# ============================================================================
# BATCH FITTING
# ============================================================================

def fit_all_galaxies(master_df, N_max=50):
    """
    Fit α_cosmo for multiple galaxies.
    
    Parameters:
    -----------
    master_df : pandas.DataFrame
        SPARC master table
    N_max : int
        Maximum number of galaxies to fit (for speed)
    
    Returns:
    --------
    results : pandas.DataFrame
        Fitted parameters per galaxy
    """
    print(f"\nFitting rotation curves (max {N_max} galaxies)...")
    print()
    
    results = []
    
    for i, row in master_df.head(N_max).iterrows():
        galaxy = row['Galaxy']
        f_z = row['f_z']
        
        print(f"[{i+1}/{N_max}] {galaxy}...", end=' ')
        
        # Download RC
        rc_df = download_rotation_curve(galaxy)
        
        if rc_df is None or len(rc_df) < 5:
            print("❌ No data")
            continue
        
        # Fit
        alpha, M_DM, chi2_red = fit_galaxy_rotation_curve(rc_df, f_z)
        
        if alpha is not None:
            print(f"α = {alpha:.3f}, χ²/dof = {chi2_red:.2f}")
            
            results.append({
                'Galaxy': galaxy,
                'Type': row['Type'],
                'Dist_Mpc': row['Dist'],
                'z_formation': row['z_formation'],
                'f_z': f_z,
                'alpha_cosmo': alpha,
                'M_DM_norm': M_DM,
                'chi2_red': chi2_red,
                'N_points': len(rc_df)
            })
        else:
            print("❌ Fit failed")
    
    return pd.DataFrame(results)

# ============================================================================
# ANALYSIS
# ============================================================================

def analyze_alpha_cosmo(results_df):
    """
    Analyze distribution of fitted α_cosmo.
    """
    print("\n" + "="*70)
    print("α_cosmo ANALYSIS")
    print("="*70)
    print()
    
    alpha = results_df['alpha_cosmo'].values
    
    # Statistics
    print(f"N galaxies: {len(alpha)}")
    print(f"Mean α_cosmo: {alpha.mean():.3f}")
    print(f"Median α_cosmo: {np.median(alpha):.3f}")
    print(f"Std dev: {alpha.std():.3f}")
    print(f"Range: {alpha.min():.3f} - {alpha.max():.3f}")
    print()
    
    # Compare to prediction
    alpha_predicted = ALPHA_COMPACT / 4  # Theory: α_cosmo ~ α/4
    print(f"Predicted α_cosmo: {alpha_predicted:.3f} (α_compact/4)")
    print()
    
    # Statistical test
    t_stat, p_value = stats.ttest_1samp(alpha, alpha_predicted)
    print(f"T-test vs prediction:")
    print(f"  t-statistic: {t_stat:.2f}")
    print(f"  p-value: {p_value:.4f}")
    
    if p_value > 0.05:
        print(f"  ✅ Consistent with prediction (p > 0.05)")
    else:
        print(f"  ⚠️  Deviates from prediction (p < 0.05)")
    print()
    
    # Quality filter (chi2_red < 3)
    good = results_df[results_df['chi2_red'] < 3]['alpha_cosmo']
    if len(good) > 0:
        print(f"Good fits only (χ²/dof < 3): N = {len(good)}")
        print(f"  Mean α_cosmo: {good.mean():.3f}")
        print(f"  Median: {good.median():.3f}")
        print()

# ============================================================================
# PLOTTING
# ============================================================================

def plot_alpha_cosmo_distribution(results_df, output_file=None):
    """
    Plot distribution of α_cosmo.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    alpha = results_df['alpha_cosmo'].values
    
    # 1. Histogram
    ax = axes[0]
    ax.hist(alpha, bins=20, alpha=0.7, edgecolor='black', color='blue')
    
    # Mean and prediction
    alpha_mean = alpha.mean()
    alpha_pred = ALPHA_COMPACT / 4
    
    ax.axvline(alpha_mean, color='red', linestyle='--', lw=2,
               label=f'Mean = {alpha_mean:.3f}')
    ax.axvline(alpha_pred, color='green', linestyle='--', lw=2,
               label=f'Predicted = {alpha_pred:.3f}')
    
    ax.set_xlabel('α_cosmo', fontsize=12)
    ax.set_ylabel('Number of Galaxies', fontsize=12)
    ax.set_title('Distribution of α_cosmo from SPARC', fontsize=14, weight='bold')
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3)
    
    # 2. α_cosmo vs galaxy type
    ax = axes[1]
    
    types_unique = results_df['Type'].unique()
    for gtype in types_unique:
        mask = results_df['Type'] == gtype
        alpha_sub = results_df[mask]['alpha_cosmo']
        ax.scatter([gtype]*len(alpha_sub), alpha_sub, alpha=0.6, s=50, label=gtype)
    
    ax.axhline(alpha_mean, color='red', linestyle='--', lw=2, label='Mean')
    ax.axhline(alpha_pred, color='green', linestyle='--', lw=2, label='Predicted')
    
    ax.set_xlabel('Galaxy Type', fontsize=12)
    ax.set_ylabel('α_cosmo', fontsize=12)
    ax.set_title('α_cosmo by Morphology', fontsize=14, weight='bold')
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✅ Figure saved to {output_file}")
    else:
        plt.show()
    
    plt.close()

# ============================================================================
# MAIN
# ============================================================================

def main():
    """
    Main execution.
    """
    print("="*70)
    print("SPARC GALAXIES - α_cosmo VALIDATION")
    print("="*70)
    print()
    
    # 1. Download master table
    master_df = download_sparc_master()
    
    if master_df is None:
        print("❌ Failed. Exiting.")
        return
    
    # 2. Calculate cosmology
    calculate_cosmology_galaxy(master_df)
    
    # 3. Fit rotation curves
    results_df = fit_all_galaxies(master_df, N_max=50)
    
    if len(results_df) == 0:
        print("❌ No successful fits. Exiting.")
        return
    
    # 4. Analyze
    analyze_alpha_cosmo(results_df)
    
    # 5. Plot
    plot_alpha_cosmo_distribution(results_df, '../06_ANALISI/figures/alpha_cosmo_sparc.png')
    
    # 6. Save
    results_df.to_csv('../04_RISULTATI/sparc_alpha_cosmo.csv', index=False)
    print("✅ Results saved to ../04_RISULTATI/sparc_alpha_cosmo.csv")
    print()

if __name__ == "__main__":
    main()
