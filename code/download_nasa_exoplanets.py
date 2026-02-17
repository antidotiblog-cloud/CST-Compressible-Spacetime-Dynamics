#!/usr/bin/env python3
"""
DOWNLOAD_NASA_EXOPLANETS.PY - NASA Exoplanet Archive Downloader
================================================================

Downloads complete exoplanet dataset from NASA Exoplanet Archive.

Retrieves:
- Confirmed planets (N ~ 5,500+)
- Stellar parameters (mass, age, metallicity, etc.)
- Orbital parameters (period, semi-major axis, eccentricity)

Output:
- CSV file with ~4,585 high-quality systems (after filtering)

API Documentation:
https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html

Author: Michele Vizzutti
Date: February 10, 2026
"""

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.constants import G, M_sun, R_sun, au
import requests
from io import StringIO
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# NASA TAP service URL
TAP_URL = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"

# Query parameters
QUERY = """
SELECT 
    pl_name,
    hostname,
    discoverymethod,
    disc_year,
    
    -- Planetary parameters
    pl_orbper,
    pl_orbsmax,
    pl_rade,
    pl_bmasse,
    pl_eqt,
    pl_orbeccen,
    
    -- Stellar parameters
    st_mass,
    st_rad,
    st_teff,
    st_logg,
    st_met,
    st_age,
    st_lum,
    st_dens,
    
    -- Errors (if available)
    st_masserr1,
    st_masserr2,
    st_ageerr1,
    st_ageerr2,
    
    -- System info
    sy_dist,
    sy_vmag

FROM 
    ps

WHERE 
    -- Basic requirements
    pl_orbper IS NOT NULL
    AND pl_orbsmax IS NOT NULL
    AND st_mass IS NOT NULL
    AND st_age IS NOT NULL
    
    -- Quality filters
    AND st_mass > 0
    AND st_age > 0
    AND st_age < 13.8
    AND pl_orbsmax > 0
    
    -- Exclude brown dwarfs
    AND pl_rade < 30

ORDER BY 
    pl_name
"""

# Output filepath
OUTPUT_FILE = '../02_DATI/exoplanets/exoplanet_data_full.csv'

# ============================================================================
# DOWNLOAD FUNCTION
# ============================================================================

def download_exoplanets(query=QUERY, url=TAP_URL, output_format='csv'):
    """
    Download exoplanet data from NASA Archive via TAP.
    
    Parameters:
    -----------
    query : str
        ADQL query string
    url : str
        TAP service URL
    output_format : str
        Output format ('csv', 'votable', 'json')
    
    Returns:
    --------
    df : pandas.DataFrame
        Downloaded data
    """
    print("Downloading from NASA Exoplanet Archive...")
    print(f"URL: {url}")
    print()
    
    # Prepare request
    params = {
        'REQUEST': 'doQuery',
        'LANG': 'ADQL',
        'FORMAT': output_format,
        'QUERY': query
    }
    
    try:
        # Send request
        response = requests.get(url, params=params, timeout=120)
        response.raise_for_status()
        
        # Parse response
        data = response.text
        df = pd.read_csv(StringIO(data))
        
        print(f"✅ Downloaded {len(df)} systems")
        print()
        
        return df
    
    except requests.exceptions.Timeout:
        print("❌ Request timed out (>120s)")
        print("   Try again or reduce query size")
        return None
    
    except requests.exceptions.RequestException as e:
        print(f"❌ Download failed: {e}")
        return None
    
    except Exception as e:
        print(f"❌ Parsing failed: {e}")
        return None

# ============================================================================
# DERIVED PARAMETERS
# ============================================================================

def calculate_cosmology(df):
    """
    Calculate cosmological parameters from stellar age.
    """
    print("Calculating cosmology...")
    
    # Age -> redshift
    age_universe = 13.8  # Gyr
    t_formation = age_universe - df['st_age']
    t_formation = np.maximum(t_formation, 0.01)
    df['z_formation'] = (age_universe / t_formation)**(2/3) - 1
    df['z_formation'] = np.maximum(df['z_formation'], 0)
    
    # H(z)/H₀
    Omega_m = 0.315
    Omega_L = 0.685
    z = df['z_formation']
    df['H_ratio'] = np.sqrt(Omega_m * (1+z)**3 + Omega_L)
    
    # Transition function
    z_trans = 30
    H_ratio = df['H_ratio']
    S = 1.0 / (1.0 + (z/z_trans)**3)
    df['f_z'] = H_ratio * S
    
    print(f"   z range: {df['z_formation'].min():.3f} - {df['z_formation'].max():.3f}")
    print(f"   H/H₀ range: {df['H_ratio'].min():.3f} - {df['H_ratio'].max():.3f}")
    print()

def calculate_velocities(df):
    """
    Calculate Keplerian and observed velocities.
    """
    print("Calculating velocities...")
    
    # Constants (SI)
    G_val = 6.674e-11  # m³/(kg·s²)
    M_sun_val = 1.989e30  # kg
    AU_val = 1.496e11  # m
    
    # Keplerian velocity
    M_kg = df['st_mass'] * M_sun_val
    a_m = df['pl_orbsmax'] * AU_val
    df['v_kepler'] = np.sqrt(G_val * M_kg / a_m)
    
    # Observed velocity (from period)
    P_seconds = df['pl_orbper'] * 86400
    df['v_observed'] = 2 * np.pi * a_m / P_seconds
    
    # Velocity ratio
    df['v_ratio'] = df['v_observed'] / df['v_kepler']
    
    # Enhancement percentage
    df['enhancement_pct'] = (df['v_ratio'] - 1) * 100
    
    print(f"   v_ratio range: {df['v_ratio'].min():.4f} - {df['v_ratio'].max():.4f}")
    print(f"   Enhancement: {df['enhancement_pct'].min():.2f}% - {df['enhancement_pct'].max():.2f}%")
    print()

def calculate_derived_all(df):
    """
    Calculate all derived parameters.
    """
    calculate_cosmology(df)
    calculate_velocities(df)
    return df

# ============================================================================
# QUALITY FILTERING
# ============================================================================

def apply_quality_filters(df):
    """
    Apply quality filters to dataset.
    
    Filters:
    --------
    1. Physical consistency (v_ratio > 0.5, < 2.0)
    2. Age range (0 < age < 13.8 Gyr, with buffer)
    3. Stellar mass range (0.08 < M < 3.0 M_☉)
    4. Orbital validity (a > 0, P > 0)
    5. Remove extreme outliers (|residual| > 5σ)
    """
    print("Applying quality filters...")
    N_initial = len(df)
    
    # 1. Physical v_ratio
    mask = (df['v_ratio'] > 0.5) & (df['v_ratio'] < 2.0)
    df = df[mask].copy()
    print(f"   Physical v_ratio: {len(df)}/{N_initial} retained")
    
    # 2. Age range
    mask = (df['st_age'] > 0) & (df['st_age'] < 13.5)
    df = df[mask].copy()
    print(f"   Age range: {len(df)}/{N_initial} retained")
    
    # 3. Stellar mass
    mask = (df['st_mass'] > 0.08) & (df['st_mass'] < 3.0)
    df = df[mask].copy()
    print(f"   Stellar mass: {len(df)}/{N_initial} retained")
    
    # 4. Orbital validity
    mask = (df['pl_orbsmax'] > 0) & (df['pl_orbper'] > 0)
    df = df[mask].copy()
    print(f"   Orbital validity: {len(df)}/{N_initial} retained")
    
    # 5. Remove extreme outliers (simple 5σ filter)
    residuals = df['v_ratio'] - 1.0
    threshold = 5 * np.std(residuals)
    mask = np.abs(residuals) < threshold
    df = df[mask].copy()
    print(f"   Outlier removal: {len(df)}/{N_initial} retained")
    
    print(f"\n   Final: {len(df)} high-quality systems")
    print()
    
    return df

# ============================================================================
# DATA SUMMARY
# ============================================================================

def print_data_summary(df):
    """Print comprehensive data summary."""
    print("="*70)
    print("DATA SUMMARY")
    print("="*70)
    print()
    
    print(f"Total systems: {len(df)}")
    print()
    
    print("Stellar parameters:")
    print(f"  Mass: {df['st_mass'].min():.2f} - {df['st_mass'].max():.2f} M_☉")
    print(f"        (median: {df['st_mass'].median():.2f}, mean: {df['st_mass'].mean():.2f})")
    print(f"  Age:  {df['st_age'].min():.2f} - {df['st_age'].max():.2f} Gyr")
    print(f"        (median: {df['st_age'].median():.2f}, mean: {df['st_age'].mean():.2f})")
    
    if 'st_met' in df.columns:
        print(f"  [Fe/H]: {df['st_met'].min():.2f} - {df['st_met'].max():.2f}")
        print(f"          (median: {df['st_met'].median():.2f})")
    
    print()
    
    print("Orbital parameters:")
    print(f"  Period: {df['pl_orbper'].min():.2f} - {df['pl_orbper'].max():.2f} days")
    print(f"          (median: {df['pl_orbper'].median():.2f})")
    print(f"  Semi-major axis: {df['pl_orbsmax'].min():.4f} - {df['pl_orbsmax'].max():.2f} AU")
    print(f"                   (median: {df['pl_orbsmax'].median():.2f})")
    print()
    
    print("Cosmology:")
    print(f"  Redshift: {df['z_formation'].min():.3f} - {df['z_formation'].max():.3f}")
    print(f"            (median: {df['z_formation'].median():.3f})")
    print(f"  H/H₀: {df['H_ratio'].min():.3f} - {df['H_ratio'].max():.3f}")
    print(f"        (median: {df['H_ratio'].median():.3f})")
    print()
    
    print("Observables:")
    print(f"  v_ratio: {df['v_ratio'].min():.4f} - {df['v_ratio'].max():.4f}")
    print(f"           (median: {df['v_ratio'].median():.4f})")
    print(f"  Enhancement: {df['enhancement_pct'].min():.2f}% - {df['enhancement_pct'].max():.2f}%")
    print(f"               (median: {df['enhancement_pct'].median():.2f}%)")
    print()
    
    print("Discovery methods:")
    if 'discoverymethod' in df.columns:
        method_counts = df['discoverymethod'].value_counts()
        for method, count in method_counts.head(5).items():
            pct = 100 * count / len(df)
            print(f"  {method}: {count} ({pct:.1f}%)")
    print()

# ============================================================================
# MAIN
# ============================================================================

def main():
    """
    Main execution function.
    """
    print("="*70)
    print("NASA EXOPLANET ARCHIVE DOWNLOADER")
    print("="*70)
    print()
    
    # 1. Download
    df = download_exoplanets()
    
    if df is None:
        print("❌ Download failed. Exiting.")
        return
    
    # 2. Calculate derived parameters
    df = calculate_derived_all(df)
    
    # 3. Apply quality filters
    df = apply_quality_filters(df)
    
    # 4. Summary
    print_data_summary(df)
    
    # 5. Save
    import os
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    df.to_csv(OUTPUT_FILE, index=False)
    
    print("="*70)
    print(f"✅ COMPLETE")
    print("="*70)
    print(f"\nDataset saved to: {OUTPUT_FILE}")
    print(f"Systems: {len(df)}")
    print(f"Columns: {len(df.columns)}")
    print()
    print("Ready for fitting with: python fit_master.py")
    print()

if __name__ == "__main__":
    main()
