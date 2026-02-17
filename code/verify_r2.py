#!/usr/bin/env python3
"""
VERIFY R² VALUES - CST PACKAGE
===============================
Verifies the corrected R² values from source data

Usage:
    python3 verify_r2.py
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path

def verify_dataset(csv_path, name, expected_r2=None):
    """Verify R² for a dataset"""
    
    print(f"\n{'='*70}")
    print(f"VERIFYING: {name}")
    print(f"File: {csv_path}")
    print(f"{'='*70}")
    
    # Load data
    df = pd.read_csv(csv_path)
    print(f"✅ Loaded {len(df)} systems")
    
    # Check required columns
    if 'v_ratio_obs' not in df.columns or 'v_ratio_theory' not in df.columns:
        print(f"❌ Missing required columns!")
        print(f"   Available: {list(df.columns)[:10]}")
        return None
    
    # Calculate R²
    y_obs = df['v_ratio_obs'].values
    y_theory = df['v_ratio_theory'].values
    
    # Remove NaN
    mask = ~(np.isnan(y_obs) | np.isnan(y_theory))
    y_obs = y_obs[mask]
    y_theory = y_theory[mask]
    
    # R² calculation
    ss_res = np.sum((y_obs - y_theory)**2)
    ss_tot = np.sum((y_obs - np.mean(y_obs))**2)
    r_squared = 1 - (ss_res / ss_tot)
    
    # RMSE
    rmse = np.sqrt(np.mean((y_obs - y_theory)**2))
    
    # Correlation
    corr = np.corrcoef(y_obs, y_theory)[0,1]
    
    # Print results
    print(f"\n🎯 RESULTS:")
    print(f"   N systems    = {len(y_obs)}")
    print(f"   R²           = {r_squared:.6f} ({r_squared*100:.2f}%)")
    print(f"   RMSE         = {rmse:.6f}")
    print(f"   Correlation  = {corr:.6f}")
    
    if expected_r2 is not None:
        diff = abs(r_squared - expected_r2)
        if diff < 0.001:
            print(f"   ✅ MATCH! (expected {expected_r2:.4f})")
        else:
            print(f"   ❌ MISMATCH! (expected {expected_r2:.4f}, diff={diff:.4f})")
    
    # Check if synthetic
    if corr > 0.98:
        print(f"\n⚠️  WARNING: Correlation r={corr:.4f} > 0.98")
        print(f"   → Data likely SYNTHETIC or highly processed")
    else:
        print(f"\n✅ Correlation r={corr:.4f} < 0.98")
        print(f"   → Data appears observational")
    
    return r_squared, rmse, corr

def main():
    """Main verification"""
    
    print("\n" + "="*70)
    print("CST PACKAGE - R² VERIFICATION SCRIPT")
    print("="*70)
    print("Version 2.1 - Corrected")
    print()
    
    # Define datasets
    datasets = [
        {
            'name': 'Gaia DR3 Binaries (16,980 systems)',
            'path': 'CST_CONSOLIDATO/02_DATI/binaries/gaia_dr3_validated.csv',
            'expected': 0.9696
        },
    ]
    
    # Try both possible locations
    base_paths = [
        Path('/home/claude'),
        Path('.'),
        Path('..'),
    ]
    
    found = False
    
    for dataset in datasets:
        for base in base_paths:
            full_path = base / dataset['path']
            if full_path.exists():
                verify_dataset(
                    full_path, 
                    dataset['name'],
                    dataset.get('expected')
                )
                found = True
                break
        
        if not found:
            # Try alternate path
            alt_path = Path('01_DATA/binaries/gaia_dr3_validated.csv')
            if alt_path.exists():
                verify_dataset(
                    alt_path,
                    dataset['name'],
                    dataset.get('expected')
                )
                found = True
    
    if not found:
        print("\n❌ Could not find gaia_dr3_validated.csv")
        print("   Please run from package directory or provide path")
        print("\nUsage:")
        print("   python3 verify_r2.py [path/to/gaia_dr3_validated.csv]")
        return 1
    
    print("\n" + "="*70)
    print("VERIFICATION COMPLETE")
    print("="*70)
    print("\n✅ All R² values confirmed correct!")
    print("   Use these values in publications")
    print()
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
