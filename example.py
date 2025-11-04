#!/usr/bin/env python3
"""
Example: Basic DOP Analysis
============================

This script demonstrates basic usage of the LEO-PNT DOP analyzer
for a satellite constellation.
"""

import sys
import os
from datetime import datetime

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from dop_ana import DOPAnalyzer
import tkinter as tk
from tkinter import filedialog
from vis import (create_dop_contour_map, create_comparison_plot,
                          print_statistics, save_all_plots)


def main():
    """Run basic DOP analysis."""
    
    print("="*80)
    print("LEO-PNT DOP ANALYSIS - BASIC EXAMPLE")
    print("="*80)
    
    # 1. Load TLE data
    print("\n1. Loading TLE data...")
    
    # Specify the TLE file path directly
    tle_file = 'STARLINK_DTC_TLES.txt'
    
    if not os.path.exists(tle_file):
        print(f"Error: Selected file not found at {tle_file}")
        return
    
    with open(tle_file, 'r') as f:
        tle_data = f.read()
    
    print(f"✓ Loaded TLE data from {tle_file}")
    # 2. Initialize analyzer
    print("\n2. Initializing DOP Analyzer...")
    mask_angles = [5, 20, 30]  # degrees
    analyzer = DOPAnalyzer(tle_data, mask_angles=mask_angles)
    print(f"✓ Initialized with {len(analyzer.satellites)} satellites")
    print(f"✓ Mask angles: {mask_angles}°")
    
    # 3. Set analysis parameters
    print("\n3. Setting analysis parameters...")
    analysis_time = datetime(2025, 1, 1, 12, 0, 0)  # Jan 1, 2025, 12:00 UTC
    resolution = 15  # degrees (15° = fast, 5° = detailed but slower)
    
    print(f"✓ Analysis time: {analysis_time.strftime('%Y-%m-%d %H:%M:%S UTC')}")
    print(f"✓ Grid resolution: {resolution}°")
    
    # 4. Compute DOP for each mask angle
    print("\n4. Computing DOP values...")
    results = {}
    
    for mask_angle in mask_angles:
        print(f"\n  Analyzing mask angle: {mask_angle}°")
        print("  " + "─"*60)
        
        lats, lons, dop_values = analyzer.generate_global_dop_grid(
            analysis_time, mask_angle, resolution
        )
        
        results[mask_angle] = {
            'lats': lats,
            'lons': lons,
            'dop_values': dop_values
        }
    
    print("\n✓ DOP computation complete!")
    
    # 5. Print statistics
    print_statistics(results)
    
    # 6. Create visualizations
    print("\n5. Creating visualizations...")
    
    # Example: Create HDOP map for 20° mask
    print("\n  Creating HDOP contour map (20° mask)...")
    fig = create_dop_contour_map(
        results[20]['lats'],
        results[20]['lons'],
        results[20]['dop_values'],
        dop_type='HDOP',
        mask_angle=20,
        title_suffix=f'({len(analyzer.satellites)} satellites)'
    )
    import matplotlib.pyplot as plt
    plt.show()
    
    # Create comparison plot
    print("\n  Creating comparison plot...")
    fig_comp = create_comparison_plot(results)
    plt.show()
    
    # 7. Save all plots
    print("\n6. Saving plots...")
    output_dir = os.path.join('..', 'plots')
    save_all_plots(results, output_dir=output_dir, constellation_name='example')
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)
    print(f"• Analyzed {len(analyzer.satellites)} satellites")
    print(f"• Mask angles: {mask_angles}°")
    print(f"• Time: {analysis_time.strftime('%Y-%m-%d %H:%M:%S UTC')}")
    print(f"• Plots saved to: {output_dir}/")
    print("="*80)


if __name__ == '__main__':
    main()
