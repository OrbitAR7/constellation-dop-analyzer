"""
Visualization utilities for DOP analysis.

Provides functions to create:
- Static contour maps using Matplotlib and Cartopy
- Summary statistics and comparison plots
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import Rectangle

import os
from datetime import datetime

def make_output_dir(constellation_name: str):
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out = os.path.abspath(os.path.join(".", "dop_plots", f"{constellation_name}_{ts}"))
    os.makedirs(out, exist_ok=True)
    print(f"Saving plots to: {out}")
    return out

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False
    print("Warning: Cartopy not available. Install with: pip install cartopy")


def create_dop_contour_map(lats, lons, dop_values, dop_type='HDOP', 
                           mask_angle=30, title_suffix='', vmax=10):
    """
    Create a global contour map of DOP values.
    
    Args:
        lats (list): Latitude values
        lons (list): Longitude values  
        dop_values (list): List of DOP dictionaries
        dop_type (str): Type of DOP to plot ('HDOP', 'VDOP', 'PDOP', 'GDOP', 'TDOP')
        mask_angle (float): Elevation mask angle used
        title_suffix (str): Additional text for title
        vmax (float): Maximum value for color scale
        
    Returns:
        matplotlib.figure.Figure: The created figure
    """
    if not HAS_CARTOPY:
        raise ImportError("Cartopy is required for map visualization")
    
    # Extract values for specified DOP type
    values = np.array([dop[dop_type] for dop in dop_values])
    
    # Cap extreme values for better visualization
    values = np.clip(values, 0, 50)
    
    # Reshape into grid
    lats_unique = sorted(set(lats))
    lons_unique = sorted(set(lons))
    
    grid = np.full((len(lats_unique), len(lons_unique)), np.nan)
    
    for i, (lat, lon, val) in enumerate(zip(lats, lons, values)):
        lat_idx = lats_unique.index(lat)
        lon_idx = lons_unique.index(lon)
        grid[lat_idx, lon_idx] = val
    
    # Create figure with map projection
    fig = plt.figure(figsize=(14, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    # Add map features
    ax.coastlines(resolution='110m', linewidth=0.5, color='black')
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, edgecolor='gray')
    ax.add_feature(cfeature.LAND, facecolor='lightgray', alpha=0.3)
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue', alpha=0.2)
    
    # Create contour plot
    levels = np.linspace(0, vmax, 21)
    
    contourf = ax.contourf(lons_unique, lats_unique, grid,
                           levels=levels,
                           cmap='RdYlGn_r',
                           extend='max',
                           transform=ccrs.PlateCarree())
    
    contour = ax.contour(lons_unique, lats_unique, grid,
                        levels=levels[::2],
                        colors='black',
                        linewidths=0.3,
                        alpha=0.4,
                        transform=ccrs.PlateCarree())
    
    # Add colorbar
    cbar = plt.colorbar(contourf, ax=ax, orientation='vertical',
                       pad=0.05, shrink=0.8)
    cbar.set_label(f'{dop_type} Value', fontsize=11, weight='bold')
    
    # Title and grid
    ax.set_title(f'Global {dop_type} Analysis - Mask Angle: {mask_angle}° {title_suffix}',
                fontsize=13, weight='bold', pad=10)
    ax.gridlines(draw_labels=True, linewidth=0.5, alpha=0.5, linestyle='--')
    
    # Add DOP interpretation legend
    legend_text = _get_dop_interpretation(dop_type)
    ax.text(0.02, 0.02, legend_text,
            transform=ax.transAxes,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
            fontsize=9, verticalalignment='bottom')
    
    plt.tight_layout()
    return fig


def _get_dop_interpretation(dop_type):
    """Get interpretation guide for DOP values."""
    interpretations = {
        'HDOP': 'Horizontal DOP:\n< 1: Ideal\n1-2: Excellent\n2-5: Good\n5-10: Moderate\n> 10: Poor',
        'VDOP': 'Vertical DOP:\n< 1: Ideal\n1-2: Excellent\n2-5: Good\n5-10: Moderate\n> 10: Poor',
        'PDOP': 'Position DOP:\n< 2: Ideal\n2-3: Excellent\n3-6: Good\n6-10: Moderate\n> 10: Poor',
        'GDOP': 'Geometric DOP:\n< 2: Ideal\n2-4: Excellent\n4-8: Good\n8-12: Moderate\n> 12: Poor',
        'TDOP': 'Time DOP:\n< 1: Ideal\n1-2: Excellent\n2-5: Good\n5-10: Moderate\n> 10: Poor'
    }
    return interpretations.get(dop_type, 'DOP Values')



def create_outage_map(lats, lons, dop_values, title="Availability (valid 3D fix)"):
    """Create a map where value=1 if PDOP is finite (valid fix) else 0."""

    valid = np.array([1.0 if np.isfinite(d['PDOP']) else 0.0 for d in dop_values])
    # Build grid
    lats_unique = sorted(set(lats))
    lons_unique = sorted(set(lons))
    grid = np.full((len(lats_unique), len(lons_unique)), np.nan)
    for lat, lon, v in zip(lats, lons, valid):
        i = lats_unique.index(lat)
        j = lons_unique.index(lon)
        grid[i, j] = v

    # Try cartopy map; fallback to plain imshow
    fig = plt.figure(figsize=(14, 8))
    try:
        import cartopy.crs as ccrs
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines(linewidth=0.5)
        cs = ax.contourf(lons_unique, lats_unique, grid, levels=[-0.1,0.5,1.1], extend='neither', transform=ccrs.PlateCarree())
        ax.set_global()
        ax.set_title(title)
    except Exception as e:
        ax = plt.gca()
        im = ax.imshow(grid, origin='lower', extent=[min(lons_unique), max(lons_unique), min(lats_unique), max(lats_unique)], aspect='auto')
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title(title + " (no cartopy)")
    return fig
def create_comparison_plot(results, dop_types=['HDOP', 'VDOP', 'PDOP', 'GDOP']):
    """
    Create a comparison plot showing DOP statistics across mask angles.
    
    Args:
        results (dict): Results dictionary with mask angles as keys
        dop_types (list): List of DOP types to compare
        
    Returns:
        matplotlib.figure.Figure: The created figure
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('DOP Analysis: Mask Angle Comparison', 
                fontsize=15, weight='bold')
    
    axes_flat = axes.flatten()
    
    for idx, dop_type in enumerate(dop_types):
        ax = axes_flat[idx]
        
        mask_angles = sorted(results.keys())
        means = []
        stds = []
        
        for mask_angle in mask_angles:
            dop_vals = results[mask_angle]['dop_values']
            values = [d[dop_type] for d in dop_vals if np.isfinite(d[dop_type])]
            
            if values:
                means.append(np.mean(values))
                stds.append(np.std(values))
            else:
                means.append(0)
                stds.append(0)
        
        x_pos = np.arange(len(mask_angles))
        colors_list = ['green', 'orange', 'red'][:len(mask_angles)]
        
        bars = ax.bar(x_pos, means, yerr=stds, capsize=5,
                     color=colors_list, alpha=0.7, edgecolor='black')
        
        ax.set_xlabel('Mask Angle (degrees)', fontsize=10, weight='bold')
        ax.set_ylabel(f'{dop_type} Value', fontsize=10, weight='bold')
        ax.set_title(f'{dop_type} vs Mask Angle', fontsize=11, weight='bold')
        ax.set_xticks(x_pos)
        ax.set_xticklabels([f'{angle}°' for angle in mask_angles])
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add value labels
        for bar, mean in zip(bars, means):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                   f'{mean:.2f}', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    return fig


def print_statistics(results):
    """
    Print summary statistics for DOP analysis.
    
    Args:
        results (dict): Results dictionary with mask angles as keys
    """
    print("\n" + "="*80)
    print("DOP ANALYSIS STATISTICS")
    print("="*80)
    
    for mask_angle in sorted(results.keys()):
        print(f"\n{'─'*80}")
        print(f"Mask Angle: {mask_angle}°")
        print(f"{'─'*80}")
        
        dop_vals = results[mask_angle]['dop_values']
        
        for dop_type in ['HDOP', 'VDOP', 'PDOP', 'GDOP', 'TDOP']:
            values = [d[dop_type] for d in dop_vals if np.isfinite(d[dop_type])]
            
            if values:
                print(f"\n{dop_type}:")
                print(f"  Mean:   {np.mean(values):6.2f}")
                print(f"  Median: {np.median(values):6.2f}")
                print(f"  Min:    {np.min(values):6.2f}")
                print(f"  Max:    {np.max(values):6.2f}")
                print(f"  Std:    {np.std(values):6.2f}")
                print(f"  P95:    {np.percentile(values, 95):6.2f}")
                print(f"  Coverage: {100*len(values)/len(dop_vals):5.1f}%")
            else:
                print(f"\n{dop_type}: No valid data")
        
        # Satellite visibility
        num_sats = [d['num_satellites'] for d in dop_vals]
        print(f"\nSatellite Visibility:")
        print(f"  Mean:   {np.mean(num_sats):6.1f}")
        print(f"  Min:    {np.min(num_sats):6.0f}")
        print(f"  Max:    {np.max(num_sats):6.0f}")
    
    print("\n" + "="*80)


def save_all_plots(results, output_dir='plots', constellation_name='constellation'):
    """
    Save all DOP plots to files.
    
    Args:
        results (dict): Results dictionary with mask angles as keys
        output_dir (str): Directory to save plots
        constellation_name (str): Name of constellation for filenames
    """
    import os
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"\nSaving plots to {output_dir}/...")
    
    # Individual DOP maps
    for mask_angle in sorted(results.keys()):
        for dop_type in ['HDOP', 'VDOP', 'PDOP', 'GDOP']:
            try:
                fig = create_dop_contour_map(
                    results[mask_angle]['lats'],
                    results[mask_angle]['lons'],
                    results[mask_angle]['dop_values'],
                    dop_type=dop_type,
                    mask_angle=mask_angle
                )
                
                filename = f'{constellation_name}_{dop_type.lower()}_mask_{mask_angle}deg.png'
                filepath = os.path.join(output_dir, filename)
                fig.savefig(filepath, dpi=300, bbox_inches='tight')
                print(f"  ✓ Saved: {filename}")
                plt.close(fig)
                
            except Exception as e:
                print(f"  ✗ Error saving {dop_type} for mask {mask_angle}°: {e}")
    
    # Comparison plot
    try:
        fig = create_comparison_plot(results)
        filepath = os.path.join(output_dir, f'{constellation_name}_comparison.png')
        fig.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"  ✓ Saved: {constellation_name}_comparison.png")
        plt.close(fig)
        
        # Availability map from the highest mask angle
        try:
            max_mask = sorted(results.keys())[-1]
            fig = create_outage_map(results[max_mask]['lats'], results[max_mask]['lons'], results[max_mask]['dop_values'], title=f"Availability (valid 3D fix) — {max_mask}° mask")
            filepath = os.path.join(output_dir, f'{constellation_name}_availability_mask_{max_mask}deg.png')
            fig.savefig(filepath, dpi=300, bbox_inches='tight')
            print(f"  ✓ Saved: {constellation_name}_availability_mask_{max_mask}deg.png")
            plt.close(fig)
        except Exception as e:
            print(f"  ✗ Error saving availability map: {e}")
    except Exception as e:
        print(f"  ✗ Error saving comparison plot: {e}")
    
    print(f"\n✓ All plots saved to {output_dir}/")