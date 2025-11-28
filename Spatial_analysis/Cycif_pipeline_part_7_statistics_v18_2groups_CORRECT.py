"""
V18 Statistics: HEALTHY vs LYNCH (2 groups)
- Bar plots with individual points
- p-value top right
- All normalizations (per_CD45, per_Total, per_100_Epi, per_mm2)
- All subsets
- Output: CSV, PNG, SVG folders
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pathlib import Path
import numpy as np

# Configuration
OUTPUT_BASE = Path('python/analysis-V18/statistics_2groups')
OUTPUT_BASE.mkdir(parents=True, exist_ok=True)

# Create subfolders
(OUTPUT_BASE / 'CSV').mkdir(exist_ok=True)
(OUTPUT_BASE / 'PNG').mkdir(exist_ok=True)
(OUTPUT_BASE / 'SVG').mkdir(exist_ok=True)

print("="*80)
print("V18 STATISTICS: HEALTHY vs LYNCH (2 GROUPS)")
print("="*80)

# Load data
df_lp = pd.read_csv('python/analysis-V18/LP/lamina_propria_metrics.csv')
df_crypt = pd.read_csv('python/analysis-V18/Crypt_IEL/crypt_level_metrics_raw.csv')

# Aggregate crypt data by sample
crypt_agg = df_crypt.groupby('sample').mean(numeric_only=True).reset_index()

# Get metadata from crypt data
metadata = df_crypt[['sample', 'group', 'mutation']].drop_duplicates()

# Merge metadata with aggregated crypt data
crypt_agg = crypt_agg.merge(metadata, on='sample', how='left')

# Ensure group is uppercase (LP already has group/mutation columns)
df_lp['group'] = df_lp['group'].str.upper()
crypt_agg['group'] = crypt_agg['group'].str.upper()

print(f"\nâœ… Data loaded:")
print(f"   LP: {len(df_lp)} samples")
print(f"   Crypt_IEL: {len(crypt_agg)} samples")

# Define metrics to analyze
NORMALIZATIONS = ['per_CD45', 'per_Total', 'per_100_Epi', 'per_mm2']
SUBSETS = ['Treg', 'CD4_T', 'CD8_T', 'T_Cells', 'B_Cells', 'NK_Cells', 
           'Macrophages', 'DN_T', 'DP_T', 'Tissue_Resident_T']

def create_bar_plot_with_points(data, metric, tissue, output_prefix):
    """Create bar plot with individual points and p-value"""
    
    # Get data for both groups
    healthy = data[data['group'] == 'HEALTHY'][metric].dropna()
    lynch = data[data['group'] == 'LYNCH'][metric].dropna()
    
    if len(healthy) == 0 or len(lynch) == 0:
        return None
    
    # Statistical test
    stat, pval = stats.mannwhitneyu(healthy, lynch, alternative='two-sided')
    
    # Create figure
    fig, ax = plt.subplots(figsize=(6, 5))
    
    # Calculate means and SEM
    groups = ['HEALTHY', 'LYNCH']
    means = [healthy.mean(), lynch.mean()]
    sems = [healthy.sem(), lynch.sem()]
    
    # Bar plot
    bars = ax.bar(groups, means, yerr=sems, capsize=5, 
                   color=['#3498db', '#e74c3c'], alpha=0.7, edgecolor='black')
    
    # Add individual points
    for i, (group_name, group_data) in enumerate([('HEALTHY', healthy), ('LYNCH', lynch)]):
        x = np.random.normal(i, 0.04, size=len(group_data))
        ax.scatter(x, group_data, alpha=0.6, s=50, color='black', zorder=3)
    
    # Add p-value top right
    sig_text = f"p={pval:.4f}"
    if pval < 0.001:
        sig_text += " ***"
    elif pval < 0.01:
        sig_text += " **"
    elif pval < 0.05:
        sig_text += " *"
    else:
        sig_text += " ns"
    
    ax.text(0.95, 0.95, sig_text, transform=ax.transAxes, 
            ha='right', va='top', fontsize=12, 
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Labels
    metric_name = metric.replace('_', ' ').title()
    ax.set_ylabel(metric_name, fontsize=12)
    ax.set_title(f'{tissue}: {metric_name}', fontsize=14, fontweight='bold')
    ax.set_xlabel('')
    
    # Grid
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    plt.tight_layout()
    
    # Save
    png_path = OUTPUT_BASE / 'PNG' / f'{output_prefix}_{metric}_{tissue}.png'
    svg_path = OUTPUT_BASE / 'SVG' / f'{output_prefix}_{metric}_{tissue}.svg'
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close()
    
    return {'metric': metric, 'tissue': tissue, 'p_value': pval, 
            'healthy_mean': healthy.mean(), 'lynch_mean': lynch.mean(),
            'healthy_n': len(healthy), 'lynch_n': len(lynch)}

# Process all metrics
results = []

print("\n" + "="*80)
print("PROCESSING METRICS")
print("="*80)

for tissue, df in [('Crypt_IEL', crypt_agg), ('LP', df_lp)]:
    print(f"\n{tissue}:")
    
    for subset in SUBSETS:
        for norm in NORMALIZATIONS:
            metric = f'{subset}_{norm}'
            
            if metric in df.columns:
                result = create_bar_plot_with_points(df, metric, tissue, 'healthy_vs_lynch')
                if result:
                    results.append(result)
                    sig = "***" if result['p_value'] < 0.001 else "**" if result['p_value'] < 0.01 else "*" if result['p_value'] < 0.05 else "ns"
                    print(f"  âœ… {metric}: p={result['p_value']:.4f} {sig}")

# Save results to CSV
df_results = pd.DataFrame(results)
df_results.to_csv(OUTPUT_BASE / 'CSV' / 'statistical_results.csv', index=False)

print("\n" + "="*80)
print("âœ… COMPLETE!")
print("="*80)
print(f"\nOutput: {OUTPUT_BASE}")
print(f"  CSV: {len(results)} results saved")
print(f"  PNG: {len(results)} plots created")
print(f"  SVG: {len(results)} plots created")
print(f"  Significant (p<0.05): {(df_results['p_value'] < 0.05).sum()}")
