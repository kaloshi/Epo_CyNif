"""
V18 Statistics SIMPLIFIED - Healthy vs MLH1 vs MSH2 (3-group)
Nur die wichtigsten Metriken, schnell und sicher
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import kruskal
from pathlib import Path

BASE_DIR = Path('python/analysis-V18')
OUTPUT_DIR = BASE_DIR / 'statistics_mutations_simple'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("="*80)
print("V18 SIMPLIFIED STATISTICS - MUTATIONS (3-GROUP)")
print("="*80)

# Load data
df_lp = pd.read_csv(BASE_DIR / 'LP' / 'lamina_propria_metrics.csv')
df_crypt = pd.read_csv(BASE_DIR / 'Crypt_IEL' / 'crypt_level_metrics_raw.csv')

# Aggregate crypts by sample
crypt_sample = df_crypt.groupby(['sample', 'mutation']).mean(numeric_only=True).reset_index()

print(f"\nâœ… Data loaded:")
print(f"   LP: {len(df_lp)} samples")
print(f"   Crypt: {len(crypt_sample)} samples (aggregated from {len(df_crypt)} crypts)")

# Key metrics to analyze
KEY_METRICS = [
    'Treg_per_CD45',
    'CD4_T_per_CD45',
    'CD8_T_per_CD45',
    'T_Cells_per_CD45',
    'CD45_per_mm2',
    'Treg_per_CD4'
]

GROUPS = ['HEALTHY', 'MLH1', 'MSH2']
GROUP_COLORS = {'HEALTHY': '#2ecc71', 'MLH1': '#e74c3c', 'MSH2': '#3498db'}

results = []

for tissue, df in [('Crypt_IEL', crypt_sample), ('LP', df_lp)]:
    print(f"\n{'='*80}")
    print(f"{tissue}")
    print("="*80)
    
    for metric in KEY_METRICS:
        if metric not in df.columns:
            continue
        
        # Get data for each group
        group_data = {}
        for group in GROUPS:
            data = df[df['mutation'] == group][metric].dropna()
            if len(data) >= 2:
                group_data[group] = data
        
        if len(group_data) < 2:
            continue
        
        # Kruskal-Wallis test (non-parametric ANOVA)
        H, p = kruskal(*group_data.values())
        
        # Simple plot
        fig, ax = plt.subplots(figsize=(8, 5))
        
        positions = list(range(1, len(group_data) + 1))
        data_list = [group_data[g] for g in GROUPS if g in group_data]
        colors_list = [GROUP_COLORS[g] for g in GROUPS if g in group_data]
        labels_list = [g for g in GROUPS if g in group_data]
        
        # Box plot
        bp = ax.boxplot(data_list, positions=positions, widths=0.5, patch_artist=True,
                        showfliers=False)
        for patch, color in zip(bp['boxes'], colors_list):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        # Points
        for i, (data, pos, color) in enumerate(zip(data_list, positions, colors_list)):
            ax.scatter([pos]*len(data), data, color=color, alpha=0.6, s=50, zorder=3)
        
        # P-value
        sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
        ax.text(np.mean(positions), ax.get_ylim()[1]*0.95, f'Kruskal-Wallis p={p:.4f} {sig}', 
                ha='center', fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        ax.set_ylabel(metric.replace('_', ' '))
        ax.set_xticks(positions)
        ax.set_xticklabels(labels_list)
        ax.set_title(f'{tissue}: {metric}')
        ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        
        # Save
        output_path = OUTPUT_DIR / f'{tissue}_{metric}.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"  âœ… {metric}: p={p:.4f} {sig}")
        
        # Save results
        result = {
            'tissue': tissue,
            'metric': metric,
            'p_value': p,
            'significant': p < 0.05
        }
        
        for group in GROUPS:
            if group in group_data:
                result[f'{group}_n'] = len(group_data[group])
                result[f'{group}_mean'] = group_data[group].mean()
                result[f'{group}_std'] = group_data[group].std()
        
        results.append(result)

# Save results
df_results = pd.DataFrame(results)
df_results.to_csv(OUTPUT_DIR / 'statistics_summary.csv', index=False)

print(f"\n{'='*80}")
print("âœ… COMPLETE!")
print(f"   Results: {OUTPUT_DIR}")
print(f"   Plots: {len(results)} created")
print(f"   Significant: {df_results['significant'].sum()}")
print("="*80)
