"""
CYCIF PIPELINE PART 7: FREQUENCY ANALYSIS V18
Version: 18.0 (V13 metrics + CSV subsets + 31px buffer + 2 compartments)
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
from shapely.geometry import shape, Point
from shapely.vectorized import contains
import warnings
warnings.filterwarnings('ignore')
from datetime import datetime

print("=" * 80)
print("CYCIF PIPELINE PART 7: FREQUENCY ANALYSIS V18")
print("=" * 80)
print(f"Version: 18.0 (V13 metrics + CSV subsets + 31px buffer + 2 compartments)")
print(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 80)

# ============================================================================
# CONFIGURATION
# ============================================================================

BASE_DIR = Path(r'C:\\Users\\researcher\data\scimap-master\data')
OUTPUT_BASE = Path('python/analysis-V18')
OUTPUT_CRYPT_IEL = OUTPUT_BASE / 'Crypt_IEL'
OUTPUT_LP = OUTPUT_BASE / 'LP'

EXCLUDED_SAMPLES = []  # V18: keine excluded samples
GROUPS_FILE = BASE_DIR / 'patient_groups.csv'
PIXEL_SIZE_UM = 0.325
BUFFER_DISTANCE_PX = 31  # V18: 31px buffer for Crypt_IEL

CD45_POPULATION = 'All_Leukocytes'
EPITHELIAL_POPULATION = 'Epithelmarker_Positive'
CD4_SUBSET = 'CD4_T'
CD8_SUBSET = 'CD8_T'
TREG_SUBSET = 'Treg'
TISSUE_RESIDENT_SUBSET = 'Tissue_Resident_T'

# V13: ALL available subsets - classified for proper metric differentiation
LEUKOCYTE_SUBSETS = [
    'All_Leukocytes', 'T_Cells', 'CD4_T', 'CD8_T', 'DN_T', 'DP_T',
    'Tissue_Resident_T', 'CD8_CD103_CD69_T', 'Treg',
    'CD4_T_naiv', 'CD8_T_naiv', 'TCRgd',
    'B_Cells', 'NK_Cells',
    'Macrophages', 'monocytes', 'Neutrophils', 'mDC',
    'prolif_Leuko'
]

NON_LEUKOCYTE_SUBSETS = [
    'Epithelmarker_Positive', 'prolif_Epi'
]

ALL_SUBSETS = LEUKOCYTE_SUBSETS + NON_LEUKOCYTE_SUBSETS

for dir_path in [OUTPUT_CRYPT_IEL, OUTPUT_LP]:
    dir_path.mkdir(parents=True, exist_ok=True)

print(f"\nâœ… Configuration loaded")
print(f"   Total subsets: {len(ALL_SUBSETS)} ({len(LEUKOCYTE_SUBSETS)} leukocyte + {len(NON_LEUKOCYTE_SUBSETS)} non-leukocyte)")

# ============================================================================
# LOAD PATIENT GROUPS
# ============================================================================

groups_df = pd.read_csv(GROUPS_FILE)
sample_to_group = dict(zip(groups_df['sample'].astype(str), groups_df['group']))
sample_to_mutation = dict(zip(groups_df['sample'].astype(str), groups_df['sub-group']))
SAMPLES = ['193', '197', '199', '201', '208', '209', '210', '212', '215', '219', '220', '221', '222']

print(f"âœ… Patient groups loaded: {len(SAMPLES)} samples")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def load_geojson(geojson_path):
    with open(geojson_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    crypt_features = []
    for feat in data.get('features', []):
        props = feat.get('properties', {})
        classification = props.get('classification', {}).get('name', '').lower()
        name = props.get('name', '').lower()
        if classification == 'crypt' or ('crypt' in name and 'non' not in name):
            crypt_features.append(feat)
    return crypt_features

def calculate_polygon_area(feature):
    poly = shape(feature['geometry'])
    area_pixels = poly.area
    area_mm2 = area_pixels * (PIXEL_SIZE_UM ** 2) / 1e6
    return area_pixels, area_mm2

def assign_cells_to_crypts_no_buffer(cells, crypt_features, x_col='X_centroid', y_col='Y_centroid'):
    cells = cells.copy()
    cells['crypt_id'] = ''
    cells['crypt_name'] = ''
    cells['crypt_index'] = -1
    
    for crypt_idx, feat in enumerate(crypt_features):
        poly = shape(feat['geometry'])
        mask = contains(poly, cells[x_col].values, cells[y_col].values)
        uuid = feat.get('id', f'unknown_{crypt_idx}')
        name = feat.get('properties', {}).get('name', f'Crypt_{crypt_idx}')
        unassigned_mask = (cells['crypt_id'] == '') & mask
        cells.loc[unassigned_mask, 'crypt_id'] = uuid
        cells.loc[unassigned_mask, 'crypt_name'] = name
        cells.loc[unassigned_mask, 'crypt_index'] = crypt_idx
    
    n_assigned = (cells['crypt_id'] != '').sum()
    print(f"   Option A: {n_assigned:,}/{len(cells):,} cells ({n_assigned/len(cells)*100:.1f}%)")
    return cells

def assign_cells_to_crypts_with_buffer(cells, crypt_features, buffer_px, x_col='X_centroid', y_col='Y_centroid'):
    cells = cells.copy()
    cells['crypt_id'] = ''
    cells['crypt_name'] = ''
    cells['crypt_index'] = -1
    cells['distance_to_crypt'] = np.inf
    
    for crypt_idx, feat in enumerate(crypt_features):
        poly = shape(feat['geometry'])
        poly_buffered = poly.buffer(buffer_px)
        mask = contains(poly_buffered, cells[x_col].values, cells[y_col].values)
        uuid = feat.get('id', f'unknown_{crypt_idx}')
        name = feat.get('properties', {}).get('name', f'Crypt_{crypt_idx}')
        
        candidate_cells = cells[mask]
        if len(candidate_cells) > 0:
            points = [Point(x, y) for x, y in zip(candidate_cells[x_col], candidate_cells[y_col])]
            distances = np.array([poly.distance(pt) for pt in points])
            for idx, dist in zip(candidate_cells.index, distances):
                if cells.loc[idx, 'crypt_id'] == '' or dist < cells.loc[idx, 'distance_to_crypt']:
                    cells.loc[idx, 'crypt_id'] = uuid
                    cells.loc[idx, 'crypt_name'] = name
                    cells.loc[idx, 'crypt_index'] = crypt_idx
                    cells.loc[idx, 'distance_to_crypt'] = dist
    
    cells = cells.drop(columns=['distance_to_crypt'])
    n_assigned = (cells['crypt_id'] != '').sum()
    print(f"   Option B: {n_assigned:,}/{len(cells):,} cells ({n_assigned/len(cells)*100:.1f}%)")
    return cells

print("âœ… Helper functions defined")

# ============================================================================
# METRIC CALCULATION (V10 - DIFFERENTIATED)
# ============================================================================

def calculate_crypt_metrics_v10(cells_in_crypt, crypt_feature):
    metrics = {}
    n_total = len(cells_in_crypt)
    metrics['n_total_cells'] = n_total
    
    if n_total == 0:
        return metrics
    
    area_pixels, area_mm2 = calculate_polygon_area(crypt_feature)
    metrics['area_pixels'] = area_pixels
    metrics['area_mm2'] = area_mm2
    
    # Count ALL subsets
    for subset in ALL_SUBSETS:
        col_name = f'is_{subset}'
        if col_name in cells_in_crypt.columns:
            metrics[f'n_{subset}'] = int(cells_in_crypt[col_name].sum())
        else:
            metrics[f'n_{subset}'] = 0
    
    n_cd45 = metrics.get(f'n_{CD45_POPULATION}', 0)
    n_epi = metrics.get(f'n_{EPITHELIAL_POPULATION}', 0)
    n_cd4 = metrics.get(f'n_{CD4_SUBSET}', 0)
    n_cd8 = metrics.get(f'n_{CD8_SUBSET}', 0)
    n_treg = metrics.get(f'n_{TREG_SUBSET}', 0)
    n_tissue_res = metrics.get(f'n_{TISSUE_RESIDENT_SUBSET}', 0)
    
    # Primary metrics
    if n_epi > 0:
        metrics['CD45_per_100_Epi'] = (n_cd45 / n_epi) * 100
        metrics['CD8_per_100_Epi'] = (n_cd8 / n_epi) * 100
        metrics['CD4_per_100_Epi'] = (n_cd4 / n_epi) * 100
        metrics['Treg_per_100_Epi'] = (n_treg / n_epi) * 100
        metrics['TissueResident_per_100_Epi'] = (n_tissue_res / n_epi) * 100
    else:
        metrics['CD45_per_100_Epi'] = 0.0
        metrics['CD8_per_100_Epi'] = 0.0
        metrics['CD4_per_100_Epi'] = 0.0
        metrics['Treg_per_100_Epi'] = 0.0
        metrics['TissueResident_per_100_Epi'] = 0.0
    
    if n_cd8 > 0:
        metrics['CD4_CD8_ratio'] = n_cd4 / n_cd8
        metrics['TissueResident_per_CD8'] = (n_tissue_res / n_cd8) * 100
    else:
        metrics['CD4_CD8_ratio'] = np.nan
        metrics['TissueResident_per_CD8'] = 0.0
    
    if n_cd4 > 0:
        metrics['Treg_per_CD4'] = (n_treg / n_cd4) * 100
    else:
        metrics['Treg_per_CD4'] = 0.0
    
    # LEUKOCYTE metrics: per_CD45, per_100_Epi, per_mm2
    for subset in LEUKOCYTE_SUBSETS:
        n_subset = metrics.get(f'n_{subset}', 0)
        
        if n_cd45 > 0:
            metrics[f'{subset}_per_CD45'] = (n_subset / n_cd45) * 100
        else:
            metrics[f'{subset}_per_CD45'] = 0.0
        
        if n_epi > 0:
            metrics[f'{subset}_per_100_Epi'] = (n_subset / n_epi) * 100
        else:
            metrics[f'{subset}_per_100_Epi'] = 0.0
        
        if area_mm2 > 0:
            metrics[f'{subset}_per_mm2'] = n_subset / area_mm2
        else:
            metrics[f'{subset}_per_mm2'] = 0.0
    
    # NON-LEUKOCYTE metrics: per_Total, per_mm2 ONLY
    for subset in NON_LEUKOCYTE_SUBSETS:
        n_subset = metrics.get(f'n_{subset}', 0)
        
        if n_total > 0:
            metrics[f'{subset}_per_Total'] = (n_subset / n_total) * 100
        else:
            metrics[f'{subset}_per_Total'] = 0.0
        
        if area_mm2 > 0:
            metrics[f'{subset}_per_mm2'] = n_subset / area_mm2
        else:
            metrics[f'{subset}_per_mm2'] = 0.0
    
    # Legacy compatibility
    metrics['Epi_per_Total'] = metrics.get('Epithelmarker_Positive_per_Total', 0.0)
    if area_mm2 > 0:
        metrics['CD45_per_mm2'] = n_cd45 / area_mm2
        metrics['CD4_per_mm2'] = n_cd4 / area_mm2
        metrics['CD8_per_mm2'] = n_cd8 / area_mm2
        metrics['Epi_per_mm2'] = n_epi / area_mm2
        metrics['Total_per_mm2'] = n_total / area_mm2
    else:
        metrics['CD45_per_mm2'] = 0.0
        metrics['CD4_per_mm2'] = 0.0
        metrics['CD8_per_mm2'] = 0.0
        metrics['Epi_per_mm2'] = 0.0
        metrics['Total_per_mm2'] = 0.0
    
    return metrics


def calculate_LP_metrics_v10(cells_lp, lp_area_mm2):
    metrics = {}
    n_total = len(cells_lp)
    metrics['n_total_cells'] = n_total
    metrics['area_mm2'] = lp_area_mm2
    
    if n_total == 0:
        return metrics
    
    # Count ALL subsets
    for subset in ALL_SUBSETS:
        col_name = f'is_{subset}'
        if col_name in cells_lp.columns:
            metrics[f'n_{subset}'] = int(cells_lp[col_name].sum())
        else:
            metrics[f'n_{subset}'] = 0
    
    n_cd45 = metrics.get(f'n_{CD45_POPULATION}', 0)
    n_epi = metrics.get(f'n_{EPITHELIAL_POPULATION}', 0)
    n_stromal = n_total - n_epi
    n_cd4 = metrics.get(f'n_{CD4_SUBSET}', 0)
    n_cd8 = metrics.get(f'n_{CD8_SUBSET}', 0)
    n_treg = metrics.get(f'n_{TREG_SUBSET}', 0)
    
    metrics['n_Stromal'] = n_stromal
    
    # Primary metrics
    if lp_area_mm2 > 0:
        metrics['CD45_per_mm2'] = n_cd45 / lp_area_mm2
        metrics['CD4_per_mm2'] = n_cd4 / lp_area_mm2
        metrics['CD8_per_mm2'] = n_cd8 / lp_area_mm2
    else:
        metrics['CD45_per_mm2'] = 0.0
        metrics['CD4_per_mm2'] = 0.0
        metrics['CD8_per_mm2'] = 0.0
    
    if n_cd8 > 0:
        metrics['CD4_CD8_ratio'] = n_cd4 / n_cd8
    else:
        metrics['CD4_CD8_ratio'] = np.nan
    
    if n_cd4 > 0:
        metrics['Treg_per_CD4'] = (n_treg / n_cd4) * 100
    else:
        metrics['Treg_per_CD4'] = 0.0
    
    # LEUKOCYTE metrics: per_CD45, per_mm2, per_Stromal
    for subset in LEUKOCYTE_SUBSETS:
        n_subset = metrics.get(f'n_{subset}', 0)
        
        if n_cd45 > 0:
            metrics[f'{subset}_per_CD45'] = (n_subset / n_cd45) * 100
        else:
            metrics[f'{subset}_per_CD45'] = 0.0
        
        if lp_area_mm2 > 0:
            metrics[f'{subset}_per_mm2'] = n_subset / lp_area_mm2
        else:
            metrics[f'{subset}_per_mm2'] = 0.0
        
        if n_stromal > 0:
            metrics[f'{subset}_per_Stromal'] = (n_subset / n_stromal) * 100
        else:
            metrics[f'{subset}_per_Stromal'] = 0.0
    
    # NON-LEUKOCYTE metrics: per_Total, per_mm2 ONLY
    for subset in NON_LEUKOCYTE_SUBSETS:
        n_subset = metrics.get(f'n_{subset}', 0)
        
        if n_total > 0:
            metrics[f'{subset}_per_Total'] = (n_subset / n_total) * 100
        else:
            metrics[f'{subset}_per_Total'] = 0.0
        
        if lp_area_mm2 > 0:
            metrics[f'{subset}_per_mm2'] = n_subset / lp_area_mm2
        else:
            metrics[f'{subset}_per_mm2'] = 0.0
    
    # Legacy compatibility
    if n_stromal > 0:
        metrics['CD45_per_Stromal'] = (n_cd45 / n_stromal) * 100
    else:
        metrics['CD45_per_Stromal'] = 0.0
    
    return metrics

print("âœ… V10 metric functions defined (differentiated for leukocytes vs non-leukocytes)")

# ============================================================================
# MAIN ANALYSIS
# ============================================================================

all_crypt_iel_metrics = []
all_lp_metrics = []

print("\n" + "=" * 80)
print("PROCESSING SAMPLES")
print("=" * 80)

for sample_id in SAMPLES:
    print(f"\n{'='*80}")
    print(f"ðŸ“Š SAMPLE {sample_id}")
    print("=" * 80)
    
    group = sample_to_group.get(sample_id, 'UNKNOWN')
    mutation = sample_to_mutation.get(sample_id, 'UNKNOWN')
    print(f"   Group: {group}, Mutation: {mutation}")
    
    sample_dir = BASE_DIR / f'sample_{sample_id}'
    parquet_path = sample_dir / 'cylinter' / 'customized_cutoff' / 'gating.parquet'
    geojson_path = sample_dir / 'crypt_segmentation' / 'simplified' / 'correct_format' / f'{sample_id} simplified.geojson'
    
    if not parquet_path.exists() or not geojson_path.exists():
        print(f"   âŒ Files not found")
        continue
    
    print(f"   ðŸ“‚ Loading data...")
    cells = pd.read_parquet(parquet_path)
    print(f"      {len(cells):,} cells")
    
    crypt_features = load_geojson(geojson_path)
    print(f"      {len(crypt_features)} crypts")
    
    if len(crypt_features) == 0:
        print(f"   âš ï¸  No crypts!")
        continue
    
    total_crypt_area_mm2 = sum(calculate_polygon_area(f)[1] for f in crypt_features)
    x_min, x_max = cells['X_centroid'].min(), cells['X_centroid'].max()
    y_min, y_max = cells['Y_centroid'].min(), cells['Y_centroid'].max()
    total_tissue_area_mm2 = ((x_max - x_min) * (y_max - y_min)) * (PIXEL_SIZE_UM ** 2) / 1e6
    lp_area_mm2 = total_tissue_area_mm2 - total_crypt_area_mm2
    
    # V18: Crypt_IEL (31px buffer)
    print(f"\n   ðŸ”¹ CRYPT_IEL (31px buffer):")
    cells_assigned = assign_cells_to_crypts_with_buffer(cells, crypt_features, BUFFER_DISTANCE_PX)
    crypts_in_sample = cells_assigned[cells_assigned['crypt_id'] != '']['crypt_id'].unique()
    
    for crypt_id in crypts_in_sample:
        cells_in_crypt = cells_assigned[cells_assigned['crypt_id'] == crypt_id]
        crypt_idx = cells_in_crypt['crypt_index'].iloc[0]
        crypt_feature = crypt_features[crypt_idx]
        crypt_name = cells_in_crypt['crypt_name'].iloc[0]
        
        metrics = calculate_crypt_metrics_v10(cells_in_crypt, crypt_feature)
        metrics['crypt_id'] = crypt_id
        metrics['crypt_name'] = crypt_name
        metrics['sample'] = sample_id
        metrics['group'] = group
        metrics['mutation'] = mutation
        metrics['compartment'] = 'Crypt_IEL'
        all_crypt_iel_metrics.append(metrics)
    
    # V18: Lamina Propria (exclude Crypt_IEL)
    print(f"   ðŸ”¹ LP:")
    cells_lp = cells_assigned[cells_assigned['crypt_id'] == '']
    print(f"      {len(cells_lp):,} cells, {lp_area_mm2:.2f} mmÂ²")
    
    lp_metrics = calculate_LP_metrics_v10(cells_lp, lp_area_mm2)
    lp_metrics['sample'] = sample_id
    lp_metrics['group'] = group
    lp_metrics['mutation'] = mutation
    all_lp_metrics.append(lp_metrics)
    
    print(f"   âœ… Complete")

print(f"\n{'='*80}")
print("âœ… ALL SAMPLES PROCESSED")
print("=" * 80)

# ============================================================================
# SAVE RESULTS
# ============================================================================

print("\n" + "=" * 80)
print("SAVING RESULTS")
print("=" * 80)

df_crypt_iel = pd.DataFrame(all_crypt_iel_metrics)
output_file = OUTPUT_CRYPT_IEL / 'crypt_level_metrics_raw.csv'
df_crypt_iel.to_csv(output_file, index=False)
print(f"\nâœ… Crypt_IEL: {len(df_crypt_iel)} rows, {len(df_crypt_iel.columns)} columns")
print(f"   {output_file}")

df_lp = pd.DataFrame(all_lp_metrics)
output_file = OUTPUT_LP / 'lamina_propria_metrics.csv'
df_lp.to_csv(output_file, index=False)
print(f"\nâœ… LP: {len(df_lp)} rows, {len(df_lp.columns)} columns")
print(f"   {output_file}")

# ============================================================================
# SUMMARY
# ============================================================================

print("\n" + "=" * 80)
print("SUMMARY STATISTICS (V13)")
print("=" * 80)

print("\nðŸ“Š PRIMARY METRICS - Crypts:")
print(f"   CD45/100 Epi: {df_crypts_A['CD45_per_100_Epi'].mean():.1f} Â± {df_crypts_A['CD45_per_100_Epi'].std():.1f}")
print(f"   CD4/CD8 ratio: {df_crypts_A['CD4_CD8_ratio'].mean():.2f} Â± {df_crypts_A['CD4_CD8_ratio'].std():.2f}")

print("\nðŸ“Š PRIMARY METRICS - LP:")
print(f"   CD45/mmÂ²: {df_lp['CD45_per_mm2'].mean():.0f} Â± {df_lp['CD45_per_mm2'].std():.0f}")
print(f"   CD4/CD8 ratio: {df_lp['CD4_CD8_ratio'].mean():.2f} Â± {df_lp['CD4_CD8_ratio'].std():.2f}")

print("\nðŸ“Š V13 IMPROVEMENTS:")
print(f"   prolif_Leuko/CD45: {df_lp['prolif_Leuko_per_CD45'].mean():.1f}% (was 100% in V12!) âœ…")
print(f"   Total subsets: {len(ALL_SUBSETS)} (was 17 in V12)")
leuko_per_cd45 = [c for c in df_lp.columns if '_per_CD45' in c]
non_leuko_per_total = [c for c in df_lp.columns if '_per_Total' in c]
print(f"   Leukocyte metrics (_per_CD45): {len(leuko_per_cd45)}")
print(f"   Non-leukocyte metrics (_per_Total): {len(non_leuko_per_total)}")

print("\n" + "=" * 80)
print("âœ… V13 ANALYSIS COMPLETE")
print("=" * 80)
print(f"\nTimestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
