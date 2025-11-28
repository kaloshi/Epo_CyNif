# Epoxy_CyNif Pipeline v1.0

🔬 **Epoxy-based Cyclic Nanobody Immunofluorescence (Epoxy_CyNif)** – End-to-end workflow for multicycle staining/imaging: illumination correction, stitching/registration/deconvolution+EDF, spillover/autofluorescence removal, segmentation, GUI-based QC (CyLinter), and spatial statistics.

📖 **[View Documentation](https://kaloshi.github.io/Epo_CyNif/)**

## TL;DR
- **Inputs:** Multicycle TIFF/CZI stacks plus marker CSV per sample.
- **Paths:** `BASE_EXPORT = C:\Users\researcher\data\Epoxy_CyNif\data\export\<sample>`.
- **Workflow:** Part 1 → 2 → 3 → 4a/4b → 5 (CyLinter QC) → 6/7 (Statistics).
- **Day 2+:** In Part 5, run cells 2–3 + extended GUI (Cell 12), delete checkpoints selectively, restart from desired module.
- **Outputs:** Segmentation/Processed/Clustering in `BASE_EXPORT`, CyLinter checkpoints/reports in `cylinter_output_prune_test/`, statistics in `python/analysis-V18/`.

## Pipeline Modules

### 1) Part 1 – Illumination/Stitching/Registration/Decon&EDF
`Illumination_correction,Stiching,Registration,Decon&EDF/Cycif_pipeline_part_1_stiching_REG_decon_EDF.ipynb`  

Illumination correction using **BaSiCpy** to remove shading artifacts, followed by stitching and registration using **Ashlar**, deconvolution (Richardson-Lucy), and Extended Depth of Field (EDF) processing. Generates fused multichannel stacks.

*References:*
- BaSiCpy: [DOI: 10.1038/ncomms14836](https://doi.org/10.1038/ncomms14836)
- Ashlar: [DOI: 10.1093/bioinformatics/btac544](https://doi.org/10.1093/bioinformatics/btac544)
- Deconvolution: [DOI: 10.1364/JOSA.62.000055](https://doi.org/10.1364/JOSA.62.000055), [DOI: 10.1086/111605](https://doi.org/10.1086/111605)

### 2) Part 2 – Spillover Removal
`Spillover&AF_removal/Cycif_pipeline_part_2_Spillover.ipynb`  

Supervised linear compensation workflow to correct for signal spillover (channel crosstalk). Compensation coefficients are estimated from single-positive ROIs using non-negative least squares and applied linearly to subtract donor signals from target channels.

*Reference:*
- Spectral Compensation: [DOI: 10.1111/j.1749-6632.1993.tb38775.x](https://doi.org/10.1111/j.1749-6632.1993.tb38775.x)

### 3) Part 3 – Autofluorescence Removal (ACE)
`Spillover&AF_removal/Cycif_pipeline_part_3_AF_REV_ACE.ipynb`  

Autofluorescence (AF) removal by subtracting non-specific signals from dedicated background channels using robust linear regression (Huber regressor), adapted from the **AutoSpill** algorithm. Includes multi-scale adaptive unsharp masking (ACE) for local contrast enhancement.

*Reference:*
- AutoSpill: [DOI: 10.1038/s41467-021-23126-8](https://doi.org/10.1038/s41467-021-23126-8)

### 4a) Segmentation
`Segmentation/Cycif_pipeline_part_4a_segmentation.ipynb`  

Nuclei and cellular boundaries segmented using **InstanSeg** with the `fluorescence_nuclei_and_cells` model. Employs a channel-invariant architecture (ChannelNet) integrating nuclear (DAPI) and membrane markers. Feature extraction uses marker-specific quantification: 95th percentile for membrane markers, mean intensity for cytoplasmic/nuclear markers.

*References:*
- InstanSeg: [DOI: 10.1101/2024.09.04.611150](https://doi.org/10.1101/2024.09.04.611150)
- Marker quantification strategy: [DOI: 10.1101/2022.06.08.495346](https://doi.org/10.1101/2022.06.08.495346)

### 4b) Batch Segmentation – Micro-SAM
`Segmentation/Cycif_pipeline_part_4b_batch_MICROSAM_ROBUST_V8.ipynb`  

Histological compartment segmentation (intestinal crypts) using **Micro-SAM** (Segment Anything Model optimized for microscopy). Applied to Laminin channel with adaptive tiling strategy. Morphological filtering selects valid crypts based on geometric criteria (area >8,000 px; eccentricity <0.9).

*Reference:*
- Micro-SAM: [DOI: 10.1038/s41592-024-02580-4](https://doi.org/10.1038/s41592-024-02580-4)

### 5) CyLinter QC (GUI)
`Cylinter/Cycif_pipeline_part_5_cylinter.ipynb`  

Interactive quality control using **CyLinter** pipeline. Validates feature data, establishes gating cutoffs, and applies combinatorial subset definitions for cell phenotyping.

- **Day 1:** Complete run with marker selection + pipeline start.  
- **Day 2+:** Extended GUI, delete checkpoints, restart from module X.

*Reference:*
- CyLinter: [DOI: 10.1038/s41592-024-02328-0](https://doi.org/10.1038/s41592-024-02328-0)

### 6/7) Spatial Analysis & Statistics (v18)
- `Spatial_analysis/Cycif_pipeline_part_6_run_data_frames.py`  
- `Spatial_analysis/Cycif_pipeline_part_7_statistics_v18_2groups_CORRECT.py`  
- `Spatial_analysis/Cycif_pipeline_part_7_statistics_healthy_vs_mutations_v18_SIMPLE.py`  

Spatial analysis contextualizes cellular phenotypes within tissue architecture. Analysis restricted to orthogonal crypt cross-sections to avoid stereological bias. Cells assigned to compartments (Crypt-associated with 31px buffer ~10µm, or Lamina Propria). Statistics use non-parametric tests (Mann-Whitney U, Kruskal-Wallis with Dunn's post-hoc, Benjamini-Hochberg correction). Patient-level aggregation prevents pseudoreplication.

*References:*
- Stereological considerations: [DOI: 10.1590/S0100-69912012000200010](https://doi.org/10.1590/S0100-69912012000200010)
- Pseudoreplication correction: [DOI: 10.1186/1471-2202-11-5](https://doi.org/10.1186/1471-2202-11-5)

## Data Layout
- **Sample root:** `BASE_EXPORT = C:\Users\researcher\data\Epoxy_CyNif\data\export\<sample>`
- **Markers:** `markers_<SAMPLE_ID>.csv` in sample folder.
- **Outputs:**  
  - Segmentation/Processed/Clustering → `BASE_EXPORT`  
  - CyLinter Checkpoints/Reports → `cylinter_output_prune_test/`  
  - Statistics → `python/analysis-V18/`

## Workflow
- **Day 1 (Initial run):** Part 1 → Part 2 → Part 3 → Part 4a/4b → Part 5 (CyLinter) → Part 6/7.
- **Day 2+ / Error correction:**  
  - Restart Part 5 kernel, run cells 2–3.  
  - Open extended GUI (Cell 12), delete checkpoint(s), select start module, run.  
  - Interactive GUIs: `selectROIs`, `setContrast`, `gating`. Markers without threshold will be re-gated.

## Requirements
```
python 3.9.x–3.11.x
pandas>=2.2.0
numpy>=1.26.0
scipy>=1.12.0
scikit-image>=0.22.0
tifffile>=2024.1.0
shapely>=2.0.0
matplotlib>=3.8.0
seaborn>=0.13.0
ipywidgets>=8.1.0
pyyaml>=6.0.0
instanseg  # Model: fluorescence_nuclei_and_cells
cylinter   # CLI/GUI
napari     # Optional: visual inspection
```

See [requirements.txt](requirements.txt) for full dependency list.

## Installation
```bash
pip install -r requirements.txt
```

## QC / Checkpoints (CyLinter)
- **Checkpoints:** `cylinter_output_prune_test/checkpoints/*.parquet`
- **Gating thresholds:** `cylinter_report.yml` (persists across runs)
- `--module X` sets start point; end point is always pipeline end; checkpoints skip completed modules.

## Key Features
- Reconstructs high-dimensional cyclic IF nanobody stacks (Epoxy_CyNif) across all cycles.
- Removes technical artifacts (illumination, spillover, autofluorescence) for reproducible quantification.
- Performs robust cell/nucleus segmentation with DAPI optimization and exports feature tables.
- Provides interactive QC (CyLinter) with checkpointing for fast day-2+ iterations.
- Computes spatial frequencies and group statistics for downstream analyses.

## Citation
If you use this pipeline, please cite the relevant tools:

### Image Processing (Part 1)
- **BaSiCpy** (Illumination correction): Peng et al. (2017) *A BaSiC tool for background and shading correction of optical microscopy images.* [DOI: 10.1038/ncomms14836](https://doi.org/10.1038/ncomms14836)
- **Ashlar** (Stitching/Registration): Muhlich et al. (2022) *Stitching and registering highly multiplexed whole-slide images of tissues and tumors using ASHLAR.* [DOI: 10.1093/bioinformatics/btac544](https://doi.org/10.1093/bioinformatics/btac544)

### Spectral Correction (Parts 2–3)
- **Spectral Compensation**: Roederer (1993) *Spectral Compensation for Flow Cytometry.* [DOI: 10.1111/j.1749-6632.1993.tb38775.x](https://doi.org/10.1111/j.1749-6632.1993.tb38775.x)
- **AutoSpill** (AF removal): Roca et al. (2021) *AutoSpill is a principled framework for rapid, robust, and accurate compensation.* [DOI: 10.1038/s41467-021-23126-8](https://doi.org/10.1038/s41467-021-23126-8)

### Segmentation (Parts 4–5)
- **InstanSeg**: Goldsborough et al. (2024) *A novel channel invariant architecture for the segmentation of cells and nuclei in multiplexed images.* [DOI: 10.1101/2024.09.04.611150](https://doi.org/10.1101/2024.09.04.611150)
- **Micro-SAM**: Archit et al. (2024) *Segment Anything for Microscopy.* [DOI: 10.1038/s41592-024-02580-4](https://doi.org/10.1038/s41592-024-02580-4)
- **CyLinter**: Baker et al. (2024) *CyLinter: Interactive quality control for multiplexed tissue imaging.* [DOI: 10.1038/s41592-024-02328-0](https://doi.org/10.1038/s41592-024-02328-0)
- **Marker quantification**: Schapiro et al. (2022) *MCMICRO: A scalable, modular image-processing pipeline.* [DOI: 10.1101/2022.06.08.495346](https://doi.org/10.1101/2022.06.08.495346)

### Statistical Methods (Parts 6–7)
- **Stereology**: Teixeira et al. (2012) *Corpuscle problem in stereology.* [DOI: 10.1590/S0100-69912012000200010](https://doi.org/10.1590/S0100-69912012000200010)
- **Pseudoreplication**: Lazic (2010) *The problem of pseudoreplication in neuroscientific studies.* [DOI: 10.1186/1471-2202-11-5](https://doi.org/10.1186/1471-2202-11-5)

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments
Third-party packages (InstanSeg, CyLinter, scikit-image, tifffile, shapely, pandas, numpy, seaborn, ipywidgets, pyyaml, napari) are used under their respective licenses (MIT/BSD/Apache).
