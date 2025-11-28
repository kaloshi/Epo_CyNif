# Epoxy_CyNif Pipeline v1.0

Epoxy-based Cyclic Nanobody Immunofluorescence (Epoxy_CyNif) – End-to-end Workflow für mehrzyklisches Staining/Imaging: Illumination-Korrektur, Stitching/Registration/Decon+EDF, Spillover/AF-Removal, Segmentierung, GUI-basierte QC (CyLinter) und räumliche Statistik.

## TL;DR
- Eingaben: Mehrzyklische TIFF/CZI-Stapel plus Marker-CSV pro Sample.
- Pfade: `BASE_EXPORT = C:\Users\researcher\data\Epoxy_CyNif\Epoxy_CyNif\data\export\<sample>`.
- Ablauf: Part 1 → 2 → 3 → 4a/4b → 5 (CyLinter QC) → 6/7 (Statistik).
- Tag 2+: In Part 5 nur Zellen 2–3 + erweiterte GUI (Cell 12), Checkpoints gezielt löschen, ab gewünschtem Modul neu starten.
- Outputs: Segmentation/Processed/Clustering in `BASE_EXPORT`, CyLinter-Checkpoints/-Reports in `cylinter_output_prune_test/`, Statistik in `python/analysis-V18/`.

## Inhalte (Module)
1) **Part 1 – Illumination/Stitching/Registration/Decon&EDF**  
   `Illumination_correction,Stiching,Registration,Decon&EDF/Cycif_pipeline_part_1_stiching_REG_decon_EDF.ipynb`  
   Illumination-Korrektur, Stitching, Registrierung, Deconvolution/EDF; erzeugt fused multichannel Stacks.

2) **Part 2 – Spillover Removal**  
   `Spillover&AF_removal/Cycif_pipeline_part_2_Spillover.ipynb`  
   Mutual-Information Spillover-Schätzung, 0–1 Normalisierung, Rücktransformation in Original-Datentyp.

3) **Part 3 – AF Removal (ACE)**  
   `Spillover&AF_removal/Cycif_pipeline_part_3_AF_REV_ACE.ipynb`  
   Autofluoreszenz-Bereinigung auf dem Stack.

4a) **Segmentation**  
   `Segmentation/Cycif_pipeline_part_4a_segmentation.ipynb`  
   `SAMPLE_ID` setzen, Marker-CSV aus `BASE_EXPORT` laden, DAPI-Optimierung (beste Preprocessing-Variante), InstanSeg Nuclei/Cells, Feature-Exports.

4b) **Batch Segmentation Helper**  
   `Segmentation/Cycif_pipeline_part_4b_batch_MICROSAM_ROBUST_V8.ipynb`  
   Batch/robuste Segmentierung für mehrere Samples.

5) **CyLinter QC (GUI)**  
   `Cylinter/Cycif_pipeline_part_5_cylinter.ipynb`  
   Tag 1: kompletter Lauf mit Marker-Auswahl + Pipeline-Start.  
   Tag 2+: erweiterte GUI, Checkpoints löschen, ab Modul X neu starten.  
   Nutzt `cylinter_config.yml` + `markers.csv` im Notebook-Ordner; Checkpoints/Reports in `cylinter_output_prune_test/`.

6/7) **Spatial Analysis & Statistics (v18)**  
   `Spatial_analysis/Cycif_pipeline_part_6_run_data_frames.py`  
   `Spatial_analysis/Cycif_pipeline_part_7_statistics_v18_2groups_CORRECT.py`  
   `Spatial_analysis/Cycif_pipeline_part_7_statistics_healthy_vs_mutations_v18_SIMPLE.py`  
   Frequenzmetriken und Gruppen-Statistiken für LP/Crypt und Mutationsgruppen; erwartet Inputs unter `python/analysis-V18/`.

## Datenlayout
- Root pro Sample: `BASE_EXPORT = C:\Users\researcher\data\Epoxy_CyNif\Epoxy_CyNif\data\export\<sample>`
- Marker: `markers_<SAMPLE_ID>.csv` im Sample-Ordner.  
  (Legacy-Fallback in Part 4a Cell 9: ggf. Kopie mit „193“ im Dateinamen bereitstellen.)
- Outputs:  
  - Segmentation/Processed/Clustering in `BASE_EXPORT`  
  - CyLinter Checkpoints/Reports: `cylinter_output_prune_test/`  
  - Statistik: `python/analysis-V18/`

## Arbeitsablauf
- **Tag 1 (Erstlauf)**: Part 1 → Part 2 → Part 3 → Part 4a/4b → Part 5 (CyLinter Tag 1) → Part 6/7.
- **Tag 2+ / Fehlerkorrektur**:  
  - Part 5 Kernel neu starten, Zellen 2–3 ausführen.  
  - Erweiterte GUI (Cell 12) öffnen, Checkpoint(s) löschen, Startmodul wählen, laufen lassen.  
  - Interactive GUIs: `selectROIs`, `setContrast`, `gating`. Marker ohne Threshold werden erneut gegated.

## Anforderungen / Pakete
- python 3.9.x–3.11.x
- pandas 2.2.x; numpy 1.26.x; scipy 1.12.x
- scikit-image 0.22.x; tifffile 2024.x; shapely 2.0.x
- matplotlib 3.8.x; seaborn 0.13.x
- ipywidgets 8.1.x; pyyaml 6.0.x
- instanseg (Model `fluorescence_nuclei_and_cells`, aktueller Release)
- cylinter (CLI/GUI, aktueller Release)
- napari (optional für visuelle Kontrolle)

## Hinweise zur QC / Checkpoints (CyLinter)
- Checkpoints: `cylinter_output_prune_test/checkpoints/*.parquet`
- Gating-Thresholds: `cylinter_report.yml` (bleiben über Läufe erhalten)
- `--module X` setzt den Startpunkt, Endpunkt ist immer das Pipeline-Ende; Checkpoints überspringen bereits erledigte Module.

## Zweck / Nutzen
- Rekonstruiert hochdimensionale Cyclic IF Nanobody-Stapel (Epoxy_CyNif) über alle Zyklen.
- Entfernt technische Artefakte (Illumination, Spillover, AF) für reproduzierbare Quantifizierung.
- Führt robuste Zell-/Kern-Segmentierung mit DAPI-Optimierung durch und exportiert Feature-Tabellen.
- Bietet interaktive QC (CyLinter) mit Checkpointing für schnelle Tag-2+ Iterationen.
- Berechnet räumliche Häufigkeiten und Gruppen-Statistiken für Downstream-Analysen.

## Rechtliches / Attribution
- Drittpakete (InstanSeg, CyLinter, scikit-image, tifffile, shapely, pandas, numpy, seaborn, ipywidgets, pyyaml, napari) entsprechend ihrer Lizenzen nennen (typisch MIT/BSD/Apache; GPL/AGPL beachten).
- Keine fremden Daten ohne Freigabe veröffentlichen; Modell-/Paper-Zitationen gemäß Tool-Doku ergänzen.
- Eigenes Lizenzfile (z.B. MIT) und Requirements/Env-Datei für GitHub beilegen.
