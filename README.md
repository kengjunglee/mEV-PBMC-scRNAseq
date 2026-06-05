# mEV-PBMC-scRNAseq

Single-cell RNA-seq analysis of human peripheral blood mononuclear cells (PBMCs)
treated *ex vivo* with **mechanically generated extracellular vesicles (mEV)**,
**canonical secreted EVs (sEV)**, or **PBS** (vehicle control).

Four donors were each split across the three conditions, hashtag-labeled (HTO),
and profiled on the 10x Genomics platform. The pipeline covers QC, HTO
demultiplexing, SCTransform normalization, CCA integration, marker-based
annotation, pseudobulk DESeq2 differential expression (mEV vs sEV), GO
enrichment, and Azimuth reference mapping.

## Background

Circulating tumor cells (CTCs) are the seeds of metastasis, yet most are
mechanically trapped and cleared in the microvasculature before entering
peripheral circulation. CTCs arrested in the microvasculature rapidly release
abundant **mechanically generated extracellular vesicles (mEVs)** that are
molecularly distinct from canonical secreted EVs (sEVs). Despite their
predominantly cytoplasmic origin, mEVs are enriched in the nuclear
damage-associated molecular pattern **HMGB1**. HMGB1-positive mEVs potently
activate monocytes, induce neutrophil extracellular trap (NET) formation, and
promote lethal thrombosis *in vivo*; HMGB1 depletion attenuates these effects.
In pancreatic cancer patients, plasma EV–associated HMGB1 is elevated and
correlates with metastasis, thrombotic events, and D-dimer.

This experiment tests the **monocyte-activating** arm of that axis directly: by
treating human PBMCs with mEVs versus sEVs versus vehicle and profiling the
single-cell transcriptional response, with particular attention to the
monocyte and dendritic-cell compartments. The pseudobulk DE and GO-enrichment
steps below focus on the mEV-vs-sEV contrast within those populations.

## Experimental design

| Sample | Donor | Condition |
|--------|-------|-----------|
| 606M | 606 | mEV |
| 606P | 606 | PBS |
| 606S | 606 | sEV |
| 610M | 610 | mEV |
| 610P | 610 | PBS |
| 610S | 610 | sEV |
| 706M | 706 | mEV |
| 706P | 706 | PBS |
| 706S | 706 | sEV |
| 746M | 746 | mEV |
| 746P | 746 | PBS |
| 746S | 746 | sEV |

Twelve libraries total: 4 donors × 3 conditions. Each library carries a Gene
Expression and an Antibody Capture (HTO) assay.

## Pipeline

| Step | Method / tool | Key parameters |
|------|---------------|----------------|
| Load | `Read10X` | `filtered_feature_bc_matrix` per sample |
| QC metrics | `PercentageFeatureSet` | mitochondrial `^MT-` |
| HTO demux | `HTODemux` | `positive.quantile = 0.99`, CLR-normalized HTO |
| QC filter | `subset` | `200 < nFeature_RNA < 8000`, `percent.mt < 5` |
| Normalization | `SCTransform` | default |
| Dim. reduction | `RunPCA`, `RunUMAP` | `dims = 1:30` |
| Integration | `IntegrateLayers` (CCA) | `normalization.method = "SCT"` |
| Clustering | `FindNeighbors` / `FindClusters` | `dims = 1:10`, `resolution = 0.8` |
| Cluster markers | `PrepSCTFindMarkers` → `FindAllMarkers` | top 10 by `avg_log2FC > 1` |
| Pseudobulk DE | `AggregateExpression` → `FindMarkers` | DESeq2, mEV vs sEV |
| Enrichment | `enrichGO` (clusterProfiler) | `OrgDb = org.Hs.eg.db`, BP |
| Reference map | `RunAzimuth` | `reference = "pbmcref"` |

## Cell-type annotation

Clusters (resolution 0.8) were annotated from canonical lineage markers into:

Memory CD4 T, Naive CD4 T, CD8 T, MAIT T, NK, NK T, B, Monocytes, and DC.

Markers used include `CD3D/E/G`, `CD4`, `IL7R`, `CCR7`, `SELL` (CD4 T);
`CD8A/B`, `GZMB`, `NKG7` (CD8 T); `SLC4A10` (MAIT); `GNLY`, `KLRF1`, `XCL1`
(NK / NK T); `MS4A1`, `CD79A/B`, `CD19` (B); `CD14`, `LYZ`, `CD68`, `FCGR3A`
(Monocytes); and `FCER1A`, `CST3`, `IL3RA` (DC).

## Repository structure

```
.
├── README.md                  # this file
├── mEV-PBMC.Rmd               # full analysis notebook
├── .gitignore
├── data/                      # NOT tracked — see "Data" below
│   └── post_aggr/<sample>/outs/count/filtered_feature_bc_matrix/
├── ref/                       # NOT tracked — Azimuth / mapping references
└── results/                   # DE tables and figures
    └── mono_deg_mEVvssEV.csv
```

## Data

Raw matrices and references are **not** included in the repo (size / access).
Place `cellranger` outputs under `data/post_aggr/<sample>/` and point the
`data_dir` variable at the top of `mEV-PBMC.Rmd` to that location. The Azimuth
PBMC reference (`pbmcref`) is downloaded at runtime by `RunAzimuth`.

## Dependencies

R (≥ 4.2) with:

- **Seurat** (core analysis), **SeuratData**, **Azimuth**
- **Matrix**, **dplyr**, **tidyverse**, **ggplot2**, **ggridges**, **patchwork**, **ggrepel**
- **DESeq2** (via Seurat `FindMarkers`, `test.use = "DESeq2"`)
- **clusterProfiler**, **org.Hs.eg.db**, **AnnotationDbi** (GO enrichment)
- **monocle3**

Install Bioconductor packages with `BiocManager::install(...)`.

## Reproduction

1. Clone the repo and place data under `data/post_aggr/`.
2. Open `mEV-PBMC.Rmd` in RStudio and set `data_dir`.
3. Knit, or run chunks sequentially. DE tables are written to `results/`.

## Notes on this version

This notebook is a cleaned-up version of the original exploratory analysis.
Functional changes from the raw notebook: the repeated per-sample QC and
`FeaturePlot` chunks were collapsed into loops; the HTO assay block (previously
outside a runnable chunk) was made runnable; hard-coded absolute paths were
replaced with a relative `data/` convention; and the `enrichGO` call was
corrected to use `OrgDb = org.Hs.eg.db` with `keyType = "SYMBOL"`. The analysis
logic, thresholds, and annotations are unchanged.
