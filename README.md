# WhelanMultimodal
This repository contains R scripts used to generate figures for the snRNA/snATAC multimodal analysis of murine esophagus undergoing recurrent allergen challenge and/or carcinogenesis, published by Lab Whelan. Please download the 10X FASTQ files or the processed .h5/.tsv/.tbi files processed by Cell Ranger ARC 2.0.1 from the GEO accession [GSE354601]. Alternatively, please use the Seurat object RDS file (GSE354601_full.rds) located under the same GEO accession number. Seurat object was generated following Seurat's WNN analysis of 10x Multiome, RNA + ATAC vignette accessible in the satijalab seurat repository (vignettes/weighted_nearest_neighbor_analysis.Rmd).

The core dependencies for these scripts:

1. **R 4.5.0** https://www.r-project.org/

2. **Seurat 5.3.0** https://satijalab.org/seurat/

3. **Monocle 3** https://cole-trapnell-lab.github.io/monocle3/

All necessary pipelines to generate each figure are included in their respective scripts. Some supplementary scripts require pipelines delineated in the main figure scripts.
