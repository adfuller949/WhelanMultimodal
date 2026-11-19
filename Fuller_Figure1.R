library(tidyverse)
library(reshape2)
library(Seurat)
library(Signac)
library(clustree)
library(ggrepel)
library(ggpubr)

### Raw data processed with 10X Genomics CellRanger ARC
### Preprocessing performed following Seurat's WNN analysis of 10x Multiome, RNA + ATAC (vignettes/weighted_nearest_neighbor_analysis.Rmd)

# Figure 1F ---------------------------------------------------------------

full <- readRDS("GSE354601_full.rds")

# RNA QC
DefaultAssay(full) <- 'RNA'
full <- subset(full, subset = nFeature_RNA < 4000 & nFeature_RNA > 200)
full <- subset(full, subset = nCount_RNA < 15000 & nCount_RNA > 500)
full <- subset(full, subset = percent.mt < 30)

# ATAC QC
DefaultAssay(full) <- 'ATAC'
full <- subset(full, subset = nCount_ATAC < 50000 & nCount_ATAC > 1000)
full <- subset(full, subset = TSS.enrichment > 3)
full <- subset(full, subset = nucleosome_signal < 2)

### Final Clustering

# RNA SCT analysis
DefaultAssay(full) <- "RNA"
full <- SCTransform(full, vars.to.regress = 'percent.mt', verbose = FALSE)
full <- RunPCA(full, verbose = FALSE)
ElbowPlot(full, ndims = 50, reduction = "pca")
full <- RunUMAP(full, dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
DefaultAssay(full) <- "ATAC"
full <- RunTFIDF(full)
full <- FindTopFeatures(full, min.cutoff = 'q0')
full <- RunSVD(full)
full <- RunUMAP(full, reduction = 'lsi', dims = 2:30, 
                reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# Multimodal analysis
full <- FindMultiModalNeighbors(full, reduction.list = list("pca", "lsi"),
                                dims.list = list(1:30, 2:30))
full <- RunUMAP(full, nn.name = "weighted.nn", reduction.name = "wnn.umap",
                reduction.key = "wnnUMAP_")
full <- FindClusters(full, graph.name = "wsnn", algorithm = 3, verbose = FALSE,
                     resolution = seq(0.1,1,0.1))
clustree(full, prefix = "wsnn_res.")

table(Idents(full), full$orig.ident)
# Remove clusters with <15 cells and recluster

# Visualize
DimPlot(full, reduction = 'umap.rna', group.by = 'Condition')+
  scale_color_manual(values=c('skyblue','darkgreen','orangered','darkviolet'))

DimPlot(full, reduction = 'umap.atac', group.by = 'Condition')+
  scale_color_manual(values=c('skyblue','darkgreen','orangered','darkviolet'))

DimPlot(full, reduction = 'wnn.umap', group.by = 'Condition')+
  scale_color_manual(values=c('skyblue','darkgreen','orangered','darkviolet'))


# Figure 1G ---------------------------------------------------------------

# Normalize RNA assay
DefaultAssay(full) <- 'RNA'
full <- NormalizeData(full)
full <- FindVariableFeatures(full, selection.method = "vst", nfeatures = 2000)
full <- ScaleData(full, vars.to.regress = 'percent.mt')

# Identify cell types
DefaultAssay(full) <- 'RNA'
knownmarkers <- c('Cdh5', 'Emcn','Pecam1','Selp','Ptprb',
                  'C3', 'Pi16', 'Col3a1', 'Dcn','Gsn',
                  'Cd53','Icos','Itk','Ikzf1','Ptprc',
                  'Fabp5','Krt13','Col17a1','Krt14','Krt5')
DotPlot(object = full, features = knownmarkers, cols = c('navy','red'))+
  coord_flip() + RotatedAxis()

# Rename clusters
Idents(full) <- "wsnn_res.0.4"
full@meta.data$BroadCellTypes <- Idents(full)
Idents(full) <- "BroadCellTypes"
full <- RenameIdents(full,
                     '6'='Epithelial', '14'='Epithelial', '8'='Epithelial',
                     '10'='Epithelial', '1'='Epithelial', '11'='Epithelial',
                     '18'='Epithelial', '5'='Epithelial', '2'='Epithelial',
                     '4'='Epithelial', '0'='Epithelial', '9'='Epithelial',
                     '3'='Epithelial', '15'='Epithelial',
                     '19'='Fibroblast', '7'='Fibroblast',
                     '12'='Endothelial',
                     '13'='Immune', '16'='Immune', '17'='Immune')
full@meta.data$BroadCellTypes <- Idents(full)

# Visualize
levels(full) <- c('Epithelial','Immune','Fibroblast','Endothelial')
knownmarkers <- c('Cdh5', 'Emcn','Pecam1','Selp','Ptprb',
                  'C3', 'Pi16', 'Col3a1', 'Dcn','Gsn',
                  'Cd53','Icos','Itk','Ikzf1','Ptprc',
                  'Fabp5','Krt13','Col17a1','Krt14','Krt5')
DotPlot(object = full, features = knownmarkers, cols = c('navy','red'))+
  coord_flip() + RotatedAxis()



# Figure 1H ---------------------------------------------------------------

# Visualize
full$BroadCellTypes <- factor(full$BroadCellTypes, levels = c('Epithelial','Immune',
                                                              'Fibroblast','Endothelial'))

mycols <- c('lightgreen','deepskyblue','orange','plum')
DimPlot(full, reduction = 'wnn.umap', group.by = 'BroadCellTypes', cols = mycols)


# Figure 1I ---------------------------------------------------------------

# Visualize
data <- full %>%
  FetchData(vars = c("BroadCellTypes", "Condition")) %>%
  group_by(Condition, BroadCellTypes) %>%
  tally() %>%
  group_by(Condition) %>%
  mutate(percent = n / sum(n) * 100)

data$BroadCellTypes <- factor(
  data$BroadCellTypes,
  levels = rev(levels(data$BroadCellTypes))
)

ggplot(data, aes(x = Condition, y = percent, fill = BroadCellTypes)) + 
  geom_bar(stat = "identity") + 
  labs(x = "Condition", y = "Percent of Cells (%)") +
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c('plum','orange','deepskyblue','lightgreen')) +
  scale_x_discrete(limits = rev) +
  coord_flip()
