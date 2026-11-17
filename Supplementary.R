library(tidyverse)
library(reshape2)
library(Seurat)
library(Signac)
library(clustree)
library(ggrepel)
library(ggpubr)
library(monocle3)
library(org.Mm.eg.db)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(SeuratWrappers)
library(monocle3)
library(babelgene)

# Supp Fig 2 --------------------------------------------------------------

# Extract metadata
meta <- full@meta.data %>%
  select(orig.ident, Condition, BroadCellTypes)

# Create binary cell type
meta$CellGroup <- ifelse(
  meta$BroadCellTypes == "Epithelial",
  "Epithelial",
  "Other"
)

# Extract metadata
meta$CellGroup <- factor(meta$CellGroup,
                         levels = c("Other", "Epithelial"))

meta$Condition <- factor(
  meta$Condition,
  levels = c("Vehicle", "rOVA", "4NQO", "rOVA+4NQO")
)

# Differential proportion test
prop <- propeller(
  clusters = meta$CellGroup,
  sample = meta$orig.ident,
  group = meta$Condition,
  transform = "asin"
)

prop


# Pairwise

pairwise_propeller <- function(meta, cond1, cond2) {
  
  meta_sub <- subset(meta, Condition %in% c(cond1, cond2))
  meta_sub$Condition <- droplevels(meta_sub$Condition)
  
  propeller(
    clusters = meta_sub$CellGroup,
    sample = meta_sub$orig.ident,
    group = meta_sub$Condition,
    transform = "asin"
  )
}

Veh_vs_rOVA <- pairwise_propeller(meta, "Vehicle", "rOVA")
Veh_vs_4NQO <- pairwise_propeller(meta, "Vehicle", "4NQO")
Veh_vs_Combo <- pairwise_propeller(meta, "Vehicle", "rOVA+4NQO")
rOVA_vs_4NQO <- pairwise_propeller(meta, "rOVA", "4NQO")
rOVA_vs_Combo <- pairwise_propeller(meta, "rOVA", "rOVA+4NQO")
x4NQO_vs_Combo <- pairwise_propeller(meta, "4NQO", "rOVA+4NQO")

# Adjusted pairwise

pairwise_results <- data.frame(
  Comparison = c(
    "Vehicle_vs_rOVA",
    "Vehicle_vs_4NQO",
    "Vehicle_vs_rOVA+4NQO",
    "rOVA_vs_4NQO",
    "rOVA_vs_rOVA+4NQO",
    "x4NQO_vs_rOVA+4NQO"
  ),
  p.value = c(
    Veh_vs_rOVA["Epithelial", "P.Value"],
    Veh_vs_4NQO["Epithelial", "P.Value"],
    Veh_vs_Combo["Epithelial", "P.Value"],
    rOVA_vs_4NQO["Epithelial", "P.Value"],
    rOVA_vs_Combo["Epithelial", "P.Value"],
    x4NQO_vs_Combo["Epithelial", "P.Value"]
  )
)

pairwise_results$FDR <- p.adjust(pairwise_results$p.value, method = "BH")

pairwise_results



# Visualize

pairwise_results$Comparison <- factor(
  pairwise_results$Comparison,
  levels = pairwise_results$Comparison
)

ggplot(pairwise_results,
       aes(x = Comparison,
           y = -log10(FDR),
           fill = FDR < 0.05)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             color = "red") +
  scale_x_discrete(
    labels = function(x) gsub("_vs_", " v ", x)
  ) +
  scale_fill_manual(values = c("grey70", "#2C7BB6"), guide = "none") +
  labs(
    x = NULL,
    y = expression(-log[10](FDR))
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Change in epithelial cell population")

### Repeat for each cell type



# Supp Fig 3 --------------------------------------------------------------

### A

clustree(full.epi.no4NQO, prefix = "wsnn_res.")

### B

DimPlot(full.epi.no4NQO, reduction = 'umap.rna', group.by = 'Condition')+
  scale_color_manual(values=c('skyblue','darkgreen'))
DimPlot(full.epi.no4NQO, reduction = 'umap.atac', group.by = 'Condition')+
  scale_color_manual(values=c('skyblue','darkgreen'))
DimPlot(full.epi.no4NQO, reduction = 'wnn.umap', group.by = 'Condition')+
  scale_color_manual(values=c('skyblue','darkgreen'))

### C

Idents(full.epi.no4NQO) <- 'predicted.id'
cpt <- prop.table(table(full.epi.no4NQO$predicted.id,
                        full.epi.no4NQO$orig.ident), margin = 2)
cpt <- reshape2::melt(cpt)
seuratmd <- full.epi.no4NQO@meta.data
Conditions <- c()
for(x in cpt$Var2){
  Conditions <- append(Conditions,seuratmd[seuratmd$orig.ident==x,"Condition"][1])
}
cpt$Condition <- Conditions
colnames(cpt) <- c('predicted.id','Sample','Proportion','Condition')
cpt$Condition <- factor(cpt$Condition, levels = c('Vehicle','rOVA'))

ggpubr::ggboxplot(cpt, x = "Condition", y = "Proportion", color = "Condition", 
                  add = 'jitter', palette = c('skyblue','darkgreen'),
                  repel = T, ylim = c(0,1)) + facet_wrap(~predicted.id) +
  stat_compare_means(method = 'wilcox', comparisons = list(c('Vehicle','rOVA')),
                     label = 'p.format',
                     label.y = c(0.9), hide.ns = F)


# Supp Fig 4 --------------------------------------------------------------

### A

VlnPlot(full.epi.no4NQO, features = c("dim_Krt131","dim_Top2a1"),
        group.by = "predicted.id", pt.size = 0, 
        cols = c('#E38900','#00C1AD','#DC71FA'))

# Statistics
df <- data.frame(
  Var = full.epi.no4NQO$dim_Krt131,
  predicted.id = full.epi.no4NQO$predicted.id,
  Sample = full.epi.no4NQO$orig.ident
) %>%
  group_by(Sample, predicted.id) %>%
  summarise(
    MeanVar = mean(Var),
    .groups = "drop"
  )

results <- pairwise.wilcox.test(
  df$MeanVar,
  df$predicted.id,
  p.adjust.method = "BH"
)

### Repeat for dim_Top2a1


### B

VlnPlot(full.epi.no4NQO, features = c("dim_Krt131","dim_Top2a1"),
        group.by = "predicted.id", split.by = 'Condition', pt.size = 0, 
        cols = c("skyblue", "darkgreen"))

# Statistics
df <- data.frame(
  Var = full.epi.no4NQO$dim_Krt131,
  Condition = full.epi.no4NQO$Condition,
  predicted.id = full.epi.no4NQO$predicted.id,
  Sample = full.epi.no4NQO$orig.ident
) %>%
  group_by(Sample, Condition, predicted.id) %>%
  summarise(
    MeanVar = mean(Var),
    .groups = "drop"
  )

results <- pseudo %>%
  filter(predicted.id %in% c("Basal", "Suprabasal", "Superficial")) %>%
  group_by(Condition) %>%
  pairwise_wilcox_test(
    MeanVar ~ predicted.id,
    p.adjust.method = "BH"
  ) %>%
  as.data.frame()

results <- results %>%
  rowwise() %>%
  mutate(
    median1 = median(
      pseudo$MeanVar[
        pseudo$Condition == Condition &
          pseudo$predicted.id == group1
      ],
      na.rm = TRUE
    ),
    median2 = median(
      pseudo$MeanVar[
        pseudo$Condition == Condition &
          pseudo$predicted.id == group2
      ],
      na.rm = TRUE
    ),
    direction = case_when(
      median1 > median2 ~ paste(group1, "higher"),
      median2 > median1 ~ paste(group2, "higher"),
      TRUE ~ "equal"
    )
  ) %>%
  ungroup()

### Repeat for dim_Top2a1

### C

DefaultAssay(full.epi.no4NQO) <- 'RNA'
full.epi.no4NQO <- CellCycleScoring(full.epi.no4NQO,
                                    s.features = str_to_title(cc.genes.updated.2019$s.genes),
                                    g2m.features =  str_to_title(cc.genes.updated.2019$g2m.genes))

Idents(full.epi.no4NQO) <- 'Phase'
cpt <- prop.table(table(full.epi.no4NQO$Phase,
                        full.epi.no4NQO$orig.ident), margin = 2)
cpt <- melt(cpt)
seuratmd <- full.epi.no4NQO@meta.data
Conditions <- c()
for(x in cpt$Var2){
  Conditions <- append(Conditions,seuratmd[seuratmd$orig.ident==x,"Condition"][1])
}
cpt$Condition <- Conditions
colnames(cpt) <- c('Phase','Sample','Proportion','Condition')
cpt$Condition <- factor(cpt$Condition, levels = c('Vehicle','rOVA'))

ggpubr::ggboxplot(cpt, x = "Condition", y = "Proportion", color = "Condition", 
                  add = 'jitter', palette = c('skyblue','darkgreen'),
                  repel = T, ylim = c(0,1)) + facet_wrap(~Phase) +
  stat_compare_means(method = 'wilcox', comparisons = list(c('Vehicle','rOVA')),
                     label = 'p.format',
                     label.y = c(0.7), hide.ns = F)


# Supp Fig 5 --------------------------------------------------------------

### A

plot_cells(m3, color_cells_by = 'pseudotime',label_leaves = F,
           label_branch_points = F, label_roots = T, 
           show_trajectory_graph = T)

### B

plot_cells(m3, color_cells_by = 'Phase',label_leaves = F,
           label_branch_points = F, label_roots = T, 
           show_trajectory_graph = F)

### C

plot_cells(m3, color_cells_by = 'CellStates',label_leaves = F,
           label_branch_points = F, label_roots = T, 
           show_trajectory_graph = F) + facet_wrap(~Condition, ncol = 1)

### D

VlnPlot(full.epi.no4NQO, features = "mod261",
        group.by = "predicted.id", split.by = 'Condition', pt.size = 0, 
        cols = c("skyblue", "darkgreen"))

# Statistics
df <- data.frame(
  Var = full.epi.no4NQO$mod261,
  Condition = full.epi.no4NQO$Condition,
  predicted.id = full.epi.no4NQO$predicted.id,
  Sample = full.epi.no4NQO$orig.ident
) %>%
  group_by(Sample, Condition, predicted.id) %>%
  summarise(
    MeanVar = mean(Var),
    .groups = "drop"
  )

results <- pseudo %>%
  filter(predicted.id %in% c("Basal", "Suprabasal", "Superficial")) %>%
  group_by(Condition) %>%
  pairwise_wilcox_test(
    MeanVar ~ predicted.id,
    p.adjust.method = "BH"
  ) %>%
  as.data.frame()

results <- results %>%
  rowwise() %>%
  mutate(
    median1 = median(
      pseudo$MeanVar[
        pseudo$Condition == Condition &
          pseudo$predicted.id == group1
      ],
      na.rm = TRUE
    ),
    median2 = median(
      pseudo$MeanVar[
        pseudo$Condition == Condition &
          pseudo$predicted.id == group2
      ],
      na.rm = TRUE
    ),
    direction = case_when(
      median1 > median2 ~ paste(group1, "higher"),
      median2 > median1 ~ paste(group2, "higher"),
      TRUE ~ "equal"
    )
  ) %>%
  ungroup()


# Supp Fig 6 --------------------------------------------------------------

### A

# Run graph test
res <- graph_test(m3, neighbor_graph="principal_graph", cores=4)
ids <- row.names(subset(res, q_value < 0.05))

# Find gene modules by condition
gene_module_df <- find_gene_modules(m3[ids,], resolution=1e-2,
                                    k = 50)
cell_group_df <- tibble::tibble(cell=row.names(colData(m3)), 
                                cell_group=colData(m3)$Condition)

# Visualize heat map
agg_mat <- aggregate_gene_expression(m3, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
p <- pheatmap::pheatmap(agg_mat,
                        scale="column", clustering_method="ward.D2", 
                        fontsize_row = 7)

### B

# Find gene modules by cell state
cell_group_df <- tibble::tibble(cell=row.names(colData(m3)), 
                                cell_group=colData(m3)$CellStates)

# Visualize heat map
agg_mat <- aggregate_gene_expression(m3, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
p <- pheatmap::pheatmap(agg_mat,
                        scale="column", clustering_method="ward.D2", 
                        fontsize_row = 7)

### C

VlnPlot(epi.human, features = 'mod261', pt.size = 0, group.by = 'condition')

# Statistics
df <- data.frame(
  Var = epi.human$mod261,
  condition = epi.human$condition,
  Sample = epi.human$orig.ident
) %>%
  group_by(Sample, condition) %>%
  summarise(
    MeanVar = mean(Var),
    .groups = "drop"
  )

results <- pairwise.wilcox.test(
  df$MeanVar,
  df$condition,
  p.adjust.method = "BH"
)

### D

VlnPlot(full.epi.no4NQO, features = "smad31",
        group.by = "predicted.id", split.by = 'Condition', pt.size = 0, 
        cols = c("skyblue", "darkgreen"))

# Statistics
df <- data.frame(
  Var = full.epi.no4NQO$smad31,
  Condition = full.epi.no4NQO$Condition,
  predicted.id = full.epi.no4NQO$predicted.id,
  Sample = full.epi.no4NQO$orig.ident
) %>%
  group_by(Sample, Condition, predicted.id) %>%
  summarise(
    MeanVar = mean(Var),
    .groups = "drop"
  )

results <- pseudo %>%
  filter(predicted.id %in% c("Basal", "Suprabasal", "Superficial")) %>%
  group_by(Condition) %>%
  pairwise_wilcox_test(
    MeanVar ~ predicted.id,
    p.adjust.method = "BH"
  ) %>%
  as.data.frame()

results <- results %>%
  rowwise() %>%
  mutate(
    median1 = median(
      pseudo$MeanVar[
        pseudo$Condition == Condition &
          pseudo$predicted.id == group1
      ],
      na.rm = TRUE
    ),
    median2 = median(
      pseudo$MeanVar[
        pseudo$Condition == Condition &
          pseudo$predicted.id == group2
      ],
      na.rm = TRUE
    ),
    direction = case_when(
      median1 > median2 ~ paste(group1, "higher"),
      median2 > median1 ~ paste(group2, "higher"),
      TRUE ~ "equal"
    )
  ) %>%
  ungroup()


# Supp Fig 7 --------------------------------------------------------------

### A

clustree(full.epi.4NQO.ds, prefix = "wsnn_res.")

### B

DimPlot(full.epi.4NQO.ds, reduction = 'umap.rna', group.by = 'Condition')+
  scale_color_manual(values=c('orangered','darkviolet'))
DimPlot(full.epi.4NQO.ds, reduction = 'umap.atac', group.by = 'Condition')+
  scale_color_manual(values=c('orangered','darkviolet'))
DimPlot(full.epi.4NQO.ds, reduction = 'wnn.umap', group.by = 'Condition')+
  scale_color_manual(values=c('orangered','darkviolet'))

### D

DefaultAssay(full.epi.4NQO.ds) <- 'RNA'
full.epi.4NQO.ds <- CellCycleScoring(full.epi.4NQO.ds,
                                    s.features = str_to_title(cc.genes.updated.2019$s.genes),
                                    g2m.features =  str_to_title(cc.genes.updated.2019$g2m.genes))

Idents(full.epi.4NQO.ds) <- 'Phase'
cpt <- prop.table(table(full.epi.4NQO.ds$Phase,
                        full.epi.4NQO.ds$orig.ident), margin = 2)
cpt <- reshape2::melt(cpt)
seuratmd <- full.epi.4NQO.ds@meta.data
Conditions <- c()
for(x in cpt$Var2){
  Conditions <- append(Conditions,seuratmd[seuratmd$orig.ident==x,"Condition"][1])
}
cpt$Condition <- Conditions
colnames(cpt) <- c('Phase','Sample','Proportion','Condition')
cpt$Condition <- factor(cpt$Condition, levels = c('4NQO','rOVA+4NQO'))

ggpubr::ggboxplot(cpt, x = "Condition", y = "Proportion", color = "Condition", 
                  add = 'jitter', palette = c('orangered','darkviolet'),
                  repel = T, ylim = c(0,1)) + facet_wrap(~Phase) +
  stat_compare_means(method = 'wilcox', comparisons = list(c('4NQO','rOVA+4NQO')),
                     label = 'p.format',
                     label.y = c(0.9), hide.ns = F)
### E

Idents(full.epi.4NQO.ds) <- 'predicted.id'
cpt <- prop.table(table(full.epi.4NQO.ds$predicted.id,
                        full.epi.4NQO.ds$orig.ident), margin = 2)
cpt <- reshape2::melt(cpt)
seuratmd <- full.epi.4NQO.ds@meta.data
Conditions <- c()
for(x in cpt$Var2){
  Conditions <- append(Conditions,seuratmd[seuratmd$orig.ident==x,"Condition"][1])
}
cpt$Condition <- Conditions
colnames(cpt) <- c('predicted.id','Sample','Proportion','Condition')
cpt$Condition <- factor(cpt$Condition, levels = c('4NQO','rOVA+4NQO'))

ggpubr::ggboxplot(cpt, x = "Condition", y = "Proportion", color = "Condition", 
                  add = 'jitter', palette = c('orangered','darkviolet'),
                  repel = T, ylim = c(0,1.05)) + facet_wrap(~predicted.id) +
  stat_compare_means(method = 'wilcox', comparisons = list(c('4NQO','rOVA+4NQO')),
                     label = 'p.format',
                     label.y = c(0.95), hide.ns = F)



# Supp Fig 8 --------------------------------------------------------------

### A

# Convert Seurat object into Monocle3 object
cds <- as.cell_data_set(full.epi.4NQO.ds)
cds@int_colData$reducedDims$UMAP <- Embeddings(full.epi.4NQO.ds, reduction = "wnn.umap")
cds <- cluster_cells(cds)
cds <- learn_graph(cds, learn_graph_control = list(rann.k = 10))
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = 'pseudotime')

# Extract metadata
exmat <- GetAssayData(full.epi.4NQO.ds, slot = 'counts', assay = 'SCT')
mdat <- as.data.frame(cbind(full.epi.4NQO.ds$nCount_SCT,full.epi.4NQO.ds$nFeature_SCT,
                            full.epi.4NQO.ds$nCount_RNA,full.epi.4NQO.ds$nFeature_RNA,
                            full.epi.4NQO.ds$nCount_ATAC,full.epi.4NQO.ds$nFeature_ATAC,
                            full.epi.4NQO.ds$percent.mt,full.epi.4NQO.ds$Sex,
                            full.epi.4NQO.ds$Phase, full.epi.4NQO.ds$Condition,
                            full.epi.4NQO.ds$CellStates, full.epi.4NQO.ds$orig.ident,
                            full.epi.4NQO.ds$predicted.id))
colnames(mdat) <- c('nCount_SCT','nFeature_SCT','nCount_RNA','nFeature_RNA','nCount_ATAC',
                    'nFeature_ATAC','percent.mt','Sex','Phase','Condition','CellStates',
                    'orig.ident', 'predicted.id')
DefaultAssay(full.epi.4NQO.ds) <- 'SCT'
genedf <- cbind(rownames(full.epi.4NQO.ds),rownames(full.epi.4NQO.ds),rownames(full.epi.4NQO.ds))
rownames(genedf) <- genedf[,1]
genedf <- genedf[,-1]
colnames(genedf) <- c('gene_short_name','gene_name')

# Generate cell data set plot
m3_4NQO <- new_cell_data_set(exmat,cell_metadata = mdat,
                        gene_metadata = genedf)
m3_4NQO <- preprocess_cds(m3_4NQO, num_dim = 100)
m3_4NQO <- align_cds(m3_4NQO, alignment_group = "orig.ident")
m3_4NQO <- reduce_dimension(m3_4NQO, umap.n_neighbors = 50)
m3_4NQO <- cluster_cells(m3_4NQO)
m3_4NQO <- learn_graph(m3_4NQO, use_partition = FALSE)

# Add metadata
pData(m3_4NQO)$CellStates <- full.epi.4NQO.ds$CellStates
pData(m3_4NQO)$Phase <- full.epi.4NQO.ds$Phase
pData(m3_4NQO)$Condition <- full.epi.4NQO.ds$Condition
pData(m3_4NQO)$predicted.id <- full.epi.4NQO.ds$predicted.id

# Visualize cell cycle to identify pseudotime origin point
plot_cells(m3_4NQO, color_cells_by = 'Phase', 
           label_groups_by_cluster = F,label_leaves = F,
           label_branch_points = F, label_cell_groups = F,
           trajectory_graph_color = 'grey45', 
           show_trajectory_graph = F)

# Order cells at pseudotime origin point
m3_4NQO <- order_cells(m3_4NQO)

# Visualize
plot_cells(m3_4NQO, color_cells_by = 'pseudotime',label_leaves = F,
           label_branch_points = F, label_roots = T, 
           show_trajectory_graph = T)

### B

plot_cells(m3_4NQO, color_cells_by = 'Phase',label_leaves = F,
           label_branch_points = F, label_roots = T, 
           show_trajectory_graph = F)

### C

plot_cells(m3_4NQO, color_cells_by = 'CellStates',label_leaves = F,
           label_branch_points = F, label_roots = T, 
           show_trajectory_graph = F) + facet_wrap(~Condition, ncol = 1)

### D

# Run graph test
res <- graph_test(m3_4NQO, neighbor_graph="principal_graph", cores=4)
ids <- row.names(subset(res, q_value < 0.05))

# Find gene modules
gene_module_df <- find_gene_modules(m3_4NQO[ids,], resolution=1e-2,
                                    k = 50)
cell_group_df <- tibble::tibble(cell=row.names(colData(m3_4NQO)), 
                                cell_group=colData(m3_4NQO)$Condition)

# Visualize heat map
agg_mat <- aggregate_gene_expression(m3_4NQO, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
p <- pheatmap::pheatmap(agg_mat,
                        scale="column", clustering_method="ward.D2", 
                        fontsize_row = 7)

### E

# Extract modules of interest and add to Seurat object
mod8 <- gene_module_df[gene_module_df$module=='8',]
mod8_genes <- mod8[[1]]

full.epi.4NQO.ds <- AddModuleScore(
  full.epi.4NQO.ds,
  features = list(mod8_genes),
  name = "mod8"
)

# Visualize
VlnPlot(full.epi.4NQO.ds, features = 'mod81', pt.size = 0) +
  scale_fill_manual(values=c('orangered','darkviolet'))

# Statistics
df <- data.frame(
  Mod = full.epi.4NQO.ds$mod81,
  Condition = full.epi.4NQO.ds$Condition,
  Sample = full.epi.4NQO.ds$orig.ident
) %>%
  group_by(Sample, Condition) %>%
  summarise(
    MeanMod = mean(Mod),
    .groups = "drop"
  )

results <- df %>%
  summarise(
    p.value = wilcox.test(
      MeanMod ~ Condition
    )$p.value
  )

### Repeat for modules 9, 15, and 16


# Supp Fig 9 --------------------------------------------------------------

### A

VlnPlot(full.epi.4NQO.ds, features = c("dim_Sox61","dim_Krt171"),
        group.by = "predicted.id", pt.size = 0, 
        cols = c('#E38900','#00C1AD','#DC71FA'))

# Statistics
df <- data.frame(
  Var = full.epi.4NQO.ds$dim_Sox61,
  predicted.id = full.epi.4NQO.ds$predicted.id,
  Sample = full.epi.4NQO.ds$orig.ident
) %>%
  group_by(Sample, predicted.id) %>%
  summarise(
    MeanVar = mean(Var),
    .groups = "drop"
  )

results <- pairwise.wilcox.test(
  df$MeanVar,
  df$predicted.id,
  p.adjust.method = "BH"
)

### Repeat for dim_Krt171

### B

VlnPlot(full.epi.4NQO.ds, features = c("dim_Sox61","dim_Krt171"),
        group.by = "CellStates", pt.size = 0,
        cols = c('gold','cornflowerblue','mediumpurple1','burlywood','aquamarine',
                 'mediumorchid1','chartreuse3','dodgerblue','hotpink1','salmon'))

### C

Idents(full.epi.4NQO.ds) <- 'Condition'
VlnPlot(full.epi.4NQO.ds, features = c('dim_Sox61','dim_Krt171'), pt.size = 0) +
  scale_fill_manual(values=c('orangered','darkviolet'))

# Statistics
df <- data.frame(
  Dim = full.epi.4NQO.ds$dim_Sox61,
  Condition = full.epi.4NQO.ds$Condition,
  Sample = full.epi.4NQO.ds$orig.ident
) %>%
  group_by(Sample, Condition) %>%
  summarise(
    MeanDim = mean(Dim),
    .groups = "drop"
  )

results <- df %>%
  summarise(
    p.value = wilcox.test(
      MeanDim ~ Condition
    )$p.value
  )

### Repeat for dim_Krt171

# Supp Fig 10 -------------------------------------------------------------

counts <- GetAssayData(full.epi.4NQO.ds, assay = "ATAC", layer = "counts")
granges <- StringToGRanges(rownames(counts), sep = c("-", "-"))

full.epi.4NQO.ds.ca <- full.epi.4NQO.ds

full.epi.4NQO.ds.ca[["ATAC"]] <- CreateChromatinAssay(
  counts = counts,
  ranges = granges,
  genome = "mm10"
)

class(full.epi.4NQO.ds.ca[["ATAC"]])

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 10090, all_versions = FALSE)
)

DefaultAssay(full.epi.4NQO.ds.ca) <- 'ATAC'

full.epi.4NQO.ds.ca <- AddMotifs(
  object = full.epi.4NQO.ds.ca,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

full.epi.4NQO.ds.ca <- RegionStats(
  object = obj,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

markers <- FindMarkers(full.epi.4NQO.ds.ca,
                       ident.1 = 'rOVA',
                       ident.2 = 'Vehicle',
                       only.pos = T,
                       assay = 'ATAC')
sig_peaks <- rownames(markers[markers$p_val_adj < 0.05, ])

meta <- full.epi.4NQO.ds.ca[["ATAC"]]@meta.features
valid_peaks <- rownames(meta)[!is.na(meta$GC.percent)]

full.epi.4NQO.ds.ca <- subset(
  full.epi.4NQO.ds.ca,
  features = valid_peaks
)
sig_peaks <- intersect(sig_peaks, valid_peaks)

enriched_motifs <- FindMotifs(
  object = full.epi.4NQO.ds.ca,
  features = sig_peaks,
  background = 80000
)

sig_motifs <- enriched_motifs[
  enriched_motifs$p.adjust < 0.05, 
]
top_sig <- sig_motifs[order(sig_motifs$p.adjust), ][1:min(10, nrow(sig_motifs)), ]

MotifPlot(
  object = full.epi.4NQO.ds.ca,
  motifs = head(rownames(top_sig))
)


# Supp Fig 12 -------------------------------------------------------------

### A

DimPlot(full.immune, reduction = 'umap.rna', group.by = 'Condition')+
  scale_color_manual(values=c('skyblue','darkgreen','orangered','darkviolet'))
DimPlot(full.immune, reduction = 'umap.atac', group.by = 'Condition')+
  scale_color_manual(values=c('skyblue','darkgreen','orangered','darkviolet'))
DimPlot(full.immune, reduction = 'wnn.umap', group.by = 'Condition')+
  scale_color_manual(values=c('skyblue','darkgreen','orangered','darkviolet'))

### B

# RNA
DefaultAssay(full.immune) <- 'RNA'
Idents(full.immune) <- 'CellStates'
topmarkers <- FindAllMarkers(full.immune, assay = "RNA")
top100 <- topmarkers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC)

conditions_avex <- AggregateExpression(full.immune, return.seurat = T, assays = 'RNA')
conditions_avex <- ScaleData(conditions_avex, assay = 'RNA', features = top100$gene)
DoHeatmap(object = conditions_avex, features = top100$gene, assay = 'RNA', 
          group.colors = c('chartreuse3','dodgerblue','gold'),
          label = T, size = 3, angle = 30, draw.lines = F)+
  guides(colour=FALSE)+ 
  theme(axis.text.y = element_blank())

# ATAC
DefaultAssay(full.immune) <- 'ATAC'
Idents(full.immune) <- 'CellStates'
topmarkers <- FindAllMarkers(full.immune, assay = "ATAC")
top100 <- topmarkers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC)

conditions_avex <- AggregateExpression(full.immune, return.seurat = T, assays = 'ATAC')
conditions_avex <- ScaleData(conditions_avex, assay = 'ATAC', features = top100$gene)
DoHeatmap(object = conditions_avex, features = top100$gene, assay = 'ATAC', 
          group.colors = c('chartreuse3','dodgerblue','gold'),
          label = T, size = 3, angle = 30, draw.lines = F)+
  guides(colour=FALSE)+ 
  theme(axis.text.y = element_blank())

### C

Idents(full.immune) <- 'CellTypes'
cpt <- prop.table(table(full.immune$CellTypes,
                        full.immune$orig.ident), margin = 2)
cpt <- melt(cpt)
seuratmd <- full.immune@meta.data
Conditions <- c()
for(x in cpt$Var2){
  Conditions <- append(Conditions,seuratmd[seuratmd$orig.ident==x,"Condition"][1])
}
cpt$Condition <- Conditions
colnames(cpt) <- c('Cluster','Sample','Proportion','Condition')
cpt$Condition <- factor(cpt$Condition, levels = c('Vehicle','rOVA','4NQO','rOVA+4NQO'))

ggpubr::ggboxplot(cpt, x = "Condition", y = "Proportion", color = "Condition", 
                  add = 'jitter', palette = c('skyblue','darkgreen','orangered','darkviolet'),
                  repel = T, ylim = c(0,1.4)) + facet_wrap(~Cluster) +
  stat_compare_means(method = 'wilcox', comparisons = list(c('Vehicle','rOVA'),
                                                           c('rOVA', '4NQO'),
                                                           c('4NQO', 'rOVA+4NQO'),
                                                           c('Vehicle','4NQO'),
                                                           c('Vehicle', 'rOVA+4NQO'),
                                                           c('rOVA','rOVA+4NQO')),
                     label = 'p.format',
                     label.y = c(0.75, 0.8, 0.85, 0.95, 1.1, 1.25), hide.ns = F)+
  rremove('x.text')+
  facet_wrap(~Cluster)