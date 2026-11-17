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

# Figure 3A ---------------------------------------------------------------

# DimHeatMap
full.epi.no4NQO <- ProjectDim(full.epi.no4NQO, reduction = 'wnn.umap', assay = "RNA")
DimHeatmap(object = full.epi.no4NQO, reduction = 'wnn.umap', projected = T, 
           dims = 1, cells = 1000, balanced = T, fast = F, nfeatures = 100)

# To reorder
plot <- DimHeatmap(object = full.epi.no4NQO, reduction = 'wnn.umap', projected = T, 
                   dims = 1, cells = 1000, balanced = T, fast = F, nfeatures = 100)
current_genes <- ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range
flipped_genes <- rev(current_genes)
plot + scale_y_discrete(limits = flipped_genes)

# Generate gene programs
dim_Krt13 <- c(
  "Krt13","Krt4","Mt4","Fabp5","Tgm3","Krtdap","Ly6d","Cyp2f2","Acsl1","Gsta4",
  "Lgals7","Pla2g4f","Adh7","Ablim1","Notch3","Calml3","Dmkn","Hsd17b2","Lgals3","Sbsn",
  "Gm94","Dlg2","4833423E24Rik","Dsc2","Ces1h","Dsg1b","Sptssb","Ctsd","Kdm5b","Scel",
  "Grhl1","S100a14","Pla2g4e","Rpl11","Eef1a1","Rab25","Evpl","Nectin4","Pkp1","Lipm",
  "Cyp4f39","Lypd3","Pycard","Spink5","Rdh12","Hivep3","Abca12","Ppl","Rps2","Them5"
)
write.csv(dim_Krt13, file = "full.epi.no4NQO_dim_Krt13.csv")

dim_Top2a <- c(
  "Cep128","Hmgb2","Birc5","Ndc80","Cdk1","Bub1b","Dlgap5","Sgo2a","Ccdc18","Cep55",
  "Cip2a","Smc2","Kif14","Ttk","Ccnf","Racgap1","Ckap2l","Gas2l3","Dnah8","Lockd",
  "Kif18a","Hmmr","Cdca2","Efcab11","Cenpa","Parpbp","Incenp","Cenpe","Depdc1a","Kif15",
  "Cenpf","Cenpp","Smc4","Kif4","Mki67","Nusap1","Kif20b","Ckap2","Tpx2","Aspm",
  "Mis18bp1","Kif23","Cit","Ect2","Prc1","Kif11","Knl1","Diaph3","Anln","Top2a"
)
write.csv(dim_Top2a, file = "full.epi.no4NQO_dim_Top2a.csv")

# Add gene programs as modules
full.epi.no4NQO <- AddModuleScore(
  full.epi.no4NQO,
  features = list(dim_Krt13),
  name = "dim_Krt13"
)

full.epi.no4NQO <- AddModuleScore(
  full.epi.no4NQO,
  features = list(dim_Top2a),
  name = "dim_Top2a"
)

# Visualize
FeaturePlot(full.epi.no4NQO, reduction = "wnn.umap", features = c("dim_Krt131","dim_Top2a1"))



# Figure 3B ---------------------------------------------------------------

VlnPlot(full.epi.no4NQO, features = c("dim_Krt131","dim_Top2a1"),
        group.by = "predicted.id", split.by = 'Condition', pt.size = 0, 
        cols = c("skyblue", "darkgreen"))

# Statistics
df <- data.frame(
  Dim = full.epi.no4NQO$dim_Krt131,
  Condition = full.epi.no4NQO$Condition,
  predicted.id = full.epi.no4NQO$predicted.id,
  Sample = full.epi.no4NQO$orig.ident
) %>%
  group_by(Sample, Condition, predicted.id) %>%
  summarise(
    MeanDim = mean(Dim),
    .groups = "drop"
  )

results <- df %>%
  filter(predicted.id %in% c("Basal", "Suprabasal", "Superficial")) %>%
  group_by(predicted.id) %>%
  summarise(
    p.value = wilcox.test(
      MeanDim ~ Condition
    )$p.value
  )

results$FDR <- p.adjust(results$p.value, method = "BH")

### Repeat for dim_Top2a1


# Figure 3C ---------------------------------------------------------------

# Convert Seurat object into Monocle3 object
cds <- as.cell_data_set(full.epi.no4NQO)
cds@int_colData$reducedDims$UMAP <- Embeddings(full.epi.no4NQO, reduction = "wnn.umap")
cds <- cluster_cells(cds)
cds <- learn_graph(cds, learn_graph_control = list(rann.k = 10))
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = 'pseudotime')

# Extract metadata
exmat <- GetAssayData(full.epi.no4NQO, slot = 'counts', assay = 'SCT')
mdat <- as.data.frame(cbind(full.epi.no4NQO$nCount_SCT,full.epi.no4NQO$nFeature_SCT,
                            full.epi.no4NQO$nCount_RNA,full.epi.no4NQO$nFeature_RNA,
                            full.epi.no4NQO$nCount_ATAC,full.epi.no4NQO$nFeature_ATAC,
                            full.epi.no4NQO$percent.mt,full.epi.no4NQO$Sex,
                            full.epi.no4NQO$Phase, full.epi.no4NQO$Condition,
                            full.epi.no4NQO$CellStates, full.epi.no4NQO$orig.ident,
                            full.epi.no4NQO$predicted.id))
colnames(mdat) <- c('nCount_SCT','nFeature_SCT','nCount_RNA','nFeature_RNA','nCount_ATAC',
                    'nFeature_ATAC','percent.mt','Sex','Phase','Condition','CellStates',
                    'orig.ident', 'predicted.id')
DefaultAssay(full.epi.no4NQO) <- 'SCT'
genedf <- cbind(rownames(full.epi.no4NQO),rownames(full.epi.no4NQO),rownames(full.epi.no4NQO))
rownames(genedf) <- genedf[,1]
genedf <- genedf[,-1]
colnames(genedf) <- c('gene_short_name','gene_name')

# Generate cell data set plot
m3 <- new_cell_data_set(exmat,cell_metadata = mdat,
                        gene_metadata = genedf)
m3 <- preprocess_cds(m3, num_dim = 100)
m3 <- align_cds(m3, alignment_group = "orig.ident")
m3 <- reduce_dimension(m3, umap.n_neighbors = 50)
m3 <- cluster_cells(m3)
m3 <- learn_graph(m3, use_partition = FALSE)

# Add metadata
pData(m3)$CellStates <- full.epi.no4NQO$CellStates
pData(m3)$Phase <- full.epi.no4NQO$Phase
pData(m3)$Condition <- full.epi.no4NQO$Condition
pData(m3)$predicted.id <- full.epi.no4NQO$predicted.id

# Visualize cell cycle to identify pseudotime origin point
plot_cells(m3, color_cells_by = 'Phase', 
           label_groups_by_cluster = F,label_leaves = F,
           label_branch_points = F, label_cell_groups = F,
           trajectory_graph_color = 'grey45', 
           show_trajectory_graph = F)

# Order cells at pseudotime origin point
m3 <- order_cells(m3)

# Visualize
plot_cells(m3, color_cells_by = 'pseudotime',label_leaves = F,
           label_branch_points = F, label_roots = T, 
           show_trajectory_graph = F) + facet_wrap(~Condition, ncol = 1)


# Figure 3D ---------------------------------------------------------------

# Run graph test
res <- graph_test(m3, neighbor_graph="principal_graph", cores=4)
ids <- row.names(subset(res, q_value < 0.05))

# Find gene modules
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

# Extract module of interest and add to Seurat object
mod26 <- gene_module_df[gene_module_df$module=='26',]
mod26_genes <- mod26[[1]]

full.epi.no4NQO <- AddModuleScore(
  full.epi.no4NQO,
  features = list(mod26_genes),
  name = "mod26"
)

# Visualize
VlnPlot(full.epi.no4NQO, features = 'mod261', pt.size = 0) +
  scale_fill_manual(values=c('skyblue','darkgreen'))

# Statistics
df <- data.frame(
  Mod = full.epi.no4NQO$mod261,
  Condition = full.epi.no4NQO$Condition,
  Sample = full.epi.no4NQO$orig.ident
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



# Figure 3E ---------------------------------------------------------------

plot_cells(m3,
           genes=gene_module_df %>% filter(module %in% c('26')),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE, label_roots = F) + facet_wrap('Condition')


# Figure 3F ---------------------------------------------------------------

VlnPlot(full.epi.no4NQO, features = 'mod261', pt.size = 0,
        split.by = 'Condition', group.by = 'predicted.id') +
  scale_fill_manual(values=c('skyblue','darkgreen')) +
  labs(title = "Module 26")

# Statistics
df <- data.frame(
  Mod = full.epi.no4NQO$mod261,
  Condition = full.epi.no4NQO$Condition,
  predicted.id = full.epi.no4NQO$predicted.id,
  Sample = full.epi.no4NQO$orig.ident
) %>%
  group_by(Sample, Condition, predicted.id) %>%
  summarise(
    MeanMod = mean(Mod),
    .groups = "drop"
  )

results <- df %>%
  filter(predicted.id %in% c("Basal", "Suprabasal", "Superficial")) %>%
  group_by(predicted.id) %>%
  summarise(
    p.value = wilcox.test(
      MeanMod ~ Condition
    )$p.value
  )

results$FDR <- p.adjust(results$p.value, method = "BH")


# Figure 3G ---------------------------------------------------------------

### Imputation from GEO GSE201153

# Unzip GEO folder, rename for Read10X
untar(
  "C:/Users/GSE201153_RAW.tar",
  exdir = "C:/Users/GSE201153_RAW"
)

base_dir <- "C:/Users/GSE201153_RAW"

files <- list.files(base_dir, full.names = TRUE)

gsm_ids <- unique(sub("_.*", "", basename(files)))

for (gsm in gsm_ids) {
  
  gsm_dir <- file.path(base_dir, gsm)
  dir.create(gsm_dir, showWarnings = FALSE)
  
  gsm_files <- files[grepl(gsm, basename(files))]
  
  for (f in gsm_files) {
    
    if (grepl("barcodes.tsv.gz$", f)) {
      file.copy(f, file.path(gsm_dir, "barcodes.tsv.gz"), overwrite = TRUE)
      
    } else if (grepl("features.tsv.gz$", f)) {
      file.copy(f, file.path(gsm_dir, "features.tsv.gz"), overwrite = TRUE)
      
    } else if (grepl("matrix.mtx.gz$", f)) {
      file.copy(f, file.path(gsm_dir, "matrix.mtx.gz"), overwrite = TRUE)
    }
  }
}


# Run Read10x
human_list <- lapply(
  gsm_ids,
  function(gsm) Read10X(data.dir = file.path(base_dir, gsm))
)


names(human_list) <- gsm_ids


# Create individual Seurat objects
seurat_list <- lapply(
  names(human_list),
  function(gsm) {
    CreateSeuratObject(
      counts = human_list[[gsm]],
      project = gsm
    )
  }
)
names(seurat_list) <- names(human_list)

# Merge
combined <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list)
)

# Add metadata
combined$sample <- sub("_.*", "", colnames(combined))

sample_meta <- data.frame(
  sample = c(
    "GSM6051736", "GSM6051737", "GSM6051738", "GSM6051739", "GSM6051740",
    "GSM6051741", "GSM6051742", "GSM6051743", "GSM6051744", "GSM6051745"
  ),
  condition = c(
    "Active", "Remiss", "Normal", "Active", "Remiss",
    "Active", "Active", "Normal", "Active", "Remiss"
  ),
  run = c(
    "run2233", "run2242", "run2272", "run2295", "run2312",
    "run2326", "run2403", "run2422", "run2441", "run2442"
  ),
  stringsAsFactors = FALSE
)

combined$sample <- sub("_.*", "", colnames(combined))

meta_map <- sample_meta
rownames(meta_map) <- meta_map$sample

combined$condition <- meta_map[combined$sample, "condition"]
combined$run <- meta_map[combined$sample, "run"]


# QC for Seurat object
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")
VlnPlot(object = combined, features = c('nFeature_RNA','nCount_RNA','percent.mt'),
        pt.size = 0, group.by = 'orig.ident')
combined <- subset(combined, subset = nFeature_RNA < 4000 & nFeature_RNA > 200)
combined <- subset(combined, subset = nCount_RNA < 15000 & nCount_RNA > 500)
combined <- subset(combined, subset = percent.mt < 25)

# Cluster
combined <- SCTransform(combined, vars.to.regress = 'percent.mt', verbose = FALSE)
combined <- RunPCA(combined, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(combined, dims = 1:30, verbose = FALSE)
combined <- FindClusters(combined, verbose = FALSE,
                         resolution = seq(0.1,1,0.1))
clustree(combined, prefix = "SCT_snn_res.")

# Visualize
DimPlot(combined, group.by = 'SCT_snn_res.0.1', split.by = 'condition')
FeaturePlot(combined, features = c('CCL5','LYZ','TPSAB1','PLVAP','LUM')) # Non-epithelial
FeaturePlot(combined, features = c('KRT5','KRT13')) # Epithelial
VlnPlot(combined, features = c('MKI67','KRT5','KRT13'), group.by = 'SCT_snn_res.0.1')

# Remove non-epithelial cells and re-cluster
epi.human <- subset(combined, subset = SCT_snn_res.0.1 %in% c('2','4','7','9','10'),
                    invert = T)
epi.human <- SCTransform(epi.human, vars.to.regress = 'percent.mt', verbose = FALSE)
epi.human <- RunPCA(epi.human, verbose = FALSE)
epi.human <- RunUMAP(epi.human, reduction = "pca", dims = 1:20)
epi.human <- FindNeighbors(epi.human, dims = 1:30, verbose = FALSE)
epi.human <- FindClusters(epi.human, verbose = FALSE,
                          resolution = seq(0.1,1,0.1))
clustree(epi.human, prefix = "SCT_snn_res.")

# Normalize RNA
DefaultAssay(epi.human) <- 'RNA'
epi.human <- NormalizeData(epi.human)
epi.human <- FindVariableFeatures(epi.human, selection.method = "vst", nfeatures = 2000)
epi.human <- ScaleData(epi.human, vars.to.regress = 'percent.mt')

# Visualize
DimPlot(epi.human, group.by = 'SCT_snn_res.0.2', split.by = 'condition')

FeaturePlot(epi.human, features = c('KRT15','TOP2A','DSP','SERPINB3','TGM3'))
VlnPlot(epi.human, features = c('KRT15','TOP2A','DSP','SERPINB3','TGM3'), 
        group.by = 'SCT_snn_res.0.2', pt.size = 0)

# Rename clusters
DefaultAssay(epi.human) <- 'RNA'
Idents(epi.human) <- 'SCT_snn_res.0.2'
epi.human@meta.data$CellStates <- Idents(epi.human)
Idents(epi.human) <- "CellStates"
epi.human <- RenameIdents(epi.human,
                          '4'='Quiescent','8'='Quiescent',
                          '1'='Proliferative',
                          '0'='Intermediate 2','7'='Intermediate 2',
                          '2'='Intermediate 3',
                          '6'='Intermediate 1',
                          '3'='Differentiated','5'='Differentiated')
epi.human@meta.data$CellStates <- Idents(epi.human)

# Re-level
epi.human$condition <- factor(
  epi.human$condition, 
  levels = c("Normal", "Active", "Remiss")
)

# Visualize
DimPlot(epi.human, group.by = 'CellStates', split.by = 'condition')


# Figure 3H ---------------------------------------------------------------

# Convert pseudotime module data from mouse to human
mod26_genes <- mod26[[1]]
mod26_human <- orthologs(genes = mod26_genes, species = "mouse", human = FALSE)
mod26_human <- mod26_human[[1]]

# Add module data
epi.human <- AddModuleScore(
  epi.human,
  features = list(mod26_human),
  name = "mod26"
)

# Visualize
VlnPlot(epi.human, features = 'mod261', pt.size = 0, group.by = 'CellStates')+
  labs(title = "Module 26")

# Statistics
df <- FetchData(epi.human, vars = c(
  "mod261",
  "CellStates",
  "orig.ident"
))

df_sample <- df %>%
  group_by(orig.ident, CellStates) %>%
  summarise(
    value = mean(mod261),
    .groups = "drop"
  )

cellstates <- unique(df_sample$CellStates)
comparisons <- combn(cellstates, 2, simplify = FALSE)

results <- lapply(comparisons, function(comp) {
  
  x <- df_sample %>%
    filter(CellStates == comp[1]) %>%
    pull(value)
  
  y <- df_sample %>%
    filter(CellStates == comp[2]) %>%
    pull(value)
  
  test <- wilcox.test(x, y)
  
  data.frame(
    group1 = comp[1],
    group2 = comp[2],
    n_group1 = length(x),
    n_group2 = length(y),
    statistic = as.numeric(test$statistic),
    p_value = test$p.value
  )
}) %>%
  bind_rows()

results <- results %>%
  mutate(
    p_adjusted_BH = p.adjust(p_value, method = "BH")
  )


# Figure 3I ---------------------------------------------------------------

# Set genome
TxDb.Mm <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Find differentially accessible regions
DefaultAssay(full.epi.no4NQO) <- 'ATAC'
Idents(full.epi.no4NQO) <- 'Condition'
topmarkers <- FindMarkers(full.epi.no4NQO,
                          ident.1 = 'rOVA',
                          ident.2 = 'Vehicle',
                          only.pos = T, 
                          assay = "ATAC")
signif <- topmarkers %>% filter(p_val_adj < 0.05)

# Reformat
signif_spl <- strsplit(signif$gene, split = "-" )
gene_spl <- matrix(unlist(signif_spl),ncol=3,byrow=T)
colnames(gene_spl) <- c('chr', 'start', 'end')
gene_spl <- data.frame(gene_spl)
signif_rm <- signif[ -c(7) ]
signif_spl <- cbind(signif_rm,gene_spl)
signif_granges <- GRanges(signif_spl)

# Annotate peaks
peakAnno <- annotatePeak(signif_granges, tssRegion=c(-3000, 3000), TxDb=TxDb.Mm)
write.csv(peakAnno@anno, file = 'peakanno.csv')

# Convert Entrez geneid to gene name
peakAnno <- read.csv('peakanno.csv')
myGenes <- peakAnno$geneId
myGenes <- as.character(myGenes)
gene <- select(org.Mm.eg.db, keys=myGenes, columns=c("SYMBOL"),keytype="ENTREZID")
peakAnno$gene <- gene
write.csv(peakAnno, file = 'peakanno_final.csv')

# Reformat for pie chart
data <- read.csv("peakanno_final.csv", stringsAsFactors = FALSE)
data <- data %>%
  mutate(annotation = case_when(
    str_detect(annotation, regex("intron 1 of", ignore_case = TRUE)) ~ "1st Intron",
    str_detect(annotation, regex("exon 1 of", ignore_case = TRUE))   ~ "1st Exon",
    str_detect(annotation, regex("intron", ignore_case = TRUE)) ~ "Other Intron",
    str_detect(annotation, regex("exon", ignore_case = TRUE))   ~ "Other Exon",
    TRUE ~ annotation
  ))
unique(data$annotation)

data$annotation <- factor(
  data$annotation,
  levels = c(
    "Promoter (<=1kb)",
    "Promoter (1-2kb)",
    "Promoter (2-3kb)",
    "5' UTR",
    "3' UTR",
    "1st Exon",
    "Other Exon",
    "1st Intron",
    "Other Intron",
    "Downstream (<=300bp)",
    "Distal Intergenic"
  )
)

data$cluster <- factor(
  data$cluster,
  levels = c("Vehicle", "rOVA")
)

# Summarize data
summary_df <- data %>%
  group_by(cluster, annotation) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(percent = count / sum(count))

# Visualize
ggplot(summary_df, aes(x = "", y = percent, fill = annotation)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  facet_wrap(~ cluster) +
  theme_void()


# Figure 3J ---------------------------------------------------------------

counts <- GetAssayData(full.epi.no4NQO, assay = "ATAC", layer = "counts")
granges <- StringToGRanges(rownames(counts), sep = c("-", "-"))

full.epi.no4NQO.ca <- full.epi.no4NQO

full.epi.no4NQO.ca[["ATAC"]] <- CreateChromatinAssay(
  counts = counts,
  ranges = granges,
  genome = "mm10"
)

class(full.epi.no4NQO.ca[["ATAC"]])

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 10090, all_versions = FALSE)
)

DefaultAssay(full.epi.no4NQO.ca) <- 'ATAC'

full.epi.no4NQO.ca <- AddMotifs(
  object = full.epi.no4NQO.ca,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

full.epi.no4NQO.ca <- RegionStats(
  object = obj,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

markers <- FindMarkers(full.epi.no4NQO.ca,
                       ident.1 = 'rOVA',
                       ident.2 = 'Vehicle',
                       only.pos = T,
                       assay = 'ATAC')
sig_peaks <- rownames(markers[markers$p_val_adj < 0.05, ])

meta <- full.epi.no4NQO.ca[["ATAC"]]@meta.features
valid_peaks <- rownames(meta)[!is.na(meta$GC.percent)]

full.epi.no4NQO.ca <- subset(
  full.epi.no4NQO.ca,
  features = valid_peaks
)
sig_peaks <- intersect(sig_peaks, valid_peaks)

enriched_motifs <- FindMotifs(
  object = full.epi.no4NQO.ca,
  features = sig_peaks,
  background = 80000
)

sig_motifs <- enriched_motifs[
  enriched_motifs$p.adjust < 0.05, 
]

MotifPlot(
  object = full.epi.no4NQO.ca,
  motifs = head(rownames(sig_motifs))
)


# Figure 3K ---------------------------------------------------------------

# Create gene target module
smad3_genes <- c('Deptor','Gli1','Gli2','Ptprk','Spon1','Tbx21')

full.epi.no4NQO <- AddModuleScore(
  full.epi.no4NQO,
  features = list(smad3_genes),
  name = "smad3"
)

# Visualize
VlnPlot(full.epi.no4NQO, features = 'smad31', pt.size = 0, group.by = 'predicted.id', 
        split.by = 'Condition') +
  scale_fill_manual(values=c('skyblue','darkgreen')) +
  labs(title = "Smad3 Targets")

# Statistics
df <- data.frame(
  Smad = full.epi.no4NQO$smad31,
  Condition = full.epi.no4NQO$Condition,
  predicted.id = full.epi.no4NQO$predicted.id,
  Sample = full.epi.no4NQO$orig.ident
) %>%
  group_by(Sample, Condition, predicted.id) %>%
  summarise(
    MeanSmad = mean(Smad),
    .groups = "drop"
  )

results <- df %>%
  filter(predicted.id %in% c("Basal", "Suprabasal", "Superficial")) %>%
  group_by(predicted.id) %>%
  summarise(
    p.value = wilcox.test(
      MeanSmad ~ Condition
    )$p.value
  )

results$FDR <- p.adjust(results$p.value, method = "BH")

