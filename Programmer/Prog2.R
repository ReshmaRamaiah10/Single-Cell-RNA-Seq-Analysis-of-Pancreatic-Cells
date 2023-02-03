# Author: Allison Choy
# BF528: Project 4
# Role: Programmer

setwd('/projectnb/bf528/users/tinman_2022/project_4/Programmer')

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

# Remove '#' to install necessary packages
#BiocManager::install("Seurat")
#BiocManager::install("tximport")
#remotes::install_github("satijalab/seurat-wrappers") # it may take a while to install

library(tidyverse)
library(Seurat)
library(tximport)
library(SeuratWrappers)
library(biomaRt)
library(ggpubr)

path <- file.path("/projectnb/bf528/project_4_scrnaseq/GSM2230760__salmon_quant/alevin/quants_mat.gz")
UMI_data <- tximport(path, type="alevin")
dim(UMI_data$counts) # [1]  60233 108832

# Looks at data to determine where the cutoff will be
counts <- UMI_data$counts
counts_per_cell <- Matrix::colSums(counts)
counts_per_gene <- Matrix::rowSums(counts)
genes_per_cell <- Matrix::colSums(counts>0) # count a gene only if it has non-zero reads mapped.
cells_per_gene <- Matrix::rowSums(counts>0) # only count cells where the gene is expressed

MIN_GENES_PER_CELL <- 50
MINI_GENES_PER_CELL <- 200
MID_GENES_PER_CELL <- 500  ## user-defined setting
MAX_GENES_PER_CELL <- 4500  ## user-defined setting

plot(sort(genes_per_cell), xlab='Cell', log='y', ylab='Sorted Genes per Cell', 
     main='Data Complexity with Ordered Genes per Cell')
abline(h=MIN_GENES_PER_CELL, col='blue')
abline(h=MINI_GENES_PER_CELL, col='magenta')
abline(h=MID_GENES_PER_CELL, col='orange')  # lower threshold
abline(h=MAX_GENES_PER_CELL, col='gold') # upper threshold

# edits out the version number in geneIDs
UMI_edit <- sub("[.][0-9]*$", "", rownames(UMI_data$counts))
UMI_data$counts@Dimnames[[1]] <- UMI_edit

# Gets the gene symbols from Ensembl
human <-useEnsembl('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl') #, mirror = 'useast')
mouse <- useEnsembl("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
human_en <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                  filters='ensembl_gene_id', 
                  values=UMI_edit,
                  mart = human)

# no mouse genes found
# using:
genes2 <- as.matrix(UMI_edit)
genes2[str_detect(genes2,"^ENSMUS"),] # character(0) for mouse IDs
head(genes2[str_detect(genes2,"^ENSG"),]) # returns the whole matrix otherwise

mouse_en <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                  filters='ensembl_gene_id', 
                  values=UMI_edit,
                  mart = mouse) # 0 observations

write.csv(human_en, "./human_ensembl2.csv", row.names = FALSE)
write.csv(mouse_en, "./mouse_ensembl.csv", row.names = FALSE)

# filter geneids that do not have a symbol
filtered_genes <- tibble(human_en) %>% filter(!hgnc_symbol == '')
dim(filtered_genes) #[1] 39472     2

# filtering counts
Ucounts <- UMI_data$counts
c2 <- Ucounts[rownames(Ucounts) %in% filtered_genes$ensembl_gene_id, ]
dim(c2) #[1]  39470 108832

c3 <- as_tibble(c2@Dimnames[[1]]) %>% left_join (filtered_genes, by=c('value' = 'ensembl_gene_id')) %>% 
  filter(!duplicated(hgnc_symbol)) %>% filter(!duplicated(value))

c4 <- c3 %>% pull(hgnc_symbol) 
dim(c3)

# went from # A tibble: 39,470 × 1 (start with filtered c2), # A tibble: 39,472 × 2 (left join),
# A tibble: 39,458 × 2 (filtered for duplication), then finally # A tibble: 39,456 × 2

c2 <- Ucounts[rownames(Ucounts) %in% c3$value, ]
dim (c2) #[1]  39456 108832, slightly less

c2@Dimnames[[1]] <- c4
pan_cells <- CreateSeuratObject(counts = c2, project = "pancreas", min.cells = 3, min.features = 50)
pan_cells # 21191 features across 43759 samples within 1 assay 
# with 50 : 21347 features across 71126 samples within 1 assay 
#Active assay: RNA (21347 features, 0 variable features)

# with cutoff of 100: 21185 features across 43753 samples within 1 assay 
#Active assay: RNA (21185 features, 0 variable features)
# with cutoff at 30: 21380 features across 82618 samples within 1 assay 
#Active assay: RNA (21380 features, 0 variable features)
# with 150: 20968 features across 23583 samples within 1 assay 
#Active assay: RNA (20968 features, 0 variable features)
# with 200: 20742 features across 13301 samples within 1 assay 
# Active assay: RNA (20742 features, 0 variable features)

pan_cells[["percent.mt"]] <- PercentageFeatureSet(pan_cells, pattern = "^MT-")

VlnPlot(pan_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + NoLegend()

plot1a <- FeatureScatter(pan_cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2a <- FeatureScatter(pan_cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1a + NoLegend() + plot2a + NoLegend()

pan_cells2 <- subset(pan_cells, subset = nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 30)
pan_cells2
# mt < 30
#21347 features across 2608 samples within 1 assay 
#Active assay: RNA (21347 features, 0 variable features)

# mt < 10%
#21185 features across 1248 samples within 1 assay 
#Active assay: RNA (21185 features, 0 variable features)

# with < 5% MT
#21347 features across 199 samples within 1 assay 
#Active assay: RNA (21347 features, 0 variable features)

# Normalization of data
#pan_cells2 <- NormalizeData(pan_cells2, normalization.method = "LogNormalize", scale.factor = 10000) # or
pan_cells2 <- NormalizeData(pan_cells2)

# Identification of highly variable features
pan_cells2 <- FindVariableFeatures(pan_cells2, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pan_cells2), 10)

# plot variable features with and without labels
plot1b <- VariableFeaturePlot(pan_cells2)
plot2b <- LabelPoints(plot = plot1b, points = top10, repel = TRUE)
plot1b + plot2b + NoLegend()

# Scaling the data before PCA
all_genes <- rownames(pan_cells2)
pan_cells2 <- ScaleData(pan_cells2, features = all_genes) # vars.to.regress = nFeature_RNA, percent.mt

#Run PCA
pan_cells2 <- RunPCA(pan_cells2, features = VariableFeatures(object = pan_cells2))

# Visualization of PCA
a <- VizDimLoadings(pan_cells2, dims = 1:2, reduction = "pca")
b <- DimPlot(pan_cells2, reduction = "pca")
c <- DimHeatmap(pan_cells2, dims = 1, cells = 500, balanced = TRUE) # doesn't save

# Dimensionality of dataset
d<- ElbowPlot(pan_cells2)

# JackStrawPlot
pan_cells2 <- JackStraw(pan_cells2, num.replicate = 100)
pan_cells2 <- ScoreJackStraw(pan_cells2, dims = 1:20)
e <- JackStrawPlot(pan_cells2, dims = 1:20)

# Clustering of Cells
pan_cells2 <- FindNeighbors(pan_cells2, dims = 1:10)
pan_cells2 <- FindClusters(pan_cells2, resolution = 0.5)

# Looking at the first 5 cells
head(Idents(pan_cells2), 5)

#UMAP/tSNE
pan_cells2 <- RunUMAP(pan_cells2, dims = 1:10)
pan_UMAP <- DimPlot(pan_cells2, reduction = "umap", label = TRUE)
pan_UMAP

saveRDS(pan_cells2, file = "/projectnb/bf528/users/tinman_2022/project_4/Programmer/tinman_proj4.rds")

# creates a tibble for props of pie chart
pan_df <- as.data.frame(table(Idents(pan_cells2))) %>% as_tibble() %>% 
  summarize(Prop = round((Freq/sum(Freq))*100, 2)) %>% mutate(Cluster = 0:8) %>% 
  mutate(yloc = cumsum(Prop)- 0.5*Prop ) %>% mutate(Label = paste0(Prop, '%'))

#ggplot attempt for pie chart
pan_pie <- ggplot(pan_df, aes(x="", y=Prop, fill=Cluster))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = 'y', start = 0) +
  theme_void() + 
  geom_text(aes(y = yloc, label = Label), color = "white", size=3) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "gold", "#009E73", "#F0E442", "#0072B2", "magenta", "#D55E00", "#CC79A7"),
                    name = "Cluster")

# actually used
plot_col <- c("#E69F00", "gold", "#F0E442", "#009E73", "#56B4E9", "#0072B2", "magenta", "#CC79A7", "#D55E00")
pie(pan_df$Prop, 
    labels = pan_df$Label, 
    border = 'white', 
    angle=c(20,90,30,10,0,100,270,45,300), 
    col = plot_col,
    main = "Relative Proportions of Cell Numbers\nPer Cluster")
legend(1.2, 1, legend = pan_df$Cluster, fill = plot_col, bty="n", title="Cluster")
