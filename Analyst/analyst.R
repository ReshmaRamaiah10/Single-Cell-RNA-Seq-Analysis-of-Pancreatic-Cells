# Reshma Ramaiah
# This script finds all the marker genes from the clusters using sample data
# and assigns them to a cell type with known gene

library(Seurat)
library(tidyverse)

# read the sample rda file and find the markers
cells <- readRDS("GSM2230760_seurat.rda")
cells.markers <- FindAllMarkers(cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#save the markers as a rds file
saveRDS(cells.markers, file = "markers.rds")

#find all markers adn create a csv file
gen_marker_statistics <- function(x){
  cells.markers[cells.markers$cluster == x, ]
}
all_markers <- map_dfr(0:12, gen_marker_statistics)
write.csv(all_markers,'/projectnb/bf528/users/tinman_2022/project_4/Analyst/all_marker_genes.csv')

#Violin plots for different cell types to compare the clusters
VlnPlot(cells, features = c("GCG")) #Alpha
VlnPlot(cells, features = c("SST")) # Delta 
VlnPlot(cells, features = c("INS")) # Beta
VlnPlot(cells, features = c("PPY")) # Gamma 
VlnPlot(cells, features = c("KRT19")) # Ductal
VlnPlot(cells, features = c("CPA1")) # Acinar
VlnPlot(cells, features = c("SDS")) # Macrophage
VlnPlot(cells, features = c('GHRL','RGS5','PDGFRA','VWF','TRAC')) #none consistent clusters

#heatmap for top10 marker genes in each cluster and it heatmap
top10 <- cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(cells, features = top10$gene) + NoLegend()

#log normalize the counts
cells <- NormalizeData(cells, normalization.method = "LogNormalize", scale.factor = 10000)

#heatmap for top2 marker genes in each cluster
top2_normalized <- cells.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
DoHeatmap(cells, features = top2_normalized$gene, slot = "counts") + scale_fill_gradientn(colors = c("white","purple"))

#feature plot for each cell type's marker gene
FeaturePlot(cells, features = c("GCG","INS",'SST','PPY','GHRL','KRT19','CPA1','RGS5','PDGFRA','VWF','SDS','TRAC'), label = TRUE)

#signing cell type to each cluster
cell_types <- c("Delta/Gamma", "Beta", "Alpha", "Delta/Acinar", "Alpha", "Ductal","Beta", "7", "Alpha","9","Ductal","Acinar","Macrophage")
names(cell_types) <- levels(cells)
cells <- RenameIdents(cells, cell_types)

# plotting umap
cell <- RunUMAP(cells, dims = 1:10)
DimPlot(cell, reduction = "umap", label = TRUE)
