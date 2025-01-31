# spatial RNAseq clustering
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggpubr)
library(clusterProfiler)
library(clustree)
library(future)
library(openxlsx)	
library(RColorBrewer)	


workPath <- "/data/"
analysisFolder <- paste0(workPath, "scripts/analysis")
name1 <- "Anterior" 
name2 <- "Posterior" 
seu.norm1 <- readRDS(paste0(analysisFolder, "/normalized_", name1, ".rds"))
seu.norm2 <- readRDS(paste0(analysisFolder, "/normalized_", name2, ".rds"))

seu.norm1@meta.data$orig.ident <- name1
seu.norm2@meta.data$orig.ident <- name2

#parallelization
plan("multisession", workers = availableCores())
options(future.globals.maxSize = 8000 * 1024^2)


#Integration analysis
# UMAP before integration
allData <- merge(seu.norm1, seu.norm2, merge.data = TRUE)
varFeatures <- VariableFeatures(seu.norm1)
varFeatures <- c(varFeatures, VariableFeatures(seu.norm2))

#remove duplicates
varFeatures <- varFeatures[!duplicated(varFeatures)]
VariableFeatures(allData) <- varFeatures

allData <- RunPCA(allData, npcs = 50, verbose = FALSE)
allData <- RunUMAP(allData, reduction = "pca", dims = 1:50)
DimPlot(allData, reduction = "pca", group.by = "orig.ident") + ggtitle("Before integration")
DimPlot(allData, reduction = "umap", group.by = "orig.ident") + ggtitle("Before integration") 

integrated <- IntegrateLayers(object = allData, method = RPCAIntegration, normalization.method = "SCT", new.reduction = "integrated.rpca", verbose = TRUE)
integrated <- RunUMAP(integrated, dims = 1:50, reduction = "integrated.rpca")
DimPlot(integrated, reduction = "umap", group.by = "orig.ident") + ggtitle("After integration")



#Dimensionality reduction
use.pcs  <- 50
ElbowPlot(integrated, ndims=use.pcs)



#Identifying clusters
#Selecting which resolution to use:
resolution_vector <- seq(0.2,2,0.2)
integrated <- FindNeighbors(integrated, reduction="pca", dims=1:use.pcs)
integrated <- FindClusters(object=integrated, resolution=resolution_vector, verbose=FALSE)

resTable <- as.data.frame(sapply(grep("res",colnames(integrated@meta.data),value = TRUE), function(x) length(unique(integrated@meta.data[,x]))))
colnames(resTable) <- "number of clusters"
resTable

for(i in seq(resolution_vector)){
  print(DimPlot(integrated, reduction = "umap", label=TRUE, group.by=paste0("SCT_snn_res.", resolution_vector[i])) + labs(title=paste0("res_", resolution_vector[i])))
  print(SpatialDimPlot(integrated, label = TRUE, combine=FALSE, label.size = 3, group.by=paste0("SCT_snn_res.", resolution_vector[i])))
}

#Plotting clustering trees
clustree(integrated, prefis="res.")



resolution <- "0.4"
Idents(integrated) <- paste0("SCT_snn_res.", resolution)
integrated$finalCluster <- Idents(integrated)

DimPlot(integrated, reduction = "umap", label=TRUE, group.by = "ident") 
SpatialDimPlot(integrated, label = TRUE, combine=FALSE, label.size = 3)


SpatialDimPlot(integrated, combine=TRUE, images=name1, cells.highlight = CellsByIdentities(object = integrated), cols.highlight = c("red", "NA"), facet.highlight = TRUE, ncol = 4)
SpatialDimPlot(integrated, combine=TRUE, images=name2, cells.highlight = CellsByIdentities(object = integrated), cols.highlight = c("red", "NA"), facet.highlight = TRUE, ncol = 4)

#Total number of cells per cluster per slice
integrated@meta.data %>% dplyr::select(orig.ident, finalCluster) %>% table() 


#Identify marker genes for each cluster

## gvg: clustering done, set default assay back to SCT for DGE analyses
DefaultAssay(integrated) <- "SCT"

#prepare object to run differential expression on SCT assay with multiple models
integrated <- PrepSCTFindMarkers(integrated)

saveRDS(integrated, paste0(analysisFolder, "/clustered_res", resolution, ".rds"))

#find marker genes of all clusters (only overexpressed genes -> only.pos=TRUE)
markers_all <- FindAllMarkers(object=integrated, only.pos=TRUE, min.pct = 0.25)
markers_all <- markers_all[markers_all$p_val_adj < 0.05,]


#only genes unique between clusters
markers_all_single <- markers_all[markers_all$gene %in% names(table(markers_all$gene))[table(markers_all$gene) == 1],]

#nb of unique marker genes per cluster
summary <- as.data.frame(table(markers_all_single$cluster))
colnames(summary) <- c("cluster", "Nb unique marker genes")
summary 

write.xlsx(markers_all_single %>% dplyr::select(cluster, everything()), paste0(analysisFolder, "/uniqueUpregulatedGenes_perCluster_res", resolution, ".xlsx"), sheetName="Sheet1", colNames=TRUE, rowNames=TRUE, append=FALSE)


#heatmap of genes by cluster for the top 8 marker genes per cluster
top <- markers_all_single %>% group_by(cluster) %>% top_n(8, avg_log2FC)
DoHeatmap(object = integrated, features = top$gene)
# Dotplot of top 8 unique marker genes between clusters:	
DotPlot(integrated, features = top$gene, dot.scale = 8) + coord_flip() +  RotatedAxis()	


# Identify all upregulated genes for each cluster
for(cluster in sort(unique(markers_all$cluster))){
  
  markers_all_cluster <- markers_all[markers_all$cluster == cluster, ]
  rownames(markers_all_cluster) <- markers_all_cluster$gene
  
  print(paste0("### Cluster: ",cluster))
  print(paste0("Markers genes for cluster ",cluster, ":"))
  write.xlsx(markers_all_cluster, paste0(analysisFolder, "/upregulatedGenes_res", resolution, "_cluster", cluster, ".xlsx"), sheetName="Sheet1", colNames=TRUE, rowNames=TRUE, append=FALSE)
  
  print(FeaturePlot(integrated, head(markers_all_cluster$gene, n=4), cols = c("lightgrey", "blue"), ncol = 2))
  print(SpatialFeaturePlot(integrated, features = head(markers_all_cluster$gene, n=4), alpha = c(0.1, 1)))
}
  


  
#Identify spatial variable features
seu.norm1 <- FindSpatiallyVariableFeatures(seu.norm1, assay = "SCT", features = rownames(GetAssayData(seu.norm1, slot = "scale.data")), selection.method = "markvariogram")
  
spatialFeatures <- SVFInfo(seu.norm1, selection.method="markvariogram", status=TRUE)
spatialFeatures <- na.omit(spatialFeatures)
spatialFeatures <- spatialFeatures %>% filter(variable == TRUE) %>% arrange(rank)
spatialFeatures <- spatialFeatures[,c("r.metric.5"), drop=FALSE]
 
write.xlsx(spatialFeatures, paste0(analysisFolder, "/spatialVariableFeatures_res", resolution, "_", names, ".xlsx"), sheetName="Sheet1", colNames=TRUE, rowNames=TRUE, append=FALSE) 
top.features <- head(SpatiallyVariableFeatures(seu.norm1, selection.method = "markvariogram"), 9)
SpatialFeaturePlot(seu.norm1, features=top.features, ncol=3, alpha = c(0.1, 1))
