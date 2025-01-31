
library(Seurat)    
library(umap)       
library(ggpubr)    

#Loading data
workPath <- "/data/"
countFolder <- paste0(workPath, "raw_data/")
name1 <- "Anterior"
name2 <- "Posterior"
analysisFolder <- paste0(workPath, "/scripts/analysis")
dir.create(analysisFolder)

seu1 <- Load10X_Spatial(paste0(countFolder, name1), filename = "filtered_feature_bc_matrix.h5", slice = name1)
seu2 <- Load10X_Spatial(paste0(countFolder, name2), filename = "filtered_feature_bc_matrix.h5", slice = name2)

###################################################################################################################################
#Quality control - EDA
FeatureScatter(seu1, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") + labs(title = name1) + NoLegend() 
VlnPlot(seu1, features = "nFeature_Spatial")  + labs(title = name1, y = "nFeature_Spatial", x = "") + NoLegend() 

FeatureScatter(seu2, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") + labs(title = name2) + NoLegend() 
VlnPlot(seu2, features = "nFeature_Spatial")  + labs(title = name2, y = "nFeature_Spatial", x = "") + NoLegend() 


#Mitochondrial reads
seu1 <- PercentageFeatureSet(seu1, pattern = "^MT-|^Mt-|^mt-", col.name = "percent.mt")
VlnPlot(seu1, features = "percent.mt") + labs(title = name1, y = "percent of mitochondrial reads", x = "") + NoLegend() 
SpatialFeaturePlot(seu1, features = "percent.mt") + labs(title = name1)

seu2 <- PercentageFeatureSet(seu2, pattern = "^MT-|^Mt-|^mt-", col.name = "percent.mt")
VlnPlot(seu2, features = "percent.mt") + labs(title = name2, y = "percent of mitochondrial reads", x = "") + NoLegend() 
SpatialFeaturePlot(seu2, features = "percent.mt") + labs(title = name2)



# Percent of ribosomal protein reads -------
seu1 <- PercentageFeatureSet(seu1, pattern = "^RP[SL]|^MRP[SL]|^Rp[sl]|^Mrp[sl]|^rp[sl]|^mrp[sl]", col.name = "percent.rp")
VlnPlot(seu1, features = "percent.rp") + labs(title = name1, y = "percent of ribosomal reads", x = "") + NoLegend() 
SpatialFeaturePlot(seu1, features = "percent.rp")  + labs(title = name1)

seu2 <- PercentageFeatureSet(seu2, pattern = "^RP[SL]|^MRP[SL]|^Rp[sl]|^Mrp[sl]|^rp[sl]|^mrp[sl]", col.name = "percent.rp")
VlnPlot(seu2, features = "percent.rp") + labs(title = name2, y = "percent of ribosomal reads", x = "") + NoLegend() 
SpatialFeaturePlot(seu2, features = "percent.rp")  + labs(title = name2)



#Quality control - filtering
# Identifying low-quality cells -------
seu.filt1 <- seu1[, seu1$nFeature_Spatial > 500 & seu1$percent.mt < 32]
seu.filt2 <- seu2[, seu2$nFeature_Spatial > 500 & seu2$percent.mt < 32]

# Gene-level QC
C <- seu.filt1[["Spatial"]]$counts
C@x <- C@x/rep.int(colSums(C), diff(C@p))
most_expressed <- order(Matrix::rowSums(C), decreasing = TRUE)[30:1]
boxplot(as.matrix(t(as.matrix(C[most_expressed, ]))), cex.axis=0.5, cex.lab=0.8 , cex = 0.1, las = 1, xlab = "% total count per cell", col = (scales::hue_pal())(30)[30:1], horizontal = TRUE)

C <- seu.filt2[["Spatial"]]$counts
C@x <- C@x/rep.int(colSums(C), diff(C@p))
most_expressed <- order(Matrix::rowSums(C), decreasing = TRUE)[30:1]
boxplot(as.matrix(t(as.matrix(C[most_expressed, ]))), cex.axis=0.5, cex.lab=0.8 , cex = 0.1, las = 1, xlab = "% total count per cell", col = (scales::hue_pal())(30)[30:1], horizontal = TRUE)


remove_genes <- c("mt-Co3", "mt-Co1", "mt-Atp6", "mt-Cytb", "mt-Co2", "mt-Nd4", "mt-Nd2", "mt-Nd1")

if(!is.null(remove_genes)){
  remove_feature <- rownames(seu.filt1) %in% remove_genes
  seu.filt1 <- seu.filt1[!remove_feature, ]
  
  C <- seu.filt1[["Spatial"]]$counts
  C@x <- C@x/rep.int(colSums(C), diff(C@p))
  most_expressed <- order(Matrix::rowSums(C), decreasing = T)[30:1]
  print(boxplot(as.matrix(t(as.matrix(C[most_expressed, ]))), cex.axis=0.5, cex.lab=0.8 , cex = 0.1, las = 1, xlab = "% total count per cell", col = (scales::hue_pal())(30)[30:1], horizontal = TRUE))
  
  
  remove_feature <- rownames(seu.filt2) %in% remove_genes
  seu.filt2 <- seu.filt2[!remove_feature, ]
  
  C <- seu.filt2[["Spatial"]]$counts
  C@x <- C@x/rep.int(colSums(C), diff(C@p))
  most_expressed <- order(Matrix::rowSums(C), decreasing = T)[30:1]
  print(boxplot(as.matrix(t(as.matrix(C[most_expressed, ]))), cex.axis=0.5, cex.lab=0.8 , cex = 0.1, las = 1, xlab = "% total count per cell", col = (scales::hue_pal())(30)[30:1], horizontal = TRUE))
  
}



# Statistics after filtering - Number of cells 
stats <- c(dim(seu1)[2], dim(seu.filt1)[2])
stats <- rbind(stats, c(dim(seu2)[2], dim(seu.filt2)[2]))
colnames(stats) <- c("Before filtering", "After filtering")
rownames(stats) <- c(name1, name2)
stats 

# Library size versus detected genes after filtering
FeatureScatter(seu.filt1, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") + labs(title = name1) + NoLegend() 
VlnPlot(seu.filt1, features = "nFeature_Spatial")  + labs(title = name1, y = "nFeature_Spatial", x = "") + NoLegend() 

FeatureScatter(seu.filt2, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") + labs(title = name2) + NoLegend() 
VlnPlot(seu.filt2, features = "nFeature_Spatial")  + labs(title = name2, y = "nFeature_Spatial", x = "") + NoLegend() 


# Percent of mitochondrial reads after filtering
VlnPlot(seu.filt1, features = "percent.mt") + labs(title = name1, y = "percent of mitochondrial reads", x = "") + NoLegend() 
SpatialFeaturePlot(seu.filt1, features = "percent.mt") + labs(title = name1)

VlnPlot(seu.filt2, features = "percent.mt") + labs(title = name2, y = "percent of mitochondrial reads", x = "") + NoLegend() 
SpatialFeaturePlot(seu.filt2, features = "percent.mt") + labs(title = name2)


# Percent of ribosomal protein reads after filtering
VlnPlot(seu.filt1, features = "percent.rp") + labs(title = name1, y = "percent of ribosomal reads", x = "") + NoLegend() 
SpatialFeaturePlot(seu.filt1, features = "percent.rp")  + labs(title = name1)

VlnPlot(seu.filt2, features = "percent.rp") + labs(title = name2, y = "percent of ribosomal reads", x = "") + NoLegend() 
SpatialFeaturePlot(seu.filt2, features = "percent.rp")  + labs(title = name2)



###################################################################################################################################
# Normalization
# before Normalization:
VlnPlot(seu.filt1, features = "nCount_Spatial", pt.size = 0.1) + labs(title=name1, y="nCount_Spatial", x="") + NoLegend() 
SpatialFeaturePlot(seu.filt1, features = "nCount_Spatial") + labs(title = name1) + theme(legend.position = "right")

VlnPlot(seu.filt2, features = "nCount_Spatial", pt.size = 0.1) + labs(title=name2, y="nCount_Spatial", x="") + NoLegend() 
SpatialFeaturePlot(seu.filt2, features = "nCount_Spatial") + labs(title = name2) + theme(legend.position = "right")

# after sctransform normalization:
# Apply sctransform normalization:
seu.norm1 <- SCTransform(seu.filt1, assay="Spatial", vars.to.regress=c("percent.mt", "percent.rp"), verbose=FALSE)
seu.norm2 <- SCTransform(seu.filt2, assay="Spatial", vars.to.regress=c("percent.mt", "percent.rp"), verbose=FALSE)
saveRDS(seu.norm1, paste0(analysisFolder, "/normalized_", name1, ".rds"))
saveRDS(seu.norm2, paste0(analysisFolder, "/normalized_", name2, ".rds"))

VlnPlot(seu.norm1, features="nCount_SCT", pt.size = 0.1) + labs(title=name1, y="nCount_Spatial", x="") + NoLegend()
SpatialFeaturePlot(seu.norm1, features="nCount_SCT") + labs(title = name1) + theme(legend.position = "right")
VlnPlot(seu.norm2, features="nCount_SCT", pt.size = 0.1) + labs(title=name2, y="nCount_Spatial", x="") + NoLegend()
SpatialFeaturePlot(seu.norm2, features="nCount_SCT") + labs(title = name2) + theme(legend.position = "right")

