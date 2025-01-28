#Compiled script for scRNA seq

library("AnnotationDbi") 
library("org.Hs.eg.db")
library("dplyr")
library(tidyr)
library(Seurat)

# Loading and processing the data

f <- read.csv("EXP0050_PCG_beforeQC.txt", sep="\t", row.names = 1)
f <- f[-1,]
#View(f[1:10,1:10])

f$Gene <- mapIds(org.Hs.eg.db,
                   keys=row.names(f),
                   column="SYMBOL", 
                   keytype="ENSEMBL", 
                   multiVals="first")

f <- f %>% select(Gene, everything())
#View(f$Gene)

f <- distinct(f,Gene, .keep_all= TRUE)
f <- f %>% drop_na(Gene)

row.names(f)<-f$Gene
f$Gene <- NULL


#Filtering data for min cells the gene is present in and min features the cell contains
pb <- CreateSeuratObject(counts = f, 
                         min.cells = 3,
                         min.features = 200)


#filtering mitochondrial genes
pb[["percent.mt"]] <- PercentageFeatureSet(pb, pattern = "^MT-")
pb[["percent.mt"]]

# Visualize QC metrics as a violin plot
pdf(file="Plot1.pdf")
VlnPlot(pb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Normalization
pb <- NormalizeData(pb)

# Variable features
pb <- FindVariableFeatures(pb, selection.method = "vst", nfeatures = 2000)

list_of_variable_features <- VariableFeatures(pb)
list_of_variable_features <- as.data.frame(list_of_variable_features)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pb), 10)


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pb)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(pb)
pb <- ScaleData(pb, features = all.genes)

pb <- RunPCA(pb, features = VariableFeatures(object = pb))
VizDimLoadings(pb, dims = 1:2, reduction = "pca")
DimPlot(pb, reduction = "pca")

DimHeatmap(pb, dims = 1:15, cells = 500, balanced = TRUE)
pb <- JackStraw(pb, num.replicate = 100)
pb <- ScoreJackStraw(pb, dims = 1:20)
JackStrawPlot(pb, dims = 1:15)
ElbowPlot(pb)

pb <- FindNeighbors(pb, dims = 1:10)
pb <- FindClusters(pb, resolution = 0.5)
head(Idents(pb), 5)
pb <- RunUMAP(pb, dims = 1:10)
DimPlot(pb, reduction = "umap")

pb.markers <- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(pb.markers,file="Markers_info.csv",
            sep=",",row.names=FALSE,col.names = TRUE,quote = FALSE)

#VlnPlot(pb, features = c("OR1A1","MS4A1","CD79A"))
VlnPlot(pb, features = c("OR1A1","RAP1A","JCHAIN"))
VlnPlot(pb, features = c("JCHAIN"))
VlnPlot(pb, features = c("HBB"))

FeaturePlot(pb, features = "JCHAIN")
FeaturePlot(pb, features = "JCHAIN", pt.size = 2)


check <- Idents(pb)
check <- as.data.frame(check)
library(data.table)

check<-setDT(check, keep.rownames = TRUE)[]
colnames(check) <- c("Cell_Id", "Cluster")

# Metadata
meta_data <- read.csv("SraRunTable.txt", sep=",")
meta_data <- meta_data[,c(1,5,21,25,27)]

f <- merge(check, meta_data, by.x="Cell_Id", by.y="Run")

