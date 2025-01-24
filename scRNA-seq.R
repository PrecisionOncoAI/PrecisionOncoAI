# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Set working directory (modify accordingly)
setwd("path/to/your/data")

# Load single-cell RNA-seq data (e.g., from 10X Genomics)
# If using 10X, provide the directory containing 'barcodes.tsv.gz', 'features.tsv.gz', and 'matrix.mtx.gz'
sc_data <- Read10X(data.dir = "path/to/10X_data")

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = sc_data, min.cells = 3, min.features = 200, project = "scRNA_seq_analysis")

# View metadata
print(seurat_obj)

# Quality control
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-") # Mitochondrial gene percentage

# Plot QC metrics
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filtering criteria (adjust thresholds as needed)
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Feature selection (highly variable genes)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scaling the data
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

# verify by elbow plot
ElbowPlot(four_10006BL)

# Perform PCA for dimensionality reduction
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Visualize PCA results
VizDimLoadings(seurat_obj, dims = 1:10, reduction = "pca")
DimPlot(seurat_obj, reduction = "pca")

#Check by Elbow plot to see how many PC needed
ElbowPlot(seurat_obj)

# Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# UMAP/t-SNE visualization
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
seurat_obj <- RunTSNE(seurat_obj, dims = 1:10)

DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(seurat_obj, reduction = "tsne", label = TRUE, pt.size = 0.5)

# Find differentially expressed genes (DEGs) in each cluster
cluster_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# View top markers
head(cluster_markers)

# Save markers
write.csv(cluster_markers, file = "cluster_markers.csv")

# Feature plot for specific marker genes
FeaturePlot(seurat_obj, features = c("CD3D", "MS4A1", "CD14")) # Modify with genes of interest

# Save Seurat object
saveRDS(seurat_obj, file = "seurat_obj.rds")

# singleR annotation
#BiocManager::install("SingleR")
library(SingleR)
#BiocManager::install("celldex")
library(celldex)

hpca.se <- celldex::HumanPrimaryCellAtlasData()

#hpca.se <- HumanPrimaryCellAtlasData()

##Need to run before find marker
singler <- SingleR(test = GetAssayData(seurat_obj, assay = "RNA",slot = "data"),clusters=Idents(seurat_obj),ref=hpca.se,assay.type.test=1,labels = hpca.se$label.main)

write.table(singler,"seurat_obj_singleR.txt")
