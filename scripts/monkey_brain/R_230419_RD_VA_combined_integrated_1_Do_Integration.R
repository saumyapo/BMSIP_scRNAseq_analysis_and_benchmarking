library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)

packageVersion('dplyr')
packageVersion('Seurat')
packageVersion('patchwork')
packageVersion('ggplot2')
packageVersion('cowplot')

# packageVersion('dplyr')
# [1] ‘1.0.7’
# packageVersion('Seurat')
# [1] ‘4.1.0’
# packageVersion('patchwork')
# [1] ‘1.1.1’
# library(ggplot2)
# packageVersion('ggplot2')
# [1] ‘3.3.5’
# packageVersion('cowplot')
# [1] ‘1.1.1’

suffix <- "Seurat"
sampleNames <- c("RD_VA_01", "RD_VA_02", "RD_VA_03", "RD_VA_04", "RD_VA_05", "RD_VA_06", 
                 "RD_VA_07", "RD_VA_08", "RD_VA_09", "RD_VA_10","RD_VA_11","RD_VA_12")
sample_Names_12 = gsub("_", "-", sampleNames)
rawDataDir_1 <- "/projectnb/czlab/A10_Rosene/Data/20230408_combinedTwo_results/"
rawDataDir_2 <- "/outs/filtered_feature_bc_matrix"
projectName <- "RD_VA_01_12"

# Seurat parameters
variable.features <- 2000
nPCs <- 50


# Load data
samples <- c()
for (sample in sampleNames) {
  name <- paste(rawDataDir_1, sample, rawDataDir_2, sep = "")
  samples <- c(samples, name)
}
names(samples) = sample_Names_12 ## No dashes should be in the sample name as that conflicts with Seurat's nomenclature
all.raw.data <- Read10X(data.dir = samples)
data <- CreateSeuratObject(counts = all.raw.data, project = projectName,  min.cells = 3, min.features = 200)#, min.cells = seuratQCMinCells, min.features = seuratQCMinGenes,)

DefaultAssay(data)

# Save raw data
saveRDS(data, file = paste("data_raw.RDS") )
# Load raw data
data <- readRDS("data_raw.RDS")
data$group <- data@meta.data$orig.ident

# split the dataset into a list of two seurat objects (stim and CTRL)
data.list = SplitObject(data, split.by = "ident")

# normalize and identify variable features for each dataset independently
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = data.list)


# Perform integration
# We then identify anchors using the FindIntegrationAnchors() function, which takes a list 
# of Seurat objects as input, and use these anchors to integrate the two datasets together with IntegrateData().
data.anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features)

# this command creates an 'integrated' data assay
data.combined <- IntegrateData(anchorset = data.anchors)

# Save raw data
saveRDS(data.combined, file = paste("data_Integrated_12_samples.RDS") )
# Load raw data
data.combined <- readRDS("data_Integrated_12_samples.RDS")

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(data.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
data.combined <- ScaleData(data.combined, verbose = FALSE)
data.combined <- RunPCA(data.combined, npcs = 50, verbose = FALSE)
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:50)

# check integrated umap
DimPlot(data.combined, reduction = "umap", group.by = "group")

# Save integrated data
saveRDS(data.combined, file = paste("data_Integrated_12_samples_umap.RDS") )
# Load integrated data
data.combined <- readRDS("data_Integrated_12_samples_umap.RDS")

# Plot and save integrated umap 
pdf(paste("Integrated_12_samples_umap.pdf"), width=10, height=8, useDingbats = FALSE)
DimPlot(data.combined, reduction = "umap", group.by = "group")
for (sample in sample_Names_12) {
  print(DimPlot(subset(data.combined,subset = (group %in% c(sample))),reduction = "umap",group.by = "group"))
}
dev.off()
