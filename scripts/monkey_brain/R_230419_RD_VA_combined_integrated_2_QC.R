## R 4.1.2
## Seurat 4.1.1
## dplyr 1.0.10
## patchwork 1.1.2
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
projectName <- "RD_VA_01_12"

# Seurat parameters
variable.features <- 2000
nPCs <- 50


# Load integrated umap data
data.combined <- readRDS("data_Integrated_12_samples_umap.RDS")

DefaultAssay(data.combined) <- "integrated"
# check integrated umap
DimPlot(data.combined, reduction = "umap", group.by = "group", raster=FALSE,pt.size =0.001)


DefaultAssay(data.combined) <- "RNA"

# tmp <- data.combined@assays$RNA@counts@Dimnames[[1]]
# tmp[grepl(tmp,'MT')]
# 
# gene_name <- as.vector(data.combined@assays$RNA@counts@Dimnames[[1]])
# gene_name[grepl('MT', gene_name, ignore.case = TRUE)]

# QC 
### Quality Control ###
### QC step 1: nFeature_RNA, nCount_RNA, Mitochondrial RNA, and Ribosome RNA 


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# percentage of Mitochondrial RNA

mt_genes <- read.csv('MTgene.txt', header=FALSE)

library(stringr)
mt_genes$V2 = str_replace_all(mt_genes$V2, " ", "") 

modify_mt_gene <- mt_genes$V2[mt_genes$V2 %in% data.combined@assays$RNA@counts@Dimnames[[1]]]
data.combined[["percent.mt"]] <- PercentageFeatureSet(data.combined, features = modify_mt_gene)

VlnPlot(data.combined, features = c("percent.mt"), ncol = 1, pt.size =0.01)

#data.combined[["percent.mt"]] <- PercentageFeatureSet(data.combined, pattern = "MT")

# percentage of Ribsome RNA
#data.combined[["percent.rp"]] <- PercentageFeatureSet(data.combined, pattern = "RP")

# VlnPlot(data.combined, features = c("percent.rp"), ncol = 1)
# VlnPlot(data.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)

# plot and save VlnPlots of nFeature_RNA, nCount_RNA, Mitochondrial RNA, and Ribosome RNA 
pdf(paste("pdf_QC_VlnPlot_nF_nC_ptMT_ptRP.pdf"), width=10, height=8, useDingbats = FALSE)
  VlnPlot(data.combined, features = c("nFeature_RNA"), ncol = 1, pt.size =0.01)
  VlnPlot(data.combined, features = c("nCount_RNA"), ncol = 1, pt.size =0.01)
  VlnPlot(data.combined, features = c("percent.mt"), ncol = 1, pt.size =0.01)
  VlnPlot(data.combined, features = c("percent.rp"), ncol = 1, pt.size =0.01)
dev.off()

# check FeatureScatter Plots of nFeature_RNA, Mitochondrial RNA, and Ribosome RNA vs nCount_RNA
ptMT_nC <- FeatureScatter(data.combined, feature1 = "nCount_RNA", feature2 = "percent.mt",raster=FALSE,pt.size = 0.1)
ptRP_nC <- FeatureScatter(data.combined, feature1 = "nCount_RNA", feature2 = "percent.mt",raster=FALSE,pt.size = 0.1)
nF_nC <- FeatureScatter(data.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",raster=FALSE,pt.size = 0.1)

ptMT_nC
ptRP_nC
nF_nC

# plot and save FeatureScatter Plots of Mitochondrial RNA vs nCount_RNA
pdf(paste("pdf_QC_ptMT_vs_nCount_RNA.pdf"), width=10, height=8, useDingbats = FALSE)
  ptMT_nC
  for (sample in sample_Names_12) {
  print(FeatureScatter(subset(data.combined,subset = (group %in% c(sample_Names_12[1]))), 
                       feature1 = "nCount_RNA", feature2 = "percent.mt",raster=FALSE,pt.size = 0.1))
  }
dev.off()

# plot and save FeatureScatter Plots of Ribosome RNA vs nCount_RNA
pdf(paste("pdf_QC_ptRP_vs_nCount_RNA.pdf"), width=10, height=8, useDingbats = FALSE)
  ptRP_nC
  for (sample in sample_Names_12) {
    print(FeatureScatter(subset(data.combined,subset = (group %in% c(sample))), 
                         feature1 = "nCount_RNA", feature2 = "percent.rp",raster=FALSE,pt.size = 0.1))
  }
dev.off()

# plot and save FeatureScatter Plots of nFeature_RNA vs nCount_RNA
pdf(paste("pdf_QC_nFeature_vs_nCount_RNA.pdf"), width=10, height=8, useDingbats = FALSE)
nF_nC
for (sample in sample_Names_12) {
  print(FeatureScatter(subset(data.combined,subset = (group %in% c(sample))), 
                       feature1 = "nCount_RNA", feature2 = "nFeature_RNA",raster=FALSE,pt.size = 0.1))
}
dev.off()

# Remove unqualified cells
data.combined <- subset(data.combined, subset = nFeature_RNA < 6000 & nCount_RNA < 30000 & percent.mt < 5 & percent.rp < 5)

# check and save VlnPlots of nFeature_RNA, nCount_RNA, Mitochondrial RNA, and Ribosome RNA after QC
pdf(paste("pdf_after_QC_VlnPlot_nF_nC_ptMT_ptRP.pdf"), width=10, height=8, useDingbats = FALSE)
  VlnPlot(data.combined, features = c("nFeature_RNA"), ncol = 1, pt.size =0.01)
  VlnPlot(data.combined, features = c("nCount_RNA"), ncol = 1, pt.size =0.01)
  VlnPlot(data.combined, features = c("percent.mt"), ncol = 1, pt.size =0.01)
  VlnPlot(data.combined, features = c("percent.rp"), ncol = 1, pt.size =0.01)
dev.off()

# QC 
### Quality Control ###
### QC step 2: Predict and Remove Doublets

###############################################
## DoubletFinder
## Predict doublets
## remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
################################################
suppressMessages(require(DoubletFinder))
## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------

DefaultAssay(data.combined) <- "RNA"

data.combined <- NormalizeData(data.combined)
data.combined <- FindVariableFeatures(data.combined, selection.method = "vst", nfeatures = 2000)
data.combined <- ScaleData(data.combined)
data.combined <- RunPCA(data.combined)
data.combined <- RunUMAP(data.combined, dims = 1:10)

DefaultAssay(data.combined) <- "integrated"
# check umap
DimPlot(data.combined, reduction = "umap", group.by = "group", raster=FALSE)

str(data.combined)

pdf(paste("pdf_QC_umap_1.pdf"), width=10, height=8, useDingbats = FALSE)
DimPlot(data.combined, reduction = "umap", group.by = "group", raster=FALSE)
for (sample in sample_Names_12) {
  print(DimPlot(subset(data.combined,subset = (group %in% c(sample))),reduction = "umap",group.by = "group", raster=FALSE))
}
dev.off()

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res <- paramSweep_v3(data.combined, PCs = 1:10, sct = FALSE)
sweep.res <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn_1028 <- find.pK(sweep.res)
## pK=0.24

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- data.combined@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(data.combined@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
data.combined <- doubletFinder_v3(data.combined, PCs = 1:10, pN = 0.25, pK = 0.24, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
data.combined <- doubletFinder_v3(data.combined, PCs = 1:10, pN = 0.25, pK = 0.24, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_4835", sct = FALSE)

# name of the DF prediction can change, so extract the correct column name.
DF.name = colnames(data.combined@meta.data)[grepl("DF.classification", colnames(data.combined@meta.data))]
cowplot::plot_grid(ncol = 2, DimPlot(data.combined, group.by = "Mouse") + NoAxes(),
                   DimPlot(data.combined, group.by = DF.name) + NoAxes())



