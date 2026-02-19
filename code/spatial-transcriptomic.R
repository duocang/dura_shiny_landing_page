#===== 1.0 set up =====
#### 1.1 environment ####
library(Seurat)     
library(SeuratWrappers) 
library(patchwork)  
library(tidyverse) 
library(grid) 
library(viridis) 


#### 1.2 prepare MS sample ####
##### 1.2.1 load dataset ####
raw_data_dir <- #define directory
setwd(raw_data_dir)
MS_sample <- #file_name
sc <- LoadXenium(MS_sample, fov = "fov")
sc$orig.ident <- "dura"
sc$sample_ID <- "MS_sample"
Idents(sc) <- "dura"

result_dir <- #define directory
setwd(result_dir)
qsave(sc, file.path(result_dir, "1.2_Xenium_MS_sample_seurat.qs"))

##### 1.2.2 quality control ####
###### 1.2.2.1 Number of Counts per Bin ####
#bins = cells
#The total number of UMIs (counts) per bin gives a sense of overall transcript abundance.
p1 <- VlnPlot(sc, features="nCount_Xenium", pt.size=0, layer="counts")
p2 <- ImageFeaturePlot(sc, features="nCount_Xenium", dark.background = FALSE, border.size = NA, size = 1)
p3 <- p1 | p2
ggsave(file.path(result_dir, "1.2.2.1_QC_number_of_counts_per_cells.pdf"), p3, width = 7, heigh = 5)
summary(sc$nCount_Xenium)

#Look for bins with very few or very many counts. Bins with <30 counts or >2000 counts are suspicious and flagged for removal.
sc$filtered_by_counts = ifelse(sc$nCount_Xenium < 30 | sc$nCount_Xenium > 2000, "filtered-out", "keep")
p <- ImageDimPlot(sc, group.by="filtered_by_counts", border.size = NA, dark.background = FALSE, size = 0.5)
ggsave(file.path(result_dir, "1.2.2.1_QC_filter_by_counts.pdf"), p, width = 5, heigh = 5)

###### 1.2.2.2 Number of Genes per Bin ####
p1 <- VlnPlot(sc, features="nFeature_Xenium", pt.size=0, layer="counts") + NoLegend()
p2 <- ImageFeaturePlot(sc, features="nFeature_Xenium", dark.background = FALSE, border.size = NA, size = 1)
p3 <- p1 | p2
ggsave(file.path(result_dir, "1.2.2.2_QC_number_of_genes_per_cells.pdf"), p3, width = 7, heigh = 5)

summary(sc$nFeature_Xenium)

sc$filtered_by_genes = ifelse(sc$nFeature_Xenium < 5 | sc$nFeature_Xenium > 1000, "filtered-out", "keep")
p1 <- ImageDimPlot(sc, group.by="filtered_by_genes", border.size = NA, dark.background = FALSE, size = 0.5)
p2 <- VlnPlot(sc, features="nFeature_Xenium", pt.size=0, layer="counts", group.by = "filtered_by_genes") + NoLegend()
p3 <- p1 | p2
ggsave(file.path(result_dir, "1.2.2.2_QC_filter_by_genes.pdf"), p3, width = 7, heigh = 5)

##### 1.2.3 Filtering #####
dim(sc)
###### 1.2.3.1 features (and Counts) #####
sc_filtered <- subset(sc, nFeature_Xenium < 5 | nFeature_Xenium > 1000, invert=TRUE) %>% suppressWarnings() 
#sc_filtered <- subset(sc_filtered, nCount_Xenium < 30 | nCount_Xenium > 2000, invert=TRUE) %>% suppressWarnings()
qsave(sc_filtered, file.path(result_dir, "1.2.3.1_Xenium_NBB_Ctrl_seurat_filtered.qs"))

p1 <- ImageDimPlot(sc, border.size = NA, dark.background = FALSE, size = 0.5) + ggtitle("not filtered")
p2 <- ImageDimPlot(sc_filtered, border.size = NA, dark.background = FALSE, size = 0.5) + ggtitle("filtered")
p3 <- p1 | p2

###### 1.2.3.2 genes ####
#Genes expressed in only a few bins can add noise. First, we count how many bins each gene is detected in:
counts <- GetAssayData(sc, layer="counts")
genes_expressed <- rowSums(counts > 0)
head(genes_expressed, 5)
summary(genes_expressed)

#Keep genes expressed in >10 bins:
genes_expressed <- genes_expressed[genes_expressed > 10]
length(genes_expressed)
#133 (= no genes need to be filtered out)

ggsave(file.path(result_dir, "1.2.3.2_filtered.pdf"), p3, width = 7, heigh = 5)


##### 1.2.4 Pre-processing ####
###### 1.2.4.1 Normalization ####
sc_pre <- NormalizeData(sc_filtered, normalization.method="LogNormalize", scale.factor=10000)

###### 1.2.4.2 Identifying Highly Variable Genes (HVGs) ####
sc_pre <- FindVariableFeatures(sc_pre, nfeatures=3000)

top20 <- VariableFeatures(sc_pre) %>% head(20)
p <- VariableFeaturePlot(sc_pre, selection.method="vst")
p <- LabelPoints(plot=p, points=top20, repel=TRUE)
p
ggsave(file.path(result_dir, "1.2.4.2_VariableFeatures_top.pdf"), p, width = 7, heigh = 5)

###### 1.2.4.3 Data scaling ####
sc_pre <- ScaleData(sc_pre, features=VariableFeatures(sc_pre))

###### 1.2.4.4 Dimensionality reduction ####
sc_pre <- RunPCA(sc_pre, reduction.name="pca", features=VariableFeatures(sc_pre), npcs=50, nfeatures.print=5)

p1 <- DimPlot(sc_pre, reduction="pca", dims=c(1, 2)) + NoLegend()
p2 <- DimPlot(sc_pre, reduction="pca", dims=c(2, 3)) + NoLegend()
p3 <- DimPlot(sc_pre, reduction="pca", dims=c(1, 3)) + NoLegend()
p <- p1 | p2 | p3
ggsave(file.path(result_dir, "1.2.4.4_PCA.pdf"), p, width = 12, heigh = 5)

ElbowPlot(sc_pre, reduction="pca", ndims=50)
VizDimLoadings(sc_pre, reduction="pca", dims=1:4, nfeatures=10, balanced=TRUE)

###### 1.2.4.5 Clustering and Visualization ####
sc_pre <- FindNeighbors(sc_pre, dims=1:30) #reduction="pca",k.param=20
qsave(sc_pre, file.path(result_dir, "1.2.4.5_seurat_pre.qs"))

resolution_check <- FindClusters(sc_pre, resolution = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)) #algorithm="leiden", random.seed=1, method="igraph"
tree <- clustree(resolution_check)
ggsave(file.path(result_dir, "1.2.4.5_ClusterTree.pdf"), tree , width = 11, heigh = 7)

sc_final <- FindClusters(sc_pre, resolution=0.4) #random.seed=1, method="igraph"

sc_final <- RunUMAP(sc_final, reduction="pca", reduction.name="umap", dims=1:30, return.model=TRUE)
sc_final$seurat_clusters_normal <- sc_final$seurat_clusters
qsave(sc_final, file.path(result_dir, "1.2.4.5_seurat_final_res0.4.qs"))

p1 <- DimPlot(sc_final, reduction="umap", label=TRUE) + NoLegend()
p2 <- ImageDimPlot(sc_final, border.size = NA, dark.background = FALSE, size = 0.5, group.by = "seurat_clusters_normal") 
p3 <- p1 | p2
ggsave(file.path(result_dir, "1.2.4.5_Clustering_UMAP.pdf"), p3 , width = 7, heigh = 5)



##### 1.2.5 ROI layer #####
###### 1.2.5.1 subset ROI ####
#subset specific part of the image
#open Xenium Explorer
#mark with the free selection tool your region of interest (roi) 
#name selection
#three points next to the name --> export selection coordinates
#change to cells --> three dots next to cell stats --> export cell stats as csv --> delete first two rows --> this is your input (has cell_ID)
roi_df <- read.csv() #cell stats from Xenium Exlorer
roi_cells <- roi_df$Cell.ID
sc_roi <- subset(sc_final, cells = roi_cells)
qsave(sc_roi, file.path(result_dir, "1.2.5_sc_roi_layer.qs"))

###### 1.2.5.2 re-clustering of subset ####
sc_roi$seurat_clusters <- NULL
DefaultAssay(sc_roi) <- "Xenium" 
### seurat ###
sc_roi <- NormalizeData(sc_roi, normalization.method = "LogNormalize", scale.factor = 10000)
sc_roi <- FindVariableFeatures(sc_roi, nfeatures = 3000)
sc_roi <- ScaleData(sc_roi, features = VariableFeatures(sc_roi))
sc_roi <- RunPCA(sc_roi, features = VariableFeatures(sc_roi), npcs = 30)
ElbowPlot(sc_roi, ndims = 30) 
sc_roi <- FindNeighbors(sc_roi, dims = 1:20) 
res_check <- FindClusters(sc_roi, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8))
tree <- clustree(res_check)
ggsave(file.path(result_dir, "1.2.5.1_roi_layer_newClustering_tree.pdf"), tree , width = 8, heigh = 10)

sc_roi <- FindClusters(sc_roi, resolution = 0.6)
sc_roi$seurat_clusters_normal <- sc_roi$seurat_clusters
qsave(sc_roi, file.path(result_dir, "1.2.5.1_sc_roi_layer_newClustering_res0.6.qs"))

### BANKSY ###
sc_roi <- RunBanksy(sc_roi,
                    assay="Xenium", assay_name="BANKSY",
                    features="variable",
                    lambda=0.2, k_geom=15,
                    verbose=TRUE)
DefaultAssay(sc_roi) = "BANKSY"
sc_roi <- RunPCA(sc_roi, reduction.name="banksy_pca", features=rownames(sc_roi[["BANKSY"]]), npcs=30, nfeatures.print=5)
sc_roi <- FindNeighbors(sc_roi, reduction="banksy_pca", dims=1:20)
resolution_check <- FindClusters(sc_roi, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) 
tree <- clustree(resolution_check)
tree
ggsave(file.path(result_dir, "1.2.5.1_roi_layer_newClusterTree_Banksy.pdf"), tree , width = 11, heigh = 7)
sc_roi <- FindClusters(sc_roi, resolution=0.4)
sc_roi <- RunUMAP(sc_roi, reduction="banksy_pca", reduction.name="banksy_umap", dims=1:30, return.model=TRUE)
sc_roi$banksy_cluster <- sc_roi$BANKSY_snn_res.0.4
qsave(sc_roi, file.path(result_dir, "1.2.5.1_sc_roi_layer_newClustering_res0.4.qs"))


#### 1.3 prepare Ctrl sample ####
##### 1.3.1 load dataset ####
raw_data_dir <- #define directory
setwd(raw_data_dir)
NBB_Ctrl <- #file_name
sc <- LoadXenium(NBB_Ctrl, fov = "fov")
sc$orig.ident <- "dura"
sc$sample_ID <- "NBB_Ctrl"
Idents(sc) <- "dura"

qsave(sc, file.path(result_dir, "1.2_Xenium_NBB_Ctrl_seurat.qs"))

result_dir <- #define directory
setwd(result_dir)

##### 1.3.2 quality control ####
###### 1.3.2.1 Number of Counts per Bin ####
bins = cells
#The total number of UMIs (counts) per bin gives a sense of overall transcript abundance.
p1 <- VlnPlot(sc, features="nCount_Xenium", pt.size=0, layer="counts")
p2 <- ImageFeaturePlot(sc, features="nCount_Xenium", dark.background = FALSE, border.size = NA, size = 1)
p3 <- p1 | p2

ggsave(file.path(result_dir, "1.3.2.1_QC_number_of_counts_per_cells.pdf"), p3, width = 7, heigh = 5)
summary(sc$nCount_Xenium)
#Look for bins with very few or very many counts. Bins with <30 counts or >2000 counts are suspicious and flagged for removal.
sc$filtered_by_counts = ifelse(sc$nCount_Xenium < 30 | sc$nCount_Xenium > 2000, "filtered-out", "keep")
p <- ImageDimPlot(sc, group.by="filtered_by_counts", border.size = NA, dark.background = FALSE, size = 0.5)
ggsave(file.path(result_dir, "1.3.2.1_QC_filter_by_counts.pdf"), p, width = 5, heigh = 5)


###### 1.3.2.2 Number of Genes per Bin ####
p1 <- VlnPlot(sc, features="nFeature_Xenium", pt.size=0, layer="counts") + NoLegend()
p2 <- ImageFeaturePlot(sc, features="nFeature_Xenium", dark.background = FALSE, border.size = NA, size = 1)
p3 <- p1 | p2
ggsave(file.path(result_dir, "1.3.2.2_QC_number_of_genes_per_cells.pdf"), p3, width = 7, heigh = 5)

summary(sc$nFeature_Xenium)

sc$filtered_by_genes = ifelse(sc$nFeature_Xenium < 5 | sc$nFeature_Xenium > 1000, "filtered-out", "keep")
p1 <- ImageDimPlot(sc, group.by="filtered_by_genes", border.size = NA, dark.background = FALSE, size = 0.5)
p2 <- VlnPlot(sc, features="nFeature_Xenium", pt.size=0, layer="counts", group.by = "filtered_by_genes") + NoLegend()
p3 <- p1 | p2
ggsave(file.path(result_dir, "1.3.2.2_QC_filter_by_genes.pdf"), p3, width = 7, heigh = 5)

##### 1.3.3 Filtering ====
dim(sc)
###### 1.3.3.1 features (and nCounts) #####
sc_filtered <- subset(sc, nFeature_Xenium < 5 | nFeature_Xenium > 1000, invert=TRUE) %>% suppressWarnings()
#sc_filtered <- subset(sc_filtered, nCount_Xenium < 30 | nCount_Xenium > 2000, invert=TRUE) %>% suppressWarnings()
qsave(sc_filtered, file.path(result_dir, "1.3.3.1_Xenium_NBB_Ctrl_seurat_filtered.qs"))

p1 <- ImageDimPlot(sc, border.size = NA, dark.background = FALSE, size = 0.5) + ggtitle("not filtered")
p2 <- ImageDimPlot(sc_filtered, border.size = NA, dark.background = FALSE, size = 0.5) + ggtitle("filtered")
p3 <- p1 | p2

###### 1.3.3.2 genes ####
#Genes expressed in only a few bins can add noise. First, we count how many bins each gene is detected in:
counts <- GetAssayData(sc, layer="counts")
genes_expressed <- rowSums(counts > 0)
head(genes_expressed, 5)
summary(genes_expressed)
#Keep genes expressed in >10 bins:
genes_expressed <- genes_expressed[genes_expressed > 10]
length(genes_expressed)
#133 (= no genes need to be filtered out)
ggsave(file.path(result_dir, "1.3.3.2_filtered.pdf"), p3, width = 7, heigh = 5)

##### 1.3.4 Pre-processing ====
###### 1.3.4.1 Normalization ####
sc_pre <- NormalizeData(sc_filtered, normalization.method="LogNormalize", scale.factor=10000)

##### 1.3.4.2 Identifying Highly Variable Genes (HVGs) ####
sc_pre <- FindVariableFeatures(sc_pre, nfeatures=3000)

top20 <- VariableFeatures(sc_pre) %>% head(20)
p <- VariableFeaturePlot(sc_pre, selection.method="vst")
p <- LabelPoints(plot=p, points=top20, repel=TRUE)
p
ggsave(file.path(result_dir, "1.3.4.2_VariableFeatures_top.pdf"), p, width = 7, heigh = 5)

##### 1.3.4.3 Data scaling ####
sc_pre <- ScaleData(sc_pre, features=VariableFeatures(sc_pre))

##### 1.3.4.4 Dimensionality reduction ####
sc_pre <- RunPCA(sc_pre, reduction.name="pca", features=VariableFeatures(sc_pre), npcs=50, nfeatures.print=5)

p1 <- DimPlot(sc_pre, reduction="pca", dims=c(1, 2)) + NoLegend()
p2 <- DimPlot(sc_pre, reduction="pca", dims=c(2, 3)) + NoLegend()
p3 <- DimPlot(sc_pre, reduction="pca", dims=c(1, 3)) + NoLegend()
p <- p1 | p2 | p3
ggsave(file.path(result_dir, "1.3.4.4_PCA.pdf"), p, width = 12, heigh = 5)

ElbowPlot(sc_pre, reduction="pca", ndims=50)
VizDimLoadings(sc_pre, reduction="pca", dims=1:4, nfeatures=10, balanced=TRUE)

##### 1.3.4.5 Clustering and Visualization ####
sc_pre <- FindNeighbors(sc_pre, dims=1:30) #reduction="pca",k.param=20
qsave(sc_pre, file.path(result_dir, "1.3.4.5_seurat_pre.qs"))
 
resolution_check <- FindClusters(sc_pre, resolution = c(0.05, 0.1, 0.2, 0.3)) #algorithm="leiden", random.seed=1, method="igraph"
tree <- clustree(resolution_check)
ggsave(file.path(result_dir, "1.3.4.5_ClusterTree.pdf"), tree , width = 11, heigh = 7)

sc_final <- FindClusters(sc_pre, resolution=0.2) #random.seed=1, method="igraph"

sc_final <- RunUMAP(sc_final, reduction="pca", reduction.name="umap", dims=1:30, return.model=TRUE)
sc_final$seurat_clusters_normal <- sc_final$seurat_clusters
qsave(sc_final, file.path(result_dir, "1.3.4.5_seurat_final_res0.2.qs"))

p1 <- DimPlot(sc_final, reduction="umap", label=TRUE) + NoLegend()
p2 <- ImageDimPlot(sc_final, border.size = NA, dark.background = FALSE, size = 0.5, group.by = "seurat_clusters_normal") 
p3 <- p1 | p2
ggsave(file.path(result_dir, "1.3.4.5_Clustering_UMAP.pdf"), p3 , width = 7, heigh = 5)

#### 1.3.5 ROI layer #### 
###### 1.3.5.1 subset ROI ####
#subset specific part of the image
#open Xenium Explorer
#mark with the free selection tool your region of interest (roi) 
#name selection
#three points next to the name --> export selection coordinates
#change to cells --> three dots next to cell stats --> export cell stats as csv --> this is your input (has cell_ID)
roi_df <- read.csv() #cell stats from Xenium Exlorer
roi_cells <- roi_df$Cell.ID
sc_roi <- subset(sc_final, cells = roi_cells)
qsave(sc_roi, file.path(result_dir, "1.3.5.1_sc_roi_layer.qs"))

##### 1.3.5.2 re-clustering of subset ####
sc_roi$seurat_clusters <- NULL
sc_roi$banksy_cluster <- NULL
sc_roi$seurat_clusters_normal <- NULL
sc_roi$BANKSY_snn_res.0.45 <- NULL

DefaultAssay(sc_roi) <- "Xenium" 
### seurat ###
sc_roi <- NormalizeData(sc_roi, normalization.method = "LogNormalize", scale.factor = 10000)
sc_roi <- FindVariableFeatures(sc_roi, nfeatures = 3000)
sc_roi <- ScaleData(sc_roi, features = VariableFeatures(sc_roi))
sc_roi <- RunPCA(sc_roi, features = VariableFeatures(sc_roi), npcs = 30)
ElbowPlot(sc_roi, ndims = 30) 
sc_roi <- FindNeighbors(sc_roi, dims = 1:20) 
res_check <- FindClusters(sc_roi, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8))
tree <- clustree(res_check)
ggsave(file.path(result_dir, "1.3.5.2_roi_layer_newClustering_tree.pdf"), tree , width = 8, heigh = 10)

sc_roi <- FindClusters(sc_roi, resolution = 0.6)
sc_roi$seurat_clusters_normal <- sc_roi$seurat_clusters
qsave(sc_roi, file.path(result_dir, "1.3.5.2_sc_roi_layer_newClustering_res0.6.qs"))

### BANKSY ###
sc_roi <- RunBanksy(sc_roi,
                    assay="Xenium", assay_name="BANKSY",
                    features="variable",
                    lambda=0.2, k_geom=15,
                    verbose=TRUE)
DefaultAssay(sc_roi) = "BANKSY"
sc_roi <- RunPCA(sc_roi, reduction.name="banksy_pca", features=rownames(sc_roi[["BANKSY"]]), npcs=30, nfeatures.print=5)
sc_roi <- FindNeighbors(sc_roi, reduction="banksy_pca", dims=1:20)
resolution_check <- FindClusters(sc_roi, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) 
tree <- clustree(resolution_check)
tree
ggsave(file.path(result_dir, "1.3.5.2_roi_layer_newClusterTree_Banksy.pdf"), tree , width = 11, heigh = 7)
sc_roi <- FindClusters(sc_roi, resolution=0.8)
sc_roi <- RunUMAP(sc_roi, reduction="banksy_pca", reduction.name="banksy_umap", dims=1:30, return.model=TRUE)
sc_roi$banksy_cluster <- sc_roi$BANKSY_snn_res.0.8
qsave(sc_roi, file.path(result_dir, "1.3.5.2_sc_roi_layer_newClustering_res0.8.qs"))


#### 1.4 load MS and Ctrl dataset ####
#two samples (MS and Ctrl) and only load roi
sc_MS <- qread("1.2.5.1_sc_roi_layer_newClustering_res0.4.qs")
sc_Ctrl <- qread("1.3.5.2_sc_roi_layer_newClustering_res0.8.qs")



#==== 2.0 Integration MS and Ctrl roi ====
DefaultAssay(sc_Ctrl) <- "Xenium"
DefaultAssay(sc_MS) <- "Xenium"
samples_list <- list(sc_Ctrl, sc_MS)

#find inegration features
features <- SelectIntegrationFeatures(object.list = samples_list, nfeatures = 3000)
#prepate for integration
samples_list <- lapply(samples_list, function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
#find integration anchors
anchors <- FindIntegrationAnchors(object.list = samples_list, anchor.features = features, dims = 1:30)
#integrate data 
sc_combined <- IntegrateData(anchorset = anchors, dims = 1:30)
qsave(sc_combined, file.path(result_dir, "2.0_integrated.qs"))


#==== 3.0 Clustering + Marker genes ====
#### 3.1 Clustering ####
DefaultAssay(sc_combined) <- "integrated"
sc_combined <- ScaleData(sc_combined, verbose = FALSE)
sc_combined <- RunPCA(sc_combined, npcs = 30)
sc_combined <- FindNeighbors(sc_combined, dims = 1:30)

res_check <- FindClusters(sc_combined, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2))
tree <- clustree(res_check)
ggsave(file.path(result_dir, "3.1_combined_seurat_clustering_res_tree.pdf"), tree , width = 8, heigh = 10)

sc_combined <- FindClusters(sc_combined, resolution = 0.8)
sc_combined <- RunUMAP(sc_combined, reduction="pca", reduction.name="umap", dims=1:30, return.model=TRUE)
sc_combined$seurat_clusters_normal <- sc_combined$seurat_clusters
qsave(sc_combined, file.path(result_dir, "3.1_combined_seurat_clustering_res0.8.qs"))

p <- DimPlot(sc_combined, reduction="umap", group.by="sample_ID")
ggsave(file.path(result_dir, "3.1_combined_DimPlot_seurat.pdf"), p , width = 5, heigh = 5)

#visualization not possible on sc_combined because in the integrated object only one image can be stored --> subset samples and map integrated_clusters back
sc_Ctrl_subset <- subset(sc_combined, subset = sample_ID == "NBB_Ctrl")
sc_MS_subset <- subset(sc_combined, subset = sample_ID == "Patho_MS")
sc_Ctrl_subset$cluster_integrated <- sc_combined$seurat_clusters_normal[colnames(sc_Ctrl_subset)]
sc_MS_subset$cluster_integrated   <- sc_combined$seurat_clusters_normal[colnames(sc_MS_subset)]

#plots
Idents(sc_combined) <- "seurat_clusters_normal"
p1 <- DimPlot(sc_combined, reduction="umap", label=TRUE, group.by="seurat_clusters_normal") + NoLegend() + ggtitle("integrated")
p2 <- DimPlot(sc_combined, reduction="umap", group.by="sample_ID") + ggtitle("grouped.by = sample")
p3 <- DimPlot(sc_combined, reduction="umap", split.by="sample_ID") + ggtitle("split.by = sample")
p4 <- DimPlot(sc_Ctrl_subset, group.by = "cluster_integrated") + ggtitle("Ctrl")
p5 <- DimPlot(sc_MS_subset, group.by = "cluster_integrated") + ggtitle("MS")
p6 <-  p1 | p2 | p3 | p4 | p5
ggsave(file.path(result_dir, "3.1_DimPlot_integrated_seurat_clusters.pdf"), p6 , width = 20, heigh = 7)
ggsave(file.path(result_dir, "3.1_DimPlot_integrated_seurat_clusters_split.by_sample.pdf"), p3 , width = 10, heigh = 7)
ggsave(file.path(result_dir, "3.1_DimPlot_integrated_seurat_clusters_group.by_sample.pdf"), p2 , width = 6, heigh = 7)

p1 <- DimPlot(sc_combined, reduction="umap", label=TRUE, group.by="seurat_clusters_normal") + NoLegend() + ggtitle("integrated")
p2 <- ImageDimPlot(sc_Ctrl_subset, border.size = NA, dark.background = FALSE, size = 1) + ggtitle("Ctrl")
p3 <- ImageDimPlot(sc_Ctrl_subset, border.size = NA, dark.background = FALSE, size = 1, split.by = "cluster_integrated") + ggtitle("Ctrl")
p4 <- ImageDimPlot(sc_MS_subset, border.size = NA, dark.background = FALSE, size = 1) + ggtitle("MS")
p5 <- ImageDimPlot(sc_MS_subset, border.size = NA, dark.background = FALSE, size = 1, split.by = "cluster_integrated") + ggtitle("MS")
p6 <- p1 | p2 / p3 | p4 / p5
ggsave(file.path(result_dir, "3.1_Spatial_integrated_seurat_clusters.pdf"), p6 , width = 12, heigh = 7)

qsave(sc_combined, file.path(result_dir, "3.1_integrated_seurat_clustering.qs"))
qsave(sc_Ctrl_subset, file.path(result_dir, "3.1_integrated_Ctrl_subset_seurat_clustering.qs"))
qsave(sc_MS_subset, file.path(result_dir, "3.1_integrated_MS_subset_seurat_clustering.qs"))

#### 3.2 marker genes ####
Idents(sc_combined) <- "seurat_clusters_normal"
sc_combined <- JoinLayers(sc_combined)
markers <- FindAllMarkers(object=sc_combined,
                          test.use="wilcox",
                          layer="data",
                          only.pos=FALSE,
                          max.cells.per.ident=1000)
markers_sorted <- markers %>%
  arrange(cluster, p_val_adj, desc(avg_log2FC))
marker_list <- split(markers_sorted, markers_sorted$cluster)
write.xlsx(marker_list, file.path(result_dir, "3.2_integrated_seurat_Cluster_Markers_seurat.xlsx"))


#==== 4.0 subset Fibro and recluster ====
#### 4.1 subset Fibro ####
Idents(Xenium_combined) <- "seurat_clusters_normal"
Xenium_combined_Fibro <- subset(Xenium_combined, seurat_clusters_normal %in% c("0", "1", "5")) #Fibro cluster

#### 4.2 reclustering ####
DefaultAssay(Xenium_combined_Fibro) <- "integrated"
Xenium_combined_Fibro <- ScaleData(Xenium_combined_Fibro, verbose = FALSE)
Xenium_combined_Fibro <- RunPCA(Xenium_combined_Fibro, npcs = 30)
Xenium_combined_Fibro <- FindNeighbors(Xenium_combined_Fibro, dims = 1:30)

res_check <- FindClusters(Xenium_combined_Fibro, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0))
tree <- clustree(res_check)
ggsave(file.path(result_dir, "4.2_Xenium_combined_Fibro_subset_restree.pdf"), tree , width = 8, heigh = 10)

Xenium_combined_Fibro <- FindClusters(Xenium_combined_Fibro, resolution = 1.4)
Xenium_combined_Fibro <- RunUMAP(Xenium_combined_Fibro, reduction="pca", reduction.name="umap", dims=1:30, return.model=TRUE)
Xenium_combined_Fibro$seurat_clusters_normal <- Xenium_combined_Fibro$seurat_clusters
qsave(Xenium_combined_Fibro, file.path(result_dir, "4.2_Xenium_combined_Fibro_subset_res1.4.qs"))

p <- DimPlot(Xenium_combined_Fibro, reduction="umap", group.by="seurat_clusters_normal")
ggsave(file.path(result_dir, "4.2_Xenium_combined_Fibro_DimPlot_res1.4.pdf"), p , width = 5, heigh = 5)

#### 4.3 visualization ####
#visualization not possible on sc_combined because in the integrated object only one image can be stored --> subset samples and map integrted_clusters back
sc_Ctrl_subset_Fibro <- subset(Xenium_combined_Fibro, subset = sample_ID == "NBB_Ctrl")
sc_MS_subset_Fibro <- subset(Xenium_combined_Fibro, subset = sample_ID == "Patho_MS")
sc_Ctrl_subset_Fibro$cluster_integrated <- Xenium_combined_Fibro$seurat_clusters_normal[colnames(sc_Ctrl_subset_Fibro)]
sc_MS_subset_Fibro$cluster_integrated   <-Xenium_combined_Fibro$seurat_clusters_normal[colnames(sc_MS_subset_Fibro)]

#plots
Idents(Xenium_combined_Fibro) <- "seurat_clusters_normal"
p1 <- DimPlot(Xenium_combined_Fibro, reduction="umap", label=TRUE, group.by="seurat_clusters_normal") + NoLegend() + ggtitle("integrated")
p2 <- DimPlot(Xenium_combined_Fibro, reduction="umap", group.by="sample_ID") + ggtitle("grouped.by = sample")
p3 <- DimPlot(Xenium_combined_Fibro, reduction="umap", split.by="sample_ID") + ggtitle("split.by = sample")
p4 <- DimPlot(sc_Ctrl_subset_Fibro, group.by = "cluster_integrated") + ggtitle("Ctrl")
p5 <- DimPlot(sc_MS_subset_Fibro, group.by = "cluster_integrated") + ggtitle("MS")
p6 <-  p1 | p2 | p3 | p4 | p5
ggsave(file.path(result_dir, "4.3_DimPlot_integrated_seurat_clusters_subset_Fibro_res1.4.pdf"), p6 , width = 20, heigh = 7)
ggsave(file.path(result_dir, "4.3_DimPlot_integrated_seurat_clusters_subset_Fibro_split.by_sample_res1.4.pdf"), p3 , width = 10, heigh = 7)
ggsave(file.path(result_dir, "4.3_DimPlot_integrated_seurat_clusters_subset_Fibro_group.by_sample_res1.4.pdf"), p2 , width = 6, heigh = 7)

p1 <- DimPlot(Xenium_combined_Fibro, reduction="umap", label=TRUE, group.by="seurat_clusters_normal") + NoLegend() + ggtitle("integrated")
p2 <- ImageDimPlot(sc_Ctrl_subset_Fibro, border.size = NA, dark.background = FALSE, size = 1) + ggtitle("Ctrl")
p3 <- ImageDimPlot(sc_Ctrl_subset_Fibro, border.size = NA, dark.background = FALSE, size = 1, split.by = "cluster_integrated") + ggtitle("Ctrl")
p4 <- ImageDimPlot(sc_MS_subset_Fibro, border.size = NA, dark.background = FALSE, size = 1) + ggtitle("MS")
p5 <- ImageDimPlot(sc_MS_subset_Fibro, border.size = NA, dark.background = FALSE, size = 1, split.by = "cluster_integrated") + ggtitle("MS")
p6 <- p1 | p2 / p3 | p4 / p5
ggsave(file.path(result_dir, "4.3_Spatial_integrated_seurat_clusters_subset_Fibro_res1.4.pdf"), p6 , width = 12, heigh = 7)

qsave(Xenium_combined_Fibro, file.path(result_dir, "4.3_integrated_seurat_subset_Fibro_clustering.qs"))
qsave(sc_Ctrl_subset_Fibro, file.path(result_dir, "4.3_integrated_Ctrl_subset_strucFbro_seurat_clustering.qs"))
qsave(sc_MS_subset_Fibro, file.path(result_dir, "4.3_integrated_MS_subset_struFibro_seurat_clustering.qs"))

#### 4.4 marker genes ####
Idents(Xenium_combined_Fibro) <- "seurat_clusters_normal"
markers <- FindAllMarkers(object=Xenium_combined_Fibro,
                          test.use="wilcox",
                          layer="data",
                          only.pos=FALSE,
                          max.cells.per.ident=1000)
markers_sorted <- markers %>%
  arrange(cluster, p_val_adj, desc(avg_log2FC))
marker_list <- split(markers_sorted, markers_sorted$cluster)
write.xlsx(marker_list, file.path(result_dir, "4.4_integrated_seurat_Cluster_Markers_subset_Fibro_res1.4.xlsx"))

#visualization #
# Path to your Excel file
file_path <- file.path(result_dir, "4.4_integrated_seurat_Cluster_Markers_subset_Fibro_res1.4.xlsx")
# List all sheets
sheets <- getSheetNames(file_path)
print(sheets)
# Read all sheets into a list
marker_list_read <- lapply(sheets, function(s){
  read.xlsx(file_path, sheet = s)
})
# Assign names to the list (cluster names)
names(marker_list_read) <- sheets

markers_combined <- do.call(rbind, marker_list_read)

markers_top = markers_combined %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::arrange(p_val_adj, avg_log2FC) %>% 
  dplyr::slice_head(n=10) %>% 
  dplyr::ungroup() %>% 
  as.data.frame()

#plot top_markers in Xenium in spatialfeaturePlots 
p <- ImageFeaturePlot(sc_Ctrl_subset_Fibro, features=markers_top$gene %>% unique(), border.size = NA, dark.background = FALSE, size = 2) 
ggsave(file.path(result_dir, "4.4_spFP_marker_genes_Xenium_subset_Fibro_in_spCtrl.pdf"), p , width = 20, heigh = 20, limitsize = F)
p <- ImageFeaturePlot(sc_MS_subset_Fibro, features=markers_top$gene %>% unique(), border.size = NA, dark.background = FALSE, size = 2) 
ggsave(file.path(result_dir, "4.4_spFP_marker_genes_Xenium_subset_Fibro_in_spMS.pdf"), p , width = 20, heigh = 20, limitsize = F)


#plot top_markers in Xenium in Xenium_combined
p <- FeaturePlot(Xenium_combined_Fibro, features=markers_top$gene %>% unique()) 
ggsave(file.path(result_dir, "4.4_FP_marker_genes_Xenium_subset_Fibro_in_Xenium.pdf"), p , width = 20, heigh = 40, limitsize = F)
p <- DotPlot(Xenium_combined_Fibro, features=markers_top$gene %>% unique()) +
  viridis::scale_color_viridis() + 
  ylab("Cluster") + xlab("") + 
  guides(size=guide_legend(order=1, title="Percent expressed"), color=guide_colorbar(title="Scaled Expression")) + 
  ggtitle("Top markers per cluster (expression scaled)") + RotatedAxis()
ggsave(file.path(result_dir, "4.4_DP_marker_genes_Xenium_subset_Fibro_in_Xenium.pdf"), p , width = 15, heigh = 10)


#===== 5.0 Fibros focus ====
###### 5.1 duraFibro + vascFibro ####
#after subsetting Fibro cluster from integrated dataset and performing reclustering, there is one cluster that does not contain fibroblasts = Cluster 6 ---> remove for Fibro focus
Xenium_duraFibro <- subset(Xenium_combined_Fibro, idents = c("0", "1", "2", "3", "4", "5", "7")) #6 = non Fibros in meningeal layer

cluster_color <- c("0" = "#66c69b", #meningeal 
                   "1" = "#ee756d", #periost
                   "2" = "#87be4d", #bordFibro
                   "3" = "#d0e7fc", #vascFibro 
                   "4" = "#d0e7fc", #vascFibro 
                   "5" = "#d0e7fc", #vascFibro 
                   "7" = "#d49005") #inner meningeal


sc_Ctrl_duraFibro <- subset(Xenium_duraFibro, subset = sample_ID == "NBB_Ctrl")
sc_MS_duraFibro <- subset(Xenium_duraFibro, subset = sample_ID == "Patho_MS")
sc_Ctrl_duraFibro$cluster_integrated <- Xenium_duraFibro$seurat_clusters_normal[colnames(sc_Ctrl_duraFibro)]
sc_MS_duraFibro$cluster_integrated   <-Xenium_duraFibro$seurat_clusters_normal[colnames(sc_MS_duraFibro)]

p1 <- DimPlot(Xenium_duraFibro, reduction="umap", label=TRUE, group.by="seurat_clusters_normal", cols = cluster_color) + NoLegend() + ggtitle("integrated")
p2 <- ImageDimPlot(sc_Ctrl_duraFibro, border.size = NA, dark.background = FALSE, size = 1.5, cols = cluster_color) + ggtitle("Ctrl")
p3 <- ImageDimPlot(sc_Ctrl_duraFibro, border.size = NA, dark.background = FALSE, size = 1, split.by = "cluster_integrated", cols = cluster_color) + ggtitle("Ctrl")
p4 <- ImageDimPlot(sc_MS_duraFibro, border.size = NA, dark.background = FALSE, size = 1.5, cols = cluster_color) + ggtitle("MS")
p5 <- ImageDimPlot(sc_MS_duraFibro, border.size = NA, dark.background = FALSE, size = 1, split.by = "cluster_integrated", cols = cluster_color) + ggtitle("MS")
p6 <- p1 | p2 / p3 | p4 / p5
ggsave(file.path(result_dir, "5.1_Spatial_ImageDimPlot_integrated_seurat_clusters_duraFibro+vascFibro.pdf"), p6 , width = 12, heigh = 7)

##### 5.2 only duraFibro #####
Xenium_duraFibro <- subset(Xenium_combined_strucFibro, idents = c("0", "1", "2", "7")) #6 = non Fibros in meningeal layer, #3-5 vascFibro

cluster_color <- c("0" = "#66c69b", #meningeal duraFibro1-3
                   "1" = "#ee756d", #durFibro3, periost
                   "2" = "#87be4d", #bordFibro
                   "7" = "#d49005") #duraFibro4, inner meningeal



sc_Ctrl_duraFibro <- subset(Xenium_duraFibro, subset = sample_ID == "NBB_Ctrl")
sc_MS_duraFibro <- subset(Xenium_duraFibro, subset = sample_ID == "Patho_MS")
sc_Ctrl_duraFibro$cluster_integrated <- Xenium_duraFibro$seurat_clusters_normal[colnames(sc_Ctrl_duraFibro)]
sc_MS_duraFibro$cluster_integrated   <-Xenium_duraFibro$seurat_clusters_normal[colnames(sc_MS_duraFibro)]

p1 <- DimPlot(Xenium_duraFibro, reduction="umap", label=TRUE, group.by="seurat_clusters_normal", cols = cluster_color) + NoLegend() + ggtitle("integrated")
p2 <- ImageDimPlot(sc_Ctrl_duraFibro, border.size = NA, dark.background = FALSE, size = 1.5, cols = cluster_color) + ggtitle("Ctrl")
p3 <- ImageDimPlot(sc_Ctrl_duraFibro, border.size = NA, dark.background = FALSE, size = 1, split.by = "cluster_integrated", cols = cluster_color) + ggtitle("Ctrl")
p4 <- ImageDimPlot(sc_MS_duraFibro, border.size = NA, dark.background = FALSE, size = 1.5, cols = cluster_color) + ggtitle("MS")
p5 <- ImageDimPlot(sc_MS_duraFibro, border.size = NA, dark.background = FALSE, size = 1, split.by = "cluster_integrated", cols = cluster_color) + ggtitle("MS")
p6 <- p1 | p2 / p3 | p4 / p5
ggsave(file.path(result_dir, "5.2_Spatial_ImageDimPlot_integrated_seurat_clusters_duraFibro.pdf"), p6 , width = 12, heigh = 7)
p2 <- ImageDimPlot(sc_Ctrl_duraFibro, border.size = NA, dark.background = FALSE, size = 15, cols = cluster_color)
ggsave(file.path(result_dir, "5.2_Spatial_ImageDimPlot_integrated_seurat_clusters_duraFibro_onlyCtrl.pdf"), p2 , width = 48, heigh = 28)
p4 <- ImageDimPlot(sc_MS_duraFibro, border.size = NA, dark.background = FALSE, size = 15, cols = cluster_color) 
ggsave(file.path(result_dir, "5.2_Spatial_ImageDimPlot_integrated_seurat_clusters_duraFibro_onlyMS.pdf"), p4 , width = 48, heigh = 28)

#export data for Xenium explorer
# Build a metadata table
build_cell_table <- function(seurat_obj, cluster_color) {
  data.frame(
    cell_id = colnames(seurat_obj),
    group   = seurat_obj$cluster_integrated,
    color   = cluster_color[as.character(seurat_obj$cluster_integrated)]
  )
}
cell_table_ctrl <- build_cell_table(sc_Ctrl_duraFibro, cluster_color)
cell_table_ms   <- build_cell_table(sc_MS_duraFibro,  cluster_color)
write_csv(cell_table_ctrl, "5.2_duraFibro_cells_CTRL.csv")
write_csv(cell_table_ms,   "5.2_duraFibro_cells_MS.csv")

