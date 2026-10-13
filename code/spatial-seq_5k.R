#===== 1.0 set up =====

#### 1.1 environment ####
#define directory
#raw_data_dir <- 
#result_dir <- 

#packages
source(file.path(result_dir, "0_global.R"))

# Colors
cal_pal50 = c("#Fa1a8e", "#009B7D", "#ff9933", "#7083a4", "#ffce45", "#015e05", 
              "#fedf25", "#d2b48c", "#bb55e1", "#6ec5ab", "#5d57af", "#143341", 
              "#761445", "#d65b5a", "#94043a", "#e7a6cd", "#204519", "#87afd1", 
              "#9b9a9d", "#f95b70", "#83c874", "#808080", "#452b5b", "#ecb100", 
              "#f46124", "#525252", "#4c84a3", "#00bfff", "#01b4c6", "#174d79", 
              "#a6a0f2", "#76facc", "#8491be", "#a32a2f", "#1c8859", "#2cc012", 
              "#35782b", "#9c6755", "#3b3960", "#eeb0a1", "#3e1e13", "#0064c3", 
              "#d81e4a", "#74646c", "#f675da", "#ffce45", "#ec7014", "#e50000", 
              "#000000", "#a4527c", "#041859")

cluster_colors <- c(
  "0"  = "#Fa1a8e",
  "1"  = "#009B7D",
  "2"  = "#ff9933",
  "3"  = "#7083a4",
  "4"  = "#ffce45",
  "5"  = "#015e05",
  "6"  = "#fedf25",
  "7"  = "#d2b48c",
  "8"  = "#bb55e1",
  "9"  = "#6ec5ab",
  "10" = "#5d57af",
  "11" = "#143341",
  "12" = "#761445",
  "13" = "#d65b5a",
  "14" = "#94043a",
  "15" = "#e7a6cd",
  "16" = "#204519",
  "17" = "#87afd1",
  "18" = "#9b9a9d",
  "19" = "#f95b70",
  "20" = "#83c874",
  "21" = "#808080",
  "22" = "#452b5b",
  "23" = "#ecb100",
  "24" = "#f46124",
  "25" = "#525252",
  "26" = "#4c84a3",
  "27" = "#00bfff",
  "28" = "#01b4c6"
)

cluster_colors_major_groups_mix <- c( "ABC" = "#9584bd",
                                      "arach_bordFibro" = "#87be4d",
                                      "duraFibro" = "#f2cf22", 
                                      "vascStroma" = "#00b8df", 
                                      "Endo" = "#27a9b5", 
                                      "Myeloid" = "#e04c73",
                                      "Bc" = "#bf72ad",
                                      "Tc" = "#dd4e10", 
                                      "mix" = "grey")


cluster_colors_snRNA <- c(
  "ABC" = "#9584bd",
  "arach_bordFibro" = "#87be4d",
  "duraFibro1" = "#38b083",
  "duraFibro2" = "#5c96d0",
  "duraFibro3" = "#ee756e",
  "duraFibro4" = "#d49005",
  "duraFibro5" = "#828cc5",
  "duraFibro6" = "#ec67a2",
  "duraFibro7" = "#3bab4b",
  "duraFibro8" = "#d88356",
  "vascFibro" = "#9ea71e",
  "mixVascStroma" = "#a27eb7",
  "Peri" = "#06b8e0",
  "vSMC1" = "#bb4d63",
  "vSMC2" = "#8e7eb9",
  "mixEndo" = "#b078b1",
  "artEndo" = "#3caa34",
  "capEndo" = "#29aab5",
  "venEndo" = "#957e25",
  "Mast" = "#e26aa5",
  "inflamMono1" = "#36b39c",
  "inflamMono2" = "#75ae2b",
  "inflamMono3" = "#12b3e8",
  "Granulo" = "#cf6ea8",
  "mixMyeloid" = "#32b6b3",
  "IFN-CAM" = "#ec67a2",
  "CAM" = "#e38803",
  "Bc" = "#bf72ad",
  "Tc" = "#dd4e10",
  "unknown" = "lightgrey")


# Set default colours for ggplot2
options(ggplot2.discrete.colour=cal_pal50)
options(ggplot2.discrete.fill=cal_pal50)
options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill="viridis")
options(ggplot2.binned.colour="viridis")
options(ggplot2.binned.fill="viridis")



#### 1.2 loading data ####
#load sample 
#sample <- LoadXenium(file.path(raw_data_dir,"output-XETxy"), fov = "fov")

#### 1.3 define samples ####
sample_name <- c("Ctrl_2019_042", "Ctrl_2022_112", "Ctrl_2023_021", "MS_2022_082", "MS_4074", "MS_2019_047", "Ctrl_2019_061", "Ctrl_2023_085", "Ctrl_2023_015", "MS_2022_103", "MS_2022_110", "MS_3867")

# Create a list of Seurat objects for each sample, with respective sample names
sample_list <- list(Ctrl_2019_042, Ctrl_2022_112, Ctrl_2023_021, MS_2022_082, MS_4074, MS_2019_047, Ctrl_2019_061, Ctrl_2023_085, Ctrl_2023_015, MS_2022_103, MS_2022_110, MS_3867
)|> setNames(sample_name)

qsave(sample_list, file.path(result_dir, "1.3_sample_list.qs"))


#==== 2.0 Quality control ====
#In spatial transcriptomics (and any sequencing-based assay), not all bins are equally informative. 
#Quality Control (QC) is used to identify and remove low-quality bins from the dataset before further analysis. 
# **Low-quality bins** might be the result of:  
#* Problems during reverse transcription, PCR amplification or library preparation  
#* Cell damage during the experiment, e.g. when the tissue already contains a lot of dead or dying cells  
#* Very few transcripts overall, e.g. in empty regions  
#* Unusually high numbers of transcripts, possibly doublets or noise  

# *Why does this matter?**  
#Low-quality bins can distort the results. For instance, they might cluster separately and be falsely interpreted as novel biological states or cell types. Performing QC early helps ensure the reliability of downstream steps like normalization, clustering, and cell type annotation.

## QC metrics at a glance
#We will use three common metrics to flag low-quality bins:
# **Number of counts (UMIs):** Measures total expression in a bin. Too low = empty; too high = potential artifacts.
# **Number of detected genes:** Captures the diversity of expression in a bin. Low gene counts might mean poor capture.
# **Mitochondrial gene content:** High levels often indicate stressed or dying cells.
#These metrics will be assessed both statistically and spatially to get an idea of the sample quality and to identify thresholds for filtering.


##### 2.1  Calculate QC metrics for each sample ####
#define filtering
for(sample_name in names(sample_list)) {
  
  message("Processing: ", sample_name)
  
  sc <- sample_list[[sample_name]]
  
  # Store sample name
  sc$sample <- sample_name
  
  # Count filtering
  sc$filtered_by_counts <- ifelse(
    sc$nCount_Xenium < 5 |  #needs to be adjusted for the sample type
      sc$nCount_Xenium > 500, #needs to be adjusted for the sample type
    "filtered-out",
    "keep"
  )
  
  # Gene filtering
  sc$filtered_by_genes <- ifelse(
    sc$nFeature_Xenium < 5 | #needs to be adjusted for the sample type
      sc$nFeature_Xenium > 500, #needs to be adjusted for the sample type
    "filtered-out",
    "keep"
  )
  
  sample_list[[sample_name]] <- sc
}

# Merge objects for combined violin plots
combined <- merge(
  x = sample_list[[1]],
  y = sample_list[-1]
)


#### 2.2 Number of Counts per Cell ####
# Combined violin plot across all samples

p_counts <- VlnPlot(combined, features = "nCount_Xenium", group.by = "sample", pt.size = 0) + RotatedAxis()

ggsave(file.path(result_dir, "2.2.0_QC_number_of_counts_per_cells.pdf"), p_counts, width = 12, height = 5)

p_counts <- VlnPlot(combined,
                    features = "nCount_Xenium",
                    group.by = "sample",
                    pt.size = 0.1
) +
  RotatedAxis()

ggsave(file.path(result_dir, "2.2.0_QC_number_of_counts_per_cells_individualCells.pdf"), p_counts, width = 12, height = 5)

###### 2.2.1 Spatial distribution #####
count_maps <- lapply(names(sample_list), function(sample_name){
  
  ImageDimPlot(
    sample_list[[sample_name]],
    group.by = "filtered_by_counts",
    cols = c(
      "keep" = "red",
      "filtered-out" = "grey"
    ),
    border.size = NA,
    dark.background = FALSE,
    size = 0.5
  ) +
    ggtitle(sample_name)
  
})

p_count_maps <- wrap_plots(count_maps, ncol = 3)

ggsave(file.path(result_dir, "2.2.1_QC_filter_by_counts_spatial.pdf"), p_count_maps, width = 16, height = 12)

###### 2.2.2 Filtered Vln plots #####
filtered_cells <- subset(combined, subset = filtered_by_counts == "filtered-out")

p_counts_filtered <- VlnPlot(filtered_cells,
                             features = "nCount_Xenium",
                             group.by = "sample",
                             pt.size = 0
) +
  RotatedAxis() +
  ggtitle("Filtered-out cells")

kept_cells <- subset(combined,
                     subset = filtered_by_counts == "keep")

p_counts_keep <- VlnPlot(kept_cells,
                         features = "nCount_Xenium",
                         group.by = "sample",
                         pt.size = 0
) +
  RotatedAxis() +
  ggtitle("Kept cells")

p_counts_summary <- p_counts_filtered / p_counts_keep

ggsave(
  file.path(result_dir, "2.2.2_QC_counts_filtering_summary.pdf"),
  p_counts_summary,
  width = 12,
  height = 8
)

p_counts_filtered <- VlnPlot(filtered_cells,
                             features = "nCount_Xenium",
                             group.by = "sample",
                             pt.size = 0.1
) +
  RotatedAxis() +
  ggtitle("Filtered-out cells")

p_counts_keep <- VlnPlot(kept_cells,
                         features = "nCount_Xenium",
                         group.by = "sample",
                         pt.size = 0.1
) +
  RotatedAxis() +
  ggtitle("Kept cells")

p_counts_summary <- p_counts_filtered / p_counts_keep

ggsave(file.path(result_dir, "2.2.2_QC_counts_filtering_summary_individualCells.pdf"), p_counts_summary, width = 12, height = 8)

##### 2.3 Number of Genes per Cell ####
p_genes <- VlnPlot(combined,
                   features = "nFeature_Xenium",
                   group.by = "sample",
                   pt.size = 0
) +
  RotatedAxis()

ggsave(file.path(result_dir, "2.3.0_QC_number_of_genes_per_cells.pdf"), p_genes, width = 12, height = 5)

p_genes <- VlnPlot(combined,
                   features = "nFeature_Xenium",
                   group.by = "sample",
                   pt.size = 0.1
) +
  RotatedAxis()

ggsave(file.path(result_dir, "2.3.0_QC_number_of_genes_per_cells_imdividualCells.pdf"), p_genes, width = 12, height = 5)

###### 2.3.1 Spatial distribution #####
gene_maps <- lapply(names(sample_list), function(sample_name){
  
  ImageDimPlot(
    sample_list[[sample_name]],
    group.by = "filtered_by_genes",
    cols = c(
      "keep" = "red",
      "filtered-out" = "grey"
    ),
    border.size = NA,
    dark.background = FALSE,
    size = 0.5
  ) +
    ggtitle(sample_name)
  
})

p_gene_maps <- wrap_plots(gene_maps, ncol = 3)

ggsave(file.path(result_dir, "2.3.1_QC_filter_by_genes_spatial.pdf"), p_gene_maps, width = 16, height = 12)

###### 2.3.2 Filtered Vln plots #####
filtered_genes <- subset(combined,
                         subset = filtered_by_genes == "filtered-out")

p_genes_filtered <- VlnPlot(filtered_genes,
                            features = "nFeature_Xenium",
                            group.by = "sample",
                            pt.size = 0
) +
  RotatedAxis() +
  ggtitle("Filtered-out cells")

kept_genes <- subset(combined,
                     subset = filtered_by_genes == "keep")

p_genes_keep <- VlnPlot(kept_genes,
                        features = "nFeature_Xenium",
                        group.by = "sample",
                        pt.size = 0
) +
  RotatedAxis() +
  ggtitle("Kept cells")

p_genes_summary <- p_genes_filtered / p_genes_keep

ggsave(
  file.path(result_dir, "2.3.2_QC_genes_filtering_summary.pdf"),
  p_genes_summary,
  width = 12,
  height = 8
)

p_genes_filtered <- VlnPlot(filtered_genes,
                            features = "nFeature_Xenium",
                            group.by = "sample",
                            pt.size = 0.1
) +
  RotatedAxis() +
  ggtitle("Filtered-out cells")

p_genes_keep <- VlnPlot(kept_genes,
                        features = "nFeature_Xenium",
                        group.by = "sample",
                        pt.size = 0.1
) +
  RotatedAxis() +
  ggtitle("Kept cells")

p_genes_summary <- p_genes_filtered / p_genes_keep

ggsave(file.path(result_dir, "2.3.2_QC_genes_filtering_summary_individualCells.pdf"), p_genes_summary, width = 12, height = 8)


#===== 3.0 Filtering ====
# transfer names
for(sample_name in names(sample_list)) {
  
  sample_list[[sample_name]]$sample <- sample_name
  
}

sample_list_before <- sample_list

#save sample_list before filtering
qsave(sample_list_before, file.path(result_dir, "3.0_samples_list_before_filtering.qs"))


##### 3.1 filtering + after filtering summary ####
# Excel workbook for gene summaries
wb <- createWorkbook()

# Filter each sample
for(sample_name in names(sample_list_before)){
  
  message("Filtering: ", sample_name)
  
  sc <- sample_list_before[[sample_name]]
  sc_filtered <- subset(
    sc,
    subset =
      nCount_Xenium >= 5 &
      nCount_Xenium <= 500 &
      nFeature_Xenium >= 5 &
      nFeature_Xenium <= 500
  ) %>%
    suppressWarnings()
  
  sc_filtered$sample <- sample_name
  
  filtered_list[[sample_name]] <- sc_filtered
  
  # Gene expression summary
  counts <- GetAssayData(
    sc,
    layer = "counts"
  )
  
  genes_expressed <- rowSums(counts > 0)
  
  gene_df <- data.frame(
    gene = names(genes_expressed),
    n_bins = genes_expressed
  )
  
  addWorksheet(
    wb,
    sheetName = substr(sample_name, 1, 31)
  )
  
  writeData(
    wb,
    sheet = substr(sample_name, 1, 31),
    x = gene_df
  )
  
  message(
    sample_name,
    ": ",
    length(genes_expressed[genes_expressed > 10]),
    " genes expressed in >10 cells"
  )
  
}

# Save filtered objects
sample_list <- filtered_list

qsave(sample_list, file.path(result_dir, "3.1_samples_list_after_filtering.qs"))

# Save gene expression workbook
saveWorkbook(wb, file.path(result_dir, "3.1_gene_expression_summary.xlsx"), overwrite = TRUE)

# Filtering summary
filtering_summary <- data.frame(
  sample = names(sample_list_before),
  cells_before = sapply(sample_list_before, ncol),
  cells_after = sapply(sample_list, ncol)
)

filtering_summary$percent_kept <-
  round(
    100 *
      filtering_summary$cells_after /
      filtering_summary$cells_before,
    1
  )

write.csv2(filtering_summary, file.path(result_dir, "3.1_filtering_summary.csv"), row.names = FALSE)


#### 3.2 after QC filtering plots ####
# Combined before/after filtering plots
filter_maps <- lapply(names(sample_list_before), function(sample_name){
  
  sc_before <- sample_list_before[[sample_name]]
  sc_after  <- sample_list[[sample_name]]
  
  p1 <- ImageDimPlot(sc_before, border.size = NA, dark.background = FALSE, size = 0.5) +
    ggtitle(paste0(sample_name, "\nBefore (n=", ncol(sc_before),")"))
  
  p2 <- ImageDimPlot(sc_after, border.size = NA, dark.background = FALSE, size = 0.5) +
    ggtitle(paste0(sample_name, "\nAfter (n=", ncol(sc_after),")"))
  
  p1 | p2
  
})

p_filter_maps <- wrap_plots(filter_maps, ncol = 2)

ggsave(file.path(result_dir, "3.2_ImageDimPlot_before_after_filtering.pdf"), p_filter_maps, width = 16, height = 24)

# Combined filtered object
combined_filtered <- merge(
  x = sample_list[[1]],
  y = sample_list[-1]
)

# Counts after filtering
p_counts <- VlnPlot(combined_filtered, features = "nCount_Xenium", group.by = "sample", pt.size = 0.05) + RotatedAxis()

ggsave(file.path(result_dir, "3.2_counts_after_filtering.pdf"), p_counts, width = 12, height = 5)

# Features after filtering
p_features <- VlnPlot(
  combined_filtered,
  features = "nFeature_Xenium",
  group.by = "sample",
  pt.size = 0.05
) +
  RotatedAxis()

ggsave(file.path(result_dir, "3.3_features_after_filtering.pdf"), p_features, width = 12, height = 5)


#===== 4.0 Integration =====
##### 4.1 Normalization #######
#SCT
for (sample_name in names(sample_list)) {
  
  sample_list[[sample_name]] <- SCTransform(
    sample_list[[sample_name]],
    assay = "Xenium",
    vst.flavor = "v2", 
    verbose = T
  )
  
}

qsave(sample_list, file.path(result_dir, "4.1_sample_list_SCT.qs"))

#### 4.2 Integration ####
features <- SelectIntegrationFeatures(object.list = sample_list, nfeatures = 3000)

#problem not all 3000 features are present in all samples --> need to filter the common ones first
#common features
common_scale_features <- Reduce(intersect, lapply(sample_list, function(x) {
  rownames(GetAssayData(x, assay = "SCT", layer = "scale.data"))}))
length(common_scale_features)

#restrict integration features
features <- intersect(features, common_scale_features)
length(features)

#PreSCt Integration 
sample_list <- PrepSCTIntegration(
  object.list = sample_list,
  anchor.features = features
)


#Run PCA for each sample
sample_list <- lapply(sample_list, function(x) {
  x <- RunPCA(x, features = features)
  return(x)
})


anchors <- FindIntegrationAnchors(
  object.list = sample_list,
  normalization.method = "SCT",
  anchor.features = features,
  dims = 1:30
)

Fibro_seurat <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT"
)

qsave(sc_combined, file.path(result_dir, "4.2_samples_combined.qs"))

#### 4.3 Clustering ####
###### 4.3.1 seurat clustering ####
DefaultAssay(sc_combined) <- "integrated"
sc_combined <- RunPCA(sc_combined, npcs = 30)
ElbowPlot(sc_combined, ndims = 50)
sc_combined <- FindNeighbors(sc_combined, dims = 1:20)

res_check <- FindClusters(sc_combined, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2))
tree <- clustree(res_check)
ggsave(file.path(result_dir, "4.3_combined_seurat_clustering_res_tree.pdf"), tree , width = 8, heigh = 10)

sc_combined <- FindClusters(sc_combined, resolution = 0.4)
sc_combined <- RunUMAP(sc_combined, reduction="pca", reduction.name="umap", dims=1:30, return.model=TRUE)
sc_combined$seurat_clusters_normal <- sc_combined$seurat_clusters
qsave(sc_combined, file.path(result_dir, "4.3.1_combined_seurat_clustering_res0.4.qs"))

Idents(sc_combined) <- "seurat_clusters_normal"
p <- DimPlot(sc_combined)
ggsave(file.path(result_dir, "4.3_DimPlot_seurat_cluster.pdf"), p , width = 5, heigh = 5)

p <- DimPlot(sc_combined, reduction="umap", split.by="sample")
ggsave(file.path(result_dir, "4.3.1_DimPlot_seurat_cluster_splitby_sample.pdf"), p , width = 20, heigh = 5)

#=== 5.0 Clustering + Marker genes =====
##### 5.1 Marker genes ######
Idents(sc_combined) <- "seurat_clusters_normal"
markers_seurat <- FindAllMarkers(object=sc_combined,
                                 test.use="wilcox",
                                 layer="data",
                                 only.pos=FALSE,
                                 max.cells.per.ident=1000)

markers_sorted <- markers_seurat %>%
  arrange(cluster, p_val_adj, desc(avg_log2FC))
marker_list <- split(markers_sorted, markers_sorted$cluster)
write.xlsx(marker_list, file.path(result_dir, "5.1_Cluster_Markers_seurat_res0.4.xlsx"))

#### 5.2 visualize marker genes ####
markers_top = markers_seurat %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::arrange(p_val_adj, avg_log2FC) %>% 
  dplyr::slice_head(n=3) %>% 
  dplyr::ungroup() %>% 
  as.data.frame()

p <- DotPlot(sc_combined, features=markers_top$gene %>% unique()) +
  viridis::scale_color_viridis() + 
  ylab("Cluster") + xlab("") + 
  guides(size=guide_legend(order=1, title="Percent expressed"), color=guide_colorbar(title="Scaled Expression")) + 
  ggtitle("Top markers per cluster (expression scaled)") + RotatedAxis()
ggsave(file.path(result_dir, "5.2_DP_top_marker_genes_seurat_res0.4.pdf"), p , width = 15, heigh = 7)

#### 5.3 marker genes snRNAseq ####
Idents(sc_combined) <- "seurat_clusters_normal"
DefaultAssay(sc_combined) <- "SCT"
genes<- read_csv(file.path(result_dir, "0_marker_genes_snRNAseq.txt"), col_names = FALSE)
p <- DotPlot(sc_combined, features=genes) +
  viridis::scale_color_viridis() + 
  ylab("Cluster") + xlab("") + 
  guides(size=guide_legend(order=1, title="Percent expressed"), color=guide_colorbar(title="Scaled Expression")) + 
  ggtitle("marker genes snRNA-seq") + RotatedAxis()
ggsave(file.path(result_dir, "5.3_DP_marker_genes_snRNAseq_seurat_res0.4.pdf"), p , width = 7, heigh = 7)
genes<- read_csv(file.path(result_dir, "0_marker_genes_snRNAseq_Fibro.txt"), col_names = FALSE)
p <- DotPlot(sc_combined, features=genes) +
  viridis::scale_color_viridis() + 
  ylab("Cluster") + xlab("") + 
  guides(size=guide_legend(order=1, title="Percent expressed"), color=guide_colorbar(title="Scaled Expression")) + 
  ggtitle("marker genes snRNA-seq") + RotatedAxis()
ggsave(file.path(result_dir, "5.3_DP_marker_genes_snRNAseq_Fibro_seurat_res0.4.pdf"), p , width = 7, heigh = 7)


#==== 6.0 map snRNA-seq annotations onto spatial Xenium data ====
#snRNA-seq data as ref
#snRNA_obj <-  our human snRNA-seq dataset

DefaultAssay(sc_combined) <- "SCT"
DefaultAssay(snRNA_obj) <- "SCT"
anchors <- FindTransferAnchors(reference = snRNA_obj, query = sc_combined, dims = 1:30, reference.reduction = "pca",normalization.method = "SCT", reference.assay = "integrated")
predictions <- TransferData(anchorset = anchors, refdata = snRNA_obj$annotated_manuscript,dims = 1:30)
sc_combined_pred <- AddMetaData(sc_combined, metadata = predictions)
#if UMAP emddings are not save, you need to reun it here RUNUMAP()
sc_combined_pred <- MapQuery(anchorset = anchors, reference = snRNA_obj, query = sc_combined_pred,refdata = list(celltype = "annotated_manuscript"), reference.reduction = "pca", reduction.model = "umap")
cluster <- as.character(unique(snRNA_obj@meta.data$annotated_manuscript))

#colors snRNA-seq dataset
cluster_colors <- c(
  "ABC" = "#9584bd",
  "arach_bordFibro" = "#87be4d",
  "duraFibro1" = "#38b083",
  "duraFibro2" = "#5c96d0",
  "duraFibro3" = "#ee756e",
  "duraFibro4" = "#d49005",
  "duraFibro5" = "#828cc5",
  "duraFibro6" = "#ec67a2",
  "duraFibro7" = "#3bab4b",
  "duraFibro8" = "#d88356",
  "vascFibro" = "#9ea71e",
  "mixVascStroma" = "#a27eb7",
  "Peri" = "#06b8e0",
  "vSMC1" = "#bb4d63",
  "vSMC2" = "#8e7eb9",
  "mixEndo" = "#b078b1",
  "artEndo" = "#3caa34",
  "capEndo" = "#29aab5",
  "venEndo" = "#957e25",
  "Mast" = "#e26aa5",
  "inflamMono1" = "#36b39c",
  "inflamMono2" = "#75ae2b",
  "inflamMono3" = "#12b3e8",
  "Granulo" = "#cf6ea8",
  "mixMyeloid" = "#32b6b3",
  "IFN-CAM" = "#ec67a2",
  "CAM" = "#e38803",
  "Bc" = "#bf72ad",
  "Tc" = "#dd4e10",
  "unknown" = "lightgrey")

# transfer labels 
#label low-confidence cells as "unknown" 
sc_combined$snRNA_annotation <- ifelse(
  sc_combined_pred$prediction.score.max < 0.3,
  "unknown",
  sc_combined_pred$predicted.id
)

#major groups
cluster_to_group <- c(
  "ABC" = "ABC",
  "arach_bordFibro" = "arach_bordFibro",
  "duraFibro1" = "duraFibro",
  "duraFibro2" = "duraFibro",
  "duraFibro3" = "duraFibro",
  "duraFibro4" = "duraFibro",
  "duraFibro5" = "duraFibro",
  "duraFibro6" = "duraFibro",
  "duraFibro7" = "duraFibro",
  "duraFibro8" = "duraFibro",
  "vascFibro" = "vascStroma",
  "mixVascStroma" = "vascStroma",
  "Peri" = "vascStroma",
  "vSMC1" = "vascStroma",
  "vSMC2" = "vascStroma",
  "mixEndo" = "Endo",
  "artEndo" = "Endo",
  "capEndo" = "Endo",
  "venEndo" = "Endo",
  "Mast" = "Myeloid",
  "inflamMono1" = "Myeloid",
  "inflamMono2" = "Myeloid",
  "inflamMono3" = "Myeloid",
  "Granulo" = "Myeloid",
  "mixMyeloid" = "Myeloid",
  "IFN-CAM" = "Myeloid",
  "CAM" = "Myeloid",
  "Bc" = "Bc",
  "Tc" = "Tc",
  "unknown" = "unknown")

cluster_colors_major_groups <- c( "ABC" = "#9584bd",
                                  "arach_bordFibro" = "#87be4d",
                                  "duraFibro" = "#f2cf22", 
                                  "vascStroma" = "#00b8df", 
                                  "Endo" = "#27a9b5", 
                                  "Myeloid" = "#e04c73",
                                  "Bc" = "#bf72ad",
                                  "Tc" = "#dd4e10")


major_group_vector <- cluster_to_group[sc_combined$snRNA_annotation]
names(major_group_vector) <- colnames(sc_combined)
sc_combined <- AddMetaData(sc_combined, metadata = major_group_vector, col.name = "snRNA_major_group")


cluster_to_group_mix <- c(
  "ABC" = "ABC",
  "arach_bordFibro" = "arach_bordFibro",
  "duraFibro1" = "duraFibro",
  "duraFibro2" = "duraFibro",
  "duraFibro3" = "duraFibro",
  "duraFibro4" = "duraFibro",
  "duraFibro5" = "duraFibro",
  "duraFibro6" = "duraFibro",
  "duraFibro7" = "duraFibro",
  "duraFibro8" = "duraFibro",
  "vascFibro" = "vascStroma",
  "mixVascStroma" = "mix",
  "Peri" = "vascStroma",
  "vSMC1" = "vascStroma",
  "vSMC2" = "vascStroma",
  "mixEndo" = "mix",
  "artEndo" = "Endo",
  "capEndo" = "Endo",
  "venEndo" = "Endo",
  "Mast" = "Myeloid",
  "inflamMono1" = "Myeloid",
  "inflamMono2" = "Myeloid",
  "inflamMono3" = "Myeloid",
  "Granulo" = "Myeloid",
  "mixMyeloid" = "mix",
  "IFN-CAM" = "Myeloid",
  "CAM" = "Myeloid",
  "Bc" = "Bc",
  "Tc" = "Tc",
  "unknown" = "unknown")

cluster_colors_major_groups_mix <- c( "ABC" = "#9584bd",
                                      "arach_bordFibro" = "#87be4d",
                                      "duraFibro" = "#f2cf22", 
                                      "vascStroma" = "#00b8df", 
                                      "Endo" = "#27a9b5", 
                                      "Myeloid" = "#e04c73",
                                      "Bc" = "#bf72ad",
                                      "Tc" = "#dd4e10", 
                                      "mix" = "grey")

major_group_vector <- cluster_to_group_mix[sc_combined$snRNA_annotation]
names(major_group_vector) <- colnames(sc_combined)
sc_combined <- AddMetaData(sc_combined, metadata = major_group_vector, col.name = "snRNA_major_group_mix")

qsave(sc_combined_pred, file.path(result_dir, "6.0_snRNA_reference_prediction.qs"))

#### 6.1 reference comparison seurat clusters #####
#plotting
DefaultAssay(sc_combined_pred) <- "SCT"
Idents(sc_combined_pred) <- "seurat_clusters_normal"
DefaultAssay(sc_combined) <- "SCT"
Idents(sc_combined) <- "seurat_clusters_normal"

p1 <- DimPlot(snRNA_obj, reduction = "umap", group.by = "annotated_manuscript", label = TRUE, label.size = 3,
              repel = TRUE, cols = cluster_colors) + NoLegend() + ggtitle("snRNA-seq") + ylim(-19,19)
p2 <- DimPlot(sc_combined_pred, reduction = "ref.umap", group.by = "seurat_clusters_normal", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Xenium") +ylim(-19,19)
p <- p1 + p2
ggsave(file.path(result_dir, "6.1_DimPlot_ref_Plot_seurat.pdf"), p , width = 7, heigh = 5)

p2 <- FeaturePlot(sc_combined_pred, features = "predicted.celltype.score", reduction = "ref.umap", label = T, repel = T, pt.size = 0.8) +
  ylim(-19,19) + ggtitle("Xenium")
p <- p1 + p2
ggsave(file.path(result_dir, "6.1_DimPlot_ref_Plot_pred_seurat.pdf"), p , width = 7, heigh = 5)

p <- DimPlot(sc_combined_pred, group.by = "predicted.celltype", label = T, cols = cluster_colors)
ggsave(file.path(result_dir, "6.1_DimPlot_snRNAseq_predicted.celltypes_seurat.pdf"), p , width = 10, heigh = 8)

p <- FeaturePlot(sc_combined_pred, features = cluster,  reduction = "umap", cols = c("lightgrey", "darkred"), ncol = 3) & theme(plot.title = element_text(size = 10))
ggsave(file.path(result_dir, "6.1_FeaturePlot_cohort1_predicted.celltypes_seurat.pdf"), p , width = 50, heigh = 100, limitsize = F)

p <- DimPlot(sc_combined, group.by = "snRNA_annotation", cols = cluster_colors)
ggsave(file.path(result_dir, "6.1_DimPlot_seurat_transferred_labels_snRNA.pdf"), p , width = 7, heigh = 5)

p <- DimPlot(sc_combined, split.by = "snRNA_annotation", cols = cluster_colors)
ggsave(file.path(result_dir, "6.1_DimPlot_seurat_splitby_transferred_labels_snRNA.pdf"), p , width = 30, heigh = 5)

p <- DimPlot(sc_combined, group.by = "snRNA_major_group", cols = cluster_colors_major_groups)
ggsave(file.path(result_dir, "6.1_DimPlot_seurat_transferred_labels_snRNA_major_groups.pdf"), p , width = 5, heigh = 5)

p <- DimPlot(sc_combined, group.by = "snRNA_major_group_mix", cols = cluster_colors_major_groups_mix)
ggsave(file.path(result_dir, "6.1_DimPlot_seurat_transferred_labels_snRNA_major_groups_mix.pdf"), p , width = 5, heigh = 5)

#### 6.2 marker genes snRNA_major_group_mix ####
Idents(sc_combined) <- "snRNA_major_group_mix"
sc_combined <- PrepSCTFindMarkers(sc_combined)
markers_seurat <- FindAllMarkers(object=sc_combined,
                                 test.use="wilcox",
                                 layer="data",
                                 only.pos=FALSE,
                                 max.cells.per.ident=1000)

markers_sorted <- markers_seurat %>%
  arrange(cluster, p_val_adj, desc(avg_log2FC))
marker_list <- split(markers_sorted, markers_sorted$cluster)
write.xlsx(marker_list, file.path(result_dir, "6.2_Cluster_Markers_major_group_mix.xlsx"))


markers_top = markers_seurat %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::arrange(p_val_adj, avg_log2FC) %>% 
  dplyr::slice_head(n=3) %>% 
  dplyr::ungroup() %>% 
  as.data.frame()

p <- DotPlot(sc_combined, features=markers_top$gene %>% unique()) +
  viridis::scale_color_viridis() + 
  ylab("Cluster") + xlab("") + 
  guides(size=guide_legend(order=1, title="Percent expressed"), color=guide_colorbar(title="Scaled Expression")) + 
  ggtitle("Top markers per cluster (expression scaled)") + RotatedAxis()
ggsave(file.path(result_dir, "6.2_DP_top_marker_genes_major_group_mix.pdf"), p , width = 10, heigh = 5)

#===== 7.0 Spatial visualization ====
split_objects <- SplitObject(sc_combined, split.by = "sample")
qsave(split_objects, file.path(result_dir, "7.0_split_objects_samples.qs"))


#seurat clusters
cluster_maps <- lapply(names(split_objects), function(sample_name) {
  
  ImageDimPlot(
    split_objects[[sample_name]],
    group.by = "seurat_clusters_normal",
    border.size = NA,
    dark.background = FALSE,
    size = 0.5
  ) +
    ggtitle(sample_name)
  
})

combined_plot <- wrap_plots(cluster_maps) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right"
  )

ggsave(
  file.path(result_dir, "7.0_ImageDimPlot_all_sc_combined_clusters_res0.4.pdf"),
  combined_plot,
  width = 12,
  height = 8
)


#snRNA-seq clustering all
cluster_maps <- lapply(names(split_objects), function(sample_name) {
  
  ImageDimPlot(
    split_objects[[sample_name]],
    group.by = "snRNA_annotation",
    border.size = NA,
    dark.background = FALSE,
    size = 0.5, 
    cols = cluster_colors
  ) +
    ggtitle(sample_name)
  
})

combined_plot <- wrap_plots(cluster_maps) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right"
  )

ggsave(
  file.path(result_dir, "7.0_ImageDimPlot_all_samples_snRNA_annotation.pdf"),
  combined_plot,
  width = 12,
  height = 8
)


#snRNA-seq clustering major group mix
cluster_maps <- lapply(names(split_objects), function(sample_name) {
  
  ImageDimPlot(
    split_objects[[sample_name]],
    group.by = "snRNA_major_group_mix",
    border.size = NA,
    dark.background = FALSE,
    size = 0.5, 
    cols = cluster_colors_major_groups_mix
  ) +
    ggtitle(sample_name)
  
})

combined_plot <- wrap_plots(cluster_maps) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right"
  )

ggsave(
  file.path(result_dir, "7.0_ImageDimPlot_all_samples_snRNA_major_group_mix.pdf"),
  combined_plot,
  width = 12,
  height = 8
)

#split.by seurat_clusters_normal
DefaultAssay(sc_combined) <- "SCT"
Idents(sc_combined) <- "seurat_clusters_normal"
cluster_levels <- levels(Idents(sc_combined))
cluster_colors_seurat <- Seurat:::DiscretePalette(length(cluster_levels))
names(cluster_colors_seurat) <- cluster_levels


cluster_maps <- lapply(names(split_objects), function(sample_name) {
  
  obj <- split_objects[[sample_name]]
  
  celltypes <- sort(unique(obj$seurat_clusters_normal))
  
  plots <- lapply(celltypes, function(ct) {
    
    obj$highlight <- ifelse(
      obj$seurat_clusters_normal == ct,
      as.character(ct),
      "Other"
    )
    
    ImageDimPlot(
      obj,
      group.by = "highlight",
      cols = c(
        "Other" = "grey85",
        cluster_colors_seurat[as.character(ct)]
      ),
      border.size = NA,
      dark.background = FALSE,
      size = 1
    ) +
      ggtitle(paste(sample_name, "-", ct))
    
  })
  
  wrap_plots(plots, ncol = 3)
})

combined_plot <- wrap_plots(cluster_maps) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right"
  )

ggsave(
  file.path(
    result_dir,
    "7.0_ImageDimPlot_all_sc_combined_cluster_highlighted.pdf"
  ),
  combined_plot,
  width = 45,
  height = 50, 
  limitsize = F
)

#==== 8.0 Meta data ====
patient_info <- openxlsx::read.xlsx(file.path(patient_dir, "xy")) #use supplement table
#add meta.data
meta <- sc_combined@meta.data
# Add cell name as a new column for subsequent reassignment
meta$cell_id <- rownames(meta)
# Merge `information from `meta` and `patient_info`
meta_new <- left_join(meta, patient_info[, c("sample_id", "patient_id", "sample", "dura_region", "condition", "MS.form", "FFPE", "in.snRNAseq.cohort", "sex", "age", "pmd_h", "years.after.diagnosis", "EDSS", "modified.rankin.scale")], by = c("sample" = "sample"))
rownames(meta_new) <- meta$cell_id
#add to seurat_object
sc_combined@meta.data <- meta_new
sc_combined@meta.data$sample_id <- as.factor(sc_combined@meta.data$sample_id)
sc_combined@meta.data$patient_id <- as.factor(sc_combined@meta.data$patient_id)
sc_combined@meta.data$sample <- as.factor(sc_combined@meta.data$sample)
sc_combined@meta.data$dura_region <- as.factor(sc_combined@meta.data$dura_region)
sc_combined@meta.data$condition <- as.factor(sc_combined@meta.data$condition)
sc_combined@meta.data$MS.form <- as.factor(sc_combined@meta.data$MS.form)
sc_combined@meta.data$FFPE <- as.factor(sc_combined@meta.data$FFPE)
sc_combined@meta.data$in.snRNAseq.cohort <- as.factor(sc_combined@meta.data$in.snRNAseq.cohort)
sc_combined@meta.data$sex <- as.factor(sc_combined@meta.data$sex)
sc_combined@meta.data$age <- as.factor(sc_combined@meta.data$age)
sc_combined@meta.data$pmd_h <- as.factor(sc_combined@meta.data$pmd_h)
sc_combined@meta.data$years.after.diagnosis <- as.factor(sc_combined@meta.data$years.after.diagnosis)
sc_combined@meta.data$EDSS <- as.factor(sc_combined@meta.data$EDSS)
sc_combined@meta.data$modified.rankin.scale <- as.factor(sc_combined@meta.data$modified.rankin.scale)

rm(meta, meta_new, patient_info)

qsave(sc_combined, file.path(result_dir, "8.0_sc_combined_meta_data_downstream.qs"))


#===== 9.0 Fibro (snRNA_Major_group_mix)  res 0.115 =====
DefaultAssay(sc_combined) <- "SCT"
Idents(sc_combined) <- "snRNA_major_group_mix"
Fibro_subset <- subset(sc_combined, idents = c("duraFibro", "ABC", "arach_bordFibro"))
DimPlot(Fibro_subset)

#reclustering
Fibro_subset@assays$integrated <- NULL
Fibro_subset@assays$SCT <- NULL

# split by sample
DefaultAssay(Fibro_subset) <- "Xenium"
Fibro_list <- SplitObject(Fibro_subset, split.by = "sample")

# SCTransform per sample
Fibro_list <- lapply(Fibro_list, function(x) {
  SCTransform(x, vst.flavor = "v2", assay = "Xenium", verbose = TRUE)
})

# re-integration and re-clustering 
features <- SelectIntegrationFeatures(object.list = Fibro_list, nfeatures = 3000)

common_scale_features <- Reduce(intersect, lapply(Fibro_list, function(x) {
  rownames(GetAssayData(x, assay = "SCT", layer = "scale.data"))}))
length(common_scale_features)
features <- intersect(features, common_scale_features)
length(features)

Fibro_list <- PrepSCTIntegration(Fibro_list, anchor.features = features) 

anchors <- FindIntegrationAnchors(Fibro_list, normalization.method = "SCT",anchor.features = features)

qsave(anchors, file.path(result_dir, "9.0_integration_anchors_Fibros.qs"))
#I ran it on the cluster
Fibro_combined <- IntegrateData(anchors, normalization.method = "SCT") |>
  RunPCA(verbose = T)
Fibro_combined <- RunPCA(Fibro_combined, verbose = TRUE)
ElbowPlot(Fibro_combined, ndims = 50)
Fibro_combined <- RunUMAP(Fibro_combined, dims = 1:20) #change dims
Fibro_combined <- FindNeighbors(Fibro_combined, dims = 1:20) #change dims
qsave(Fibro_combined, file.path(result_dir, "9.0_Fibro_subset_integrated_before_clustering.qs"))

Fibro_combined <- qread(file.path(result_dir, "9.0_Fibro_subset_integrated_before_clustering.qs"))

Fibro_seurat <- FindClusters(Fibro_combined, resolution = 0.115) #change resolution
Fibro_seurat <- RunUMAP(Fibro_seurat, dims = 1:30, return.model = T) #mreturnmodle true to save umap projection
p <- DimPlot(Fibro_seurat, label = T) 
ggsave(file.path(result_dir, "9.0_Fibro_DimPlot_res0.115.pdf"), p, width = 6, height = 6)

p <- DimPlot(Fibro_seurat, split.by = "seurat_clusters")
ggsave(file.path(result_dir, "9.0_Fibro_DimPlot_cluster_res0.115.pdf"), p, width =15, height = 6)

qsave(Fibro_seurat, file.path(result_dir, "13.0_Fibro_subset_integrated_clustering_res0.115.qs"))


split_objects <- SplitObject(Fibro_seurat, split.by = "sample")
qsave(split_objects, file.path(result_dir, "9.0_Fibro_seurat-duraFibro_wo_mix_split_objects_samples.qs"))

#seurat clusters
cluster_maps <- lapply(names(split_objects), function(sample_name) {
  
  ImageDimPlot(
    split_objects[[sample_name]],
    group.by = "seurat_clusters",
    border.size = NA,
    dark.background = FALSE,
    size = 0.5
  ) +
    ggtitle(sample_name)
  
})

combined_plot <- wrap_plots(cluster_maps) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right"
  )

ggsave(
  file.path(result_dir, "9.0_ImageDimPlot_duraFibro_wo_mix_seurat_clusters_0.115.pdf"),
  combined_plot,
  width = 12,
  height = 8
)

cluster_maps <- lapply(names(split_objects), function(sample_name) {
  
  ImageDimPlot(
    split_objects[[sample_name]],
    split.by = "seurat_clusters",
    border.size = NA,
    dark.background = FALSE,
    size = 0.5
  ) +
    ggtitle(sample_name)
  
})

combined_plot <- wrap_plots(cluster_maps) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right"
  )

ggsave(
  file.path(result_dir, "9.0_ImageDimPlot_duraFibro_wo_mix_splitbyseurat_clusters_0.115.pdf"),
  combined_plot,
  width = 12,
  height = 12
)

#split.by 
DefaultAssay(Fibro_seurat) <- "SCT"
Idents(Fibro_seurat) <- "seurat_clusters"
cluster_levels <- levels(Idents(Fibro_seurat))
cluster_colors_seurat <- Seurat:::DiscretePalette(length(cluster_levels))
names(cluster_colors_seurat) <- cluster_levels


cluster_maps <- lapply(names(split_objects), function(sample_name) {
  
  obj <- split_objects[[sample_name]]
  
  celltypes <- sort(unique(obj$seurat_clusters))
  
  plots <- lapply(celltypes, function(ct) {
    
    obj$highlight <- ifelse(
      obj$seurat_clusters == ct,
      as.character(ct),
      "Other"
    )
    
    ImageDimPlot(
      obj,
      group.by = "highlight",
      cols = c(
        "Other" = "grey85",
        cluster_colors_seurat[as.character(ct)]
      ),
      border.size = NA,
      dark.background = FALSE,
      size = 1
    ) +
      ggtitle(paste(sample_name, "-", ct))
    
  })
  
  wrap_plots(plots, ncol = 3)
})

combined_plot <- wrap_plots(cluster_maps) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right"
  )

ggsave(
  file.path(
    result_dir,
    "9.0_ImageDimPlot_Fibro_cluster_highlighted_0.115.pdf"
  ),
  combined_plot,
  width = 45,
  height = 30, 
  limitsize = F
)


##### 9.1 Marker genes ######
Fibro_seurat <- PrepSCTFindMarkers(Fibro_seurat)
qsave(Fibro_seurat, file.path(result_dir, "9.1_Fibro_subset_integrated_clustering_res0.115_downstream_analysis.qs"))

Idents(Fibro_seurat) <- "seurat_clusters"
markers_seurat <- FindAllMarkers(object=Fibro_seurat,
                                 test.use="wilcox",
                                 layer="data",
                                 only.pos=FALSE,
                                 max.cells.per.ident=1000)

markers_sorted <- markers_seurat %>%
  arrange(cluster, p_val_adj, desc(avg_log2FC))
marker_list <- split(markers_sorted, markers_sorted$cluster)
write.xlsx(marker_list, file.path(result_dir, "9.1_Cluster_Markers_seurat_res0.115.xlsx"))

markers_top = markers_seurat %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::arrange(p_val_adj, avg_log2FC) %>% 
  dplyr::slice_head(n=3) %>% 
  dplyr::ungroup() %>% 
  as.data.frame()

p <- DotPlot(Fibro_seurat, features=markers_top$gene %>% unique()) +
  viridis::scale_color_viridis() + 
  ylab("Cluster") + xlab("") + 
  guides(size=guide_legend(order=1, title="Percent expressed"), color=guide_colorbar(title="Scaled Expression")) + 
  ggtitle("Top markers per cluster (expression scaled)") + RotatedAxis()
ggsave(file.path(result_dir, "9.1_DP_top_marker_genes_seurat_res0.115.pdf"), p , width = 15, heigh = 7)

#genes of interest snRNA-seq 
genes <- c("ANGPT1", "CFH", "FGF7", "SVEP1", "PI16", "SLC4A4", "KCNMA1")

p <- FeaturePlot(Fibro_seurat, features = genes)
ggsave(file.path(result_dir, "9.1_FP_Fibro_reclustered_genes_of_interest.pdf"), p, width = 10, height = 15)


#colored
genes <- c("ANGPT1", "CFH", "FGF7", "SVEP1", "PI16", "KCNMA1") #"SLC4A4" (not included)

expr <- FetchData(Fibro_seurat, vars = genes)

Fibro_seurat$marker <- apply(expr, 1, function(x){
  
  if(all(x == 0))
    return("None")
  
  genes[which.max(x)]
  
})

Fibro_seurat$marker <- factor(
  Fibro_seurat$marker,
  levels = c(genes, "None")
)

cols <- c(
  ANGPT1 = "#ee756d",
  CFH = "#e31a1c",
  FGF7 = "#66c69b",
  SVEP1 = "#a6cee3",
  PI16 = "#d49005",
  #SLC4A4 = "#87be4d",
  KCNMA1 = "#87be4d",
  None = "grey90"
)

DimPlot(
  Fibro_seurat,
  group.by = "marker",
  cols = cols,
  raster = FALSE
)

#make colored plots bigger
# Extract UMAP coordinates
umap_df <- Embeddings(
  Fibro_seurat,
  reduction = "umap"
) %>%
  as.data.frame() %>%
  rownames_to_column("cell")

# Add marker metadata
umap_df$marker <- Fibro_seurat$marker[
  match(
    umap_df$cell,
    colnames(Fibro_seurat)
  )
]

# Make sure factor levels are preserved
umap_df$marker <- factor(
  umap_df$marker,
  levels = c(genes, "None")
)


p <- ggplot() +
  geom_point(
    data = umap_df %>%
      filter(marker == "None"),
    aes(
      x = umap_1,
      y = umap_2
    ),
    color = "grey90",
    size = 0.3,
    alpha = 0.8
  ) +
  geom_point(
    data = umap_df %>%
      filter(marker != "None"),
    aes(
      x = umap_1,
      y = umap_2,
      color = marker
    ),
    size = 1.2,
    alpha = 0.9
  ) +
  scale_color_manual(
    values = cols[names(cols) != "None"],
    drop = FALSE
  ) +
  coord_fixed() +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  )

p

ggsave(file.path(result_dir, "9.1_DimPlot_Fibro_reclustered_genes_of_interest_highlighted.pdf"), p, width = 5, height = 5)


#ImageFeaturePlot
pdf(
  file.path(result_dir, "ß.1_ImageFeaturePlot_Fibro_all_genes.pdf"),
  width = 14,
  height = 12
)

for (gene in genes) {
  
  gene_maps <- lapply(names(split_objects), function(sample_name) {
    
    ImageFeaturePlot(
      split_objects[[sample_name]],
      features = gene,
      border.size = NA,
      dark.background = FALSE,
      size = 1.5,
      cols = c("grey95", "red")
    ) +
      ggtitle(sample_name) +
      theme_void() +
      theme(
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        legend.position = "right",
        plot.margin = margin(2, 2, 2, 2)
      )
  })
  
  gene_plot <- wrap_plots(gene_maps, ncol = 3) +
    plot_annotation(
      title = paste0("Fibroblast expression: ", gene),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
      )
    )
  
  print(gene_plot)
}

#### 9.2 DE analysis individually ####
DefaultAssay(Fibro_seurat) <- "SCT"

Fibro_seurat <- PrepSCTFindMarkers(
  Fibro_seurat,
  assay = "SCT",
  verbose = TRUE
)

DefaultAssay(Fibro_seurat) <- "SCT"

Fibro_seurat$annotated <- Fibro_seurat$seurat_clusters #change later
Idents(Fibro_seurat) <- "annotated"
de_nested_list <- list()

comparisons <- list(
  c("MS", "Ctrl")
)

for (comp in comparisons) {
  ident_1 <- comp[1]
  ident_2 <- comp[2]
  comp_name <- paste0(ident_1, "_vs_", ident_2)
  
  message("Running comparison: ", comp_name)
  de_nested_list[[comp_name]] <- list()
  
  for (i in levels(Fibro_seurat$annotated)) {
    message("  Processing cluster: ", i)
    
    cells_in_cluster <- WhichCells(Fibro_seurat, idents = i)
    cluster_meta <- Fibro_seurat@meta.data[cells_in_cluster, ]
    
    group_counts <- table(cluster_meta$condition)
    
    cluster_name <- paste0("Cluster_", i)
    
    if (all(c(ident_1, ident_2) %in% names(group_counts)) &&
        all(group_counts[c(ident_1, ident_2)] >= 3)) {
      
      de <- tryCatch({
        FindMarkers(
          Fibro_seurat,
          ident.1 = ident_1,
          ident.2 = ident_2,
          subset.ident = i,
          group.by = "condition",
          recorrect_umi = FALSE
        )
      }, error = function(e) {
        message("    ERROR in cluster ", i, ": ", conditionMessage(e))
        return(NULL)
      })
      
      if (!is.null(de) && nrow(de) > 0) {
        de$gene <- rownames(de)
        de <- de[order(de$p_val_adj, -de$avg_log2FC), ]
        de_nested_list[[comp_name]][[cluster_name]] <- de
      } else {
        message("    No DE genes found for cluster ", i)
        de_nested_list[[comp_name]][[cluster_name]] <- data.frame(
          Note = "No DE genes found despite sufficient cells"
        )
      }
      
    } else {
      message("    Skipping cluster ", i, " — too few cells in one or both groups")
      de_nested_list[[comp_name]][[cluster_name]] <- data.frame(
        Note = paste("Too few cells in", ident_1, "or", ident_2)
      )
    }
  }
}

saveRDS(de_nested_list, "9.2_de_nested_list_SCT_condition.rds")

# Save each comparison as its own Excel file
for (comp_name in names(de_nested_list)) {
  file_name <- file.path(
    result_dir,
    paste0("9.2_DE_SCT_", comp_name, ".xlsx")
  )
  write.xlsx(de_nested_list[[comp_name]], file = file_name)
  message("Saved: ", file_name)
}

##### 9.3 DE duraFibro merged ####
Idents(Fibro_seurat) <- "condition"

DefaultAssay(Fibro_seurat) <- "SCT"
DE_duraFibro_overall <- FindMarkers(
  Fibro_seurat,
  ident.1 = "MS",
  ident.2 = "Ctrl",
  group.by = "condition",
  recorrect_umi = FALSE
)

DE_duraFibro_overall <- DE_duraFibro_overall %>%
  rownames_to_column(var = "gene") %>%
  arrange(p_val_adj, desc(avg_log2FC))

write_xlsx(DE_duraFibro_overall, file.path(result_dir, "9.3_DE_duraFibro_overall.xlsx"))


###### 9.3.1 GO duraFibro DE ######
set.seed(42)

read_all_sheets <- function(path) {
  sheet_names <- excel_sheets(path)
  sheets <- lapply(sheet_names, function(s) read_excel(path, sheet = s))
  names(sheets) <- sheet_names
  return(sheets)
}
DE_MS_vs_Ctrl <- read_all_sheets(file.path(result_dir, "9.3_DE_duraFibro_overall.xlsx"))

# Function for enrichment and plotting
run_enrichment_and_plotting <- function(de_df, gene_list, direction, cluster, comparison) {
  if (length(gene_list) < 5) {
    message("Not enough genes for ", direction, "regulated ", cluster, " in ", comparison)
    return(NULL)
  }
  ego <- enrichGO(
    gene          = gene_list,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "ALL",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  if (is.null(ego) || nrow(ego) == 0) {
    message("No enrichment results for ", direction, "regulated ", cluster, " in ", comparison)
    return(NULL)
  }
  # Save GO table
  file_base <- paste0(file.path(result_dir, "9.3.1_GO_DE_duraFibro_"), direction, "", cluster, "", comparison)
  write.xlsx(as.data.frame(ego), file = paste0(file_base, ".xlsx"))
  # Plots
  pdf(paste0(file_base, ".pdf"), width = 8, height = 6)
  tryCatch({
    print(dotplot(ego) + ggtitle(paste("Dotplot -", direction, "-", cluster, "-", comparison)))
    print("dotplot")
    print(barplot(ego) + ggtitle(paste("Dotplot -", direction, "-", cluster, "-", comparison)))
    print("barplot")
    sim <- pairwise_termsim(ego)
    print(emapplot(sim) + ggtitle(paste("Enrichment Map -", direction, "-", cluster, "-", comparison)))
    print("emapplot")
    print(emapplot(sim, group = T) + ggtitle(paste("Cluster Map -", direction, "-", cluster, "-", comparison)))
    print("emapplot_cluster")
    de_df <- as.data.frame(de_df)
    rownames(de_df) <- de_df$gene
    fc_vector <- de_df[gene_list, "avg_log2FC"]
    names(fc_vector) <- gene_list
    print(cnetplot(ego, categorySize = "pvalue", foldChange = fc_vector) +
            ggtitle(paste("Cnetplot -", direction, "-", cluster, "-", comparison)))
    print("cnetplot")
  }, error = function(e) {
    message("Plotting failed for ", direction, "regulated ", cluster, " in ", comparison, ": ", e$message)
  })
  dev.off()
}

#Main loop
de_df <- DE_MS_vs_Ctrl[["Sheet1"]]

# Upregulated genes
gene_list_up <- de_df %>%
  filter(avg_log2FC > 0.5 & p_val_adj < 0.001) %>%
  pull(gene)

run_enrichment_and_plotting(de_df, gene_list_up, "Up", "All_duraFibros", "MS_vs_Ctrl")

# Downregulated genes
gene_list_down <- de_df %>%
  filter(avg_log2FC < -0.5 & p_val_adj < 0.001) %>%
  pull(gene)

run_enrichment_and_plotting(de_df, gene_list_down, "Down", "All_duraFibros", "MS_vs_Ctrl")


##### 9.3.2 volcanoPlot duraFibro DE top DE genes label #####
de_file <- file.path(result_dir, "9.3_DE_duraFibro_overall.xlsx") 

all_clusters <- getSheetNames(de_file)

#Loop over clusters 
for (cluster in all_clusters) {
  
  message("Processing cluster: ", cluster)
  
  # Read cluster sheet
  de_data <- read.xlsx(de_file, sheet = cluster)
  volcano_genes <- de_data
  
  #remove MT, RPL, RPS genes and Lnc 
  volcano_genes <- volcano_genes %>% filter(!grepl("^MT|^RPL|^RPS|^ZNF|orf|^LINC|^ENSG|\\.", gene))
  #Define sig. up/down-reg genes
  volcano_genes$DE <- "NO"
  volcano_genes$DE[volcano_genes$avg_log2FC > 0.5 & volcano_genes$p_val_adj < 0.001] <- "Up" 
  volcano_genes$DE[volcano_genes$avg_log2FC < -0.5 & volcano_genes$p_val_adj < 0.001] <- "Down" 
  #sign. up & down
  volcano_genes$sig <- ifelse(volcano_genes$DE == "NO", "", volcano_genes$gene)
  sig <- volcano_genes$sig
  volcano_genes$sig <- ifelse(volcano_genes$gene %in% sig, volcano_genes$gene, NA)
  volcano_genes$DE <- factor(volcano_genes$DE, levels = c("Down", "NO", "Up"))
  volcano_genes$logP <- -log10(volcano_genes$p_val_adj)
  volcano_genes$logP <- pmin(volcano_genes$logP, 50)
  
  # Count up- and down-regulated genes
  n_up <- sum(volcano_genes$DE == "Up", na.rm = TRUE)
  n_down <- sum(volcano_genes$DE == "Down", na.rm = TRUE)
  
  # Filter only significant DE genes
  sig_genes <- volcano_genes %>% 
    filter(DE != "NO") 
  # Rank first by adjusted p-value, then by absolute log2FC
  top_genes <- sig_genes %>%
    arrange(p_val_adj, desc(abs(avg_log2FC))) %>%  # first p-value, then log2FC
    slice_head(n = 400) %>%                         # top 500 genes
    pull(gene)
  
  # Update sig column to label only top 400 genes
  volcano_genes$interest <- ifelse(volcano_genes$gene %in% top_genes, volcano_genes$gene, NA)
  
  p <- ggplot(volcano_genes, aes(x = avg_log2FC, y = logP, color = DE, label = sig)) +
    geom_vline(xintercept = c(-0.5, 0.5), color = "gray40", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.001), color = "gray40", linetype = "dashed") + 
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = c("Down" = "black", "NO" = "grey70", "Up" = "#3d70b5"),
                       labels = c("Down", "ns", "Up")) +
    geom_text_repel(max.overlaps = 50, size = 3) +
    theme_classic(base_size = 14) +
    labs(
      x = "log Fold Change",
      y = expression(-log[10]("adjusted p-value")),
      color = "Regulation"
    ) +
    annotate("text", x = max(volcano_genes$avg_log2FC), y = max(volcano_genes$logP),
             label = paste0("Up: ", n_up), hjust = 1, vjust = -0.5, color = "#3d70b5", size = 5) +
    annotate("text", x = min(volcano_genes$avg_log2FC), y = max(volcano_genes$logP),
             label = paste0("Down: ", n_down), hjust = 0, vjust = -0.5, color = "black", size = 5)
  
  # Save per cluster
  out_file <- paste0(file.path(result_dir, "9.3.2_VolcanoPlot_duraFibro_"), cluster, "_DE.pdf")
  ggsave(out_file, plot = p, width = 10, height = 7)
}

rm(de_data, p, sig_genes, volcano_genes, all_clusters, cluster, de_file, n_down, n_up, out_file, sig, top_genes)


#==== 10.0 Score ====

#### 10.1 fibrosis snRNA-seq score ####
fibrosis <- readr::read_csv(
  file.path(result_dir, "fibrosis.txt"), #genes used in snRNA-seq
  col_names = FALSE,
  show_col_types = FALSE
) %>%
  dplyr::pull(X1) %>%
  unique() %>%
  na.omit()

#check which genes are present
fibrosis_present <- intersect(fibrosis, rownames(sc_combined))
fibrosis_missing <- setdiff(fibrosis, rownames(sc_combined))
fibrosis_present
fibrosis_missing
length(fibrosis_present)

##### 10.1.1 all cluster #####
DefaultAssay(sc_combined) <- "SCT"
fibrosis_score <- AddModuleScore(sc_combined, features = list(fibrosis_present), name = "fibrosisScore", nbin = 5)
FeaturePlot(fibrosis_score, features ="fibrosisScore1", split.by = "condition")
VlnPlot(fibrosis_score, features ="fibrosisScore1", split.by = "condition",  cols = c("grey", "darkblue"))
Idents(fibrosis_score) <- "condition"
DotPlot(fibrosis_score, features ="fibrosisScore1")

df <- FetchData(fibrosis_score, vars = c("fibrosisScore1", "condition", "seurat_clusters", "sample"))
df_avg <- df %>%
  group_by(seurat_clusters, condition, sample) %>%
  summarise(fibrosisScore1 = mean(fibrosisScore1, na.rm = TRUE), .groups = "drop")
df_collapsed <- df_avg %>%
  group_by(sample, condition) %>%
  summarise(fibrosisScore1 = mean(fibrosisScore1), .groups = "drop")

pdf(file.path(result_dir, "10.1.1_Xenium_5k_dura_fibrosis_score_MS_vs_Ctrl_all_cluster.pdf"), width = 4, height = 6)
ggplot(df_collapsed, aes(x = condition, y = fibrosisScore1, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.6, size = 2, width = 0.2) +
  stat_compare_means(label = "p.signif", method = "wilcox.test") +
  scale_fill_manual(values = c("Ctrl" = "darkgrey", "MS" = "#3d70b5")) +
  theme_classic() +
  labs(title = "fibrosis score all cluster", x = "Condition", y = "average fibrosis score per sample")
dev.off()

##### 10.1.2 Fibro cluster #####
DefaultAssay(Fibro_seurat) <- "SCT"
Idents(Fibro_seurat) <- "seurat_clusters"

fibrosis_score <- AddModuleScore(Fibro_seurat, features = list(fibrosis_present), name = "fibrosisScore", nbin = 5)
FeaturePlot(fibrosis_score, features ="fibrosisScore1", split.by = "condition")
VlnPlot(fibrosis_score, features ="fibrosisScore1", split.by = "condition",  cols = c("grey", "darkblue"))
Idents(fibrosis_score) <- "condition"
DotPlot(fibrosis_score, features ="fibrosisScore1")

df <- FetchData(fibrosis_score, vars = c("fibrosisScore1", "condition", "seurat_clusters", "sample"))
df_avg <- df %>%
  group_by(seurat_clusters, condition, sample) %>%
  summarise(fibrosisScore1 = mean(fibrosisScore1, na.rm = TRUE), .groups = "drop")
df_collapsed <- df_avg %>%
  group_by(sample, condition) %>%
  summarise(fibrosisScore1 = mean(fibrosisScore1), .groups = "drop")

pdf(file.path(result_dir, "10.1.2_Xenium_5k_dura_fibrosis_score_MS_vs_Ctrl_duraFibro.pdf"), width = 4, height = 6)
ggplot(df_collapsed, aes(x = condition, y = fibrosisScore1, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.6, size = 2, width = 0.2) +
  stat_compare_means(label = "p.signif", method = "wilcox.test") +
  scale_fill_manual(values = c("Ctrl" = "darkgrey", "MS" = "#3d70b5")) +
  theme_classic() +
  labs(title = "fibrosis score all cluster", x = "Condition", y = "average fibrosis score per sample")
dev.off()


#### 10.2 5k specific scores ####
Idents(sc_combined) <- "snRNA_major_group_mix"
Tc_subset <- subset(sc_combined, idents = "Tc")
Bc_subset <- subset(sc_combined, idents = "Bc")
Myeloid_subset <- subset(sc_combined, idents = "Myeloid")


gene_sets <- list(
  ECM = c("POSTN","TNC","COL4A1","COL4A2","COL5A1","COL5A2","LOX","LOXL1","LOXL2","SERPINE1","THBS1","THBS2","HSPG2","ELN","COMP"),
  activated_fibroblasts = c("POSTN","FAP","PDGFRA","PDGFRB","CD248","LRRC15","CTHRC1","PRRX1","TWIST1","TWIST2"),
  TGFB_fibrosis = c("TGFB1","TGFB2","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SERPINE1","THBS1","THBS2"),
  fibroblast_identity = c("DPT","PI16","PDGFRA","FBLN2","FBLN5","MFAP5","SFRP1","SFRP4"),
  fibrosis_score = c("POSTN","FAP","CTHRC1","COL5A1","COL5A2","LOX","SERPINE1","THBS1","PDGFRA","TGFB1"),
  leukocyte_activation_general = c("CD44","CD69","CD38","ICOS","TNFRSF9","CD40","CD83","CD86"),
  activated_Tcells = c("CD69","ICOS","TNFRSF9","IFNG","TNF","LTA","PRF1","GZMB","GZMK","GZMA"),
  activated_APCs = c("CD80","CD86","CD83","CIITA","TAP1","TAP2","ICAM1"),
  activated_myeloid = c("IL1B","TNF","CXCL9","CXCL10","CXCL11","PTGS2","CCL19","CCL20"),
  leukocyte_activation_score = c("CD69","CD38","ICOS","TNF","IFNG","CD83","CD86","IL1B","CXCL10","ICAM1"),
  typeI_IFN = c("IFI35","IFI44L","IFIH1","IFIT1","IFIT2","IFIT3","MX1","OAS1","OAS3","RSAD2","PLSCR1","SIGLEC1","HERC6","SP110","SP140"),
  IFN_signaling = c("STAT1","STAT2","IRF1","IRF3","IRF5","IRF7","IRF8","JAK1","JAK2","TYK2","IFNAR1","IFNAR2","IFNGR1","IFNGR2"),
  IFNg_response = c("CXCL9","CXCL10","CXCL11","GBP1","GBP5","IRF1","STAT1","TAP1","TAP2","CIITA"),
  typeI_IFN_score = c("IFI44L","IFIT1","IFIT2","IFIT3","MX1","OAS1","OAS3","RSAD2","STAT1","IRF7"),
  IFNg_score = c("CXCL9","CXCL10","CXCL11","GBP1","GBP5","STAT1","IRF1","TAP1","TAP2","CIITA")
)


make_score_plot <- function(obj, genes, score_name, title_text) {
  DefaultAssay(obj) <- "SCT"
  
  genes_use <- intersect(genes, rownames(obj))
  
  obj <- AddModuleScore(
    obj,
    features = list(genes_use),
    name = score_name,
    nbin = 5
  )
  
  score_col <- paste0(score_name, "1")
  
  df <- FetchData(
    obj,
    vars = c(score_col, "condition", "seurat_clusters", "sample")
  )
  
  colnames(df)[colnames(df) == score_col] <- "score"
  
  df_avg <- df %>%
    group_by(seurat_clusters, condition, sample) %>%
    summarise(score = mean(score, na.rm = TRUE), .groups = "drop")
  
  df_collapsed <- df_avg %>%
    group_by(sample, condition) %>%
    summarise(score = mean(score, na.rm = TRUE), .groups = "drop")
  
  ggplot(df_collapsed, aes(x = condition, y = score, fill = condition)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(alpha = 0.6, size = 2, width = 0.2) +
    stat_compare_means(label = "p.signif", method = "wilcox.test") +
    scale_fill_manual(values = c("Ctrl" = "darkgrey", "MS" = "#3d70b5")) +
    theme_classic() +
    labs(
      title = title_text,
      x = "Condition",
      y = paste0("average ", title_text, " per sample")
    )
}

plot_all_scores_to_pdf <- function(obj, object_name, output_file) {
  pdf(output_file, width = 4, height = 6)
  
  for (sig_name in names(gene_sets)) {
    p <- make_score_plot(
      obj = obj,
      genes = gene_sets[[sig_name]],
      score_name = paste0(sig_name, "_Score"),
      title_text = paste(sig_name, object_name)
    )
    
    print(p)
  }
  
  dev.off()
}

plot_all_scores_to_pdf(
  obj = sc_combined,
  object_name = "all clusters",
  output_file = file.path(result_dir, "10.2_Xenium_5k_all_cluster_all_module_scores.pdf")
)

plot_all_scores_to_pdf(
  obj = Fibro_seurat,
  object_name = "duraFibro",
  output_file = file.path(result_dir, "10.2_Xenium_5k_duraFibro_all_module_scores.pdf")
)


plot_all_scores_to_pdf(
  obj = Tc_subset,
  object_name = "Tc",
  output_file = file.path(result_dir, "10.2_Xenium_5k_Tc_all_module_scores.pdf")
)

plot_all_scores_to_pdf(
  obj = Tc_subset,
  object_name = "Bc",
  output_file = file.path(result_dir, "10.2_Xenium_5k_Bc_all_module_scores.pdf")
)

plot_all_scores_to_pdf(
  obj = Myeloid_subset,
  object_name = "Myeloid",
  output_file = file.path(result_dir, "10.2_Xenium_5k_Myeloid_all_module_scores.pdf")
)


#==== 11.0 Cell interactions ======
##### 11.1 NichetNet snRNA_major_group_mix remove not well defined cluster #####
#convert Seurat to single cell experiment and define meta data
sce <- Seurat::as.SingleCellExperiment(sc_combined, assay = "SCT")
sce <- alias_to_symbol_SCE(sce, "human") |> makenames_SCE()
#remove not well-defined clusters
exclude_celltypes <- c("mix", "unknown")
sce <- sce[, !colData(sce)[, celltype_id] %in% exclude_celltypes]

sample_id <- "sample"
group_id <- "condition"
celltype_id <- "snRNA_major_group_mix" 
covariates <- NA
batches <- NA

#define organism
organism = "human"

if(organism == "human"){
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  lr_network = lr_network |> dplyr::rename(ligand = from, receptor = to) |> distinct(ligand, receptor) |> mutate(ligand = make.names(ligand), receptor = make.names(receptor))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) |> make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) |> make.names()
} else if(organism == "mouse"){
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  lr_network = lr_network |> dplyr::rename(ligand = from, receptor = to) |> distinct(ligand, receptor) |> mutate(ligand = make.names(ligand), receptor = make.names(receptor))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) |> make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) |> make.names()
}

#change names because _ and dots are not allowed
colData(sce)[, celltype_id] <- factor(gsub("-", "_", colData(sce)[, celltype_id]))
levels(colData(sce)[, celltype_id])

#Sender and receiver cell types need to be defined. Here all cell types in the dataset are used
senders_oi <- SummarizedExperiment::colData(sce)[,celltype_id] |> unique()
receivers_oi <- SummarizedExperiment::colData(sce)[,celltype_id] |> unique()

sce <- sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% c(senders_oi, receivers_oi)]

min_cells <- 10

abundance_info = get_abundance_info(
  sce = sce, 
  sample_id = sample_id, 
  group_id = group_id, 
  celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, 
  receivers_oi = receivers_oi, 
  batches = batches)

## check for condition-specific cell types
abundance_df_summarized = abundance_info$abundance_data %>% mutate(keep = as.logical(keep)) %>% group_by(group_id, celltype_id) %>% summarise(samples_present = sum((keep)))
celltypes_absent_one_condition = abundance_df_summarized %>% filter(samples_present == 0) %>% pull(celltype_id) %>% unique() # find truly condition-specific cell types by searching for cell types truely absent in at least one condition
celltypes_present_one_condition = abundance_df_summarized %>% filter(samples_present >= 2) %>% pull(celltype_id) %>% unique() # require presence in at least 2 samples of one group so it is really present in at least one condition
condition_specific_celltypes = intersect(celltypes_absent_one_condition, celltypes_present_one_condition)

condition_specific_celltypes_senders = condition_specific_celltypes %>% intersect(senders_oi)
condition_specific_celltypes_receivers = condition_specific_celltypes %>% intersect(receivers_oi)

fraction_cutoff <- 0.05
min_sample_prop <- 0.5
frq_list = get_frac_exprs(sce = sce, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, batches = batches, min_cells = min_cells, fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)

abundance_expression_info <- process_abundance_expression_info(
  sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, min_cells = min_cells,
  senders_oi = union(senders_oi, condition_specific_celltypes_senders),
  receivers_oi = union(receivers_oi, condition_specific_celltypes_receivers),
  lr_network = lr_network, batches = batches, frq_list = frq_list, abundance_info = abundance_info)

#If you compare against two groups, you should divide by 2, if you compare against three groups, you should divide by 3 etcetera
contrasts_oi <- c("'M-(S+A)/2','S-(M+A)/2','A-(S+M)/2'")
#if compare against two groups
contrasts_oi <- c("'MS-Ctrl','Ctrl-MS'") 
contrast_tbl <- tibble(contrast = c("MS-Ctrl","Ctrl-MS"), group = c("MS","Ctrl"))

assay_oi_pb <- "counts"
fun_oi_pb <- "sum"
de_method_oi <- "edgeR" 
findMarkers <- FALSE

DE_info = get_DE_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, batches = batches, covariates = covariates, contrasts_oi = contrasts_oi, min_cells = min_cells, 
                      assay_oi_pb = assay_oi_pb,
                      fun_oi_pb = fun_oi_pb,
                      de_method_oi = de_method_oi,
                      findMarkers = findMarkers, 
                      expressed_df = frq_list$expressed_df)

#check p-value distribution#should be quite straight horizontal line
DE_info$hist_pvals

#if it looks ok, run "empirical_pval == FALSE 
empirical_pval = FALSE
#if it looks suboptimal run empirical_pval == TRUE and check quality
empirical_pval = TRUE

if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
}
DE_info_emp$hist_pvals_emp

if(empirical_pval == FALSE){
  celltype_de = DE_info$celltype_de$de_output_tidy
} else {
  celltype_de = DE_info_emp$de_output_tidy_emp %>% dplyr::select(-p_val, -p_adj) %>% dplyr::rename(p_val = p_emp, p_adj = p_adj_emp)
}

sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)

logFC_threshold = 0.50
p_val_threshold = 0.05
fraction_cutoff = 0.05


# p_val_adj = TRUE #In case of more samples per group + sufficient high number of DE genes per group-celltype (> 20-50), use p_val_adj = TRUE
p_val_adj = T

top_n_target = 400

verbose = TRUE
cores_system = 8
n.cores = min(cores_system, sender_receiver_de$receiver %>% unique() %>% length()) # use one core per receiver cell type

ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(get_ligand_activities_targets_DEgenes(
  receiver_de = celltype_de,
  receivers_oi = receivers_oi,
  ligand_target_matrix = ligand_target_matrix,
  logFC_threshold = logFC_threshold,
  p_val_threshold = p_val_threshold,
  p_val_adj = T,
  top_n_target = top_n_target,
  verbose = T, 
  n.cores = n.cores
)))

prioritizing_weights_DE = c("de_ligand" = 1,
                            "de_receptor" = 1)
prioritizing_weights_activity = c("activity_scaled" = 2)

prioritizing_weights_expression_specificity = c("exprs_ligand" = 2,
                                                "exprs_receptor" = 2)

prioritizing_weights_expression_sufficiency = c("frac_exprs_ligand_receptor" = 1)

prioritizing_weights_relative_abundance = c( "abund_sender" = 0,
                                             "abund_receiver" = 0)

prioritizing_weights = c(prioritizing_weights_DE, 
                         prioritizing_weights_activity, 
                         prioritizing_weights_expression_specificity,
                         prioritizing_weights_expression_sufficiency, 
                         prioritizing_weights_relative_abundance)

sender_receiver_tbl = sender_receiver_de |> dplyr::distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) |> tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] |> tibble::as_tibble() |> dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] |> tibble::as_tibble() |> dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group")
}

scenario <- "regular" # all prioritization criteria will be weighted equally
ligand_activity_down <- FALSE # use only upregulatory ligand activities to prioritize

prioritization_tables = suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = scenario, 
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = ligand_activity_down 
))

lr_target_prior_cor = lr_target_prior_cor_inference(prioritization_tables$group_prioritization_tbl$receiver |> unique(), abundance_expression_info, celltype_de, grouping_tbl, prioritization_tables, ligand_target_matrix, logFC_threshold = logFC_threshold, p_val_threshold = p_val_threshold, p_val_adj = p_val_adj)

#save previous calculations as one large list
multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
) 
multinichenet_output = make_lite_output(multinichenet_output)

saveRDS(multinichenet_output, file.path(result_dir, "11.1_multinichenet.rds"))

#### Visualization Circos plot 
#define colors based on our cluster_color 
cluster_colors <- c( "ABC" = "#9584bd",
                     "arach_bordFibro" = "#87be4d",
                     "duraFibro" = "#f2cf22", 
                     "vascStroma" = "#00b8df", 
                     "Endo" = "#27a9b5", 
                     "Myeloid" = "#e04c73",
                     "Bc" = "#bf72ad",
                     "Tc" = "#dd4e10", 
                     "mix" = "grey",
                     "unknown" = "grey")

#create plot for the top100 predictions
prioritized_tbl_oi_all = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 100, rank_per_group = FALSE)

prioritized_tbl_oi = multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender <- cluster_colors[senders_receivers]
colors_receiver <- cluster_colors[senders_receivers]

pdf(file.path(result_dir, "11.1_top100_circos_plot_snRNA_major_group_mix_removed_cluster.pdf"), width = 10, height = 10)
circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)
dev.off()

#split by receiver
# Open multi-page PDF
pdf(file.path(result_dir, "11.1_circos_split_by_receiver_SnRNA_major_group_mix_removed_cluster.pdf"), width = 10, height = 10)

unique_receivers <- sort(unique(prioritized_tbl_oi$receiver))

for(receiver in unique_receivers) {
  
  # Subset top interactions for this receiver
  subset_tbl <- prioritized_tbl_oi %>% filter(receiver == !!receiver)
  if(nrow(subset_tbl) == 0) next
  
  # Create circos object
  circos_list <- make_circos_group_comparison(subset_tbl, colors_sender, colors_receiver)
  
  circos.clear()
  circos.par(gap.degree = 5)
  
  # Initialize sectors
  cell_types <- circos_list$cell_types
  xlim_list <- lapply(cell_types, function(x) c(0,1))
  names(xlim_list) <- cell_types
  for(sector in cell_types){
    circos.initialize(factors = sector, xlim = xlim_list[[sector]])
  }
  
  # Draw arcs and labels
  cell_colors <- circos_list$colors
  for(sector in cell_types){
    circos.trackPlotRegion(factors = sector, y = c(0,1),
                           bg.col = cell_colors[sector], bg.border = NA,
                           panel.fun = function(x, y){
                             circos.text(CELL_META$xcenter, CELL_META$ycenter,
                                         CELL_META$sector.index,
                                         facing = "clockwise", niceFacing = TRUE, cex = 0.8)
                           })
  }
  
  # Draw links
  if(!is.null(circos_list$links) && nrow(circos_list$links) > 0){
    for(i in seq_len(nrow(circos_list$links))){
      link <- circos_list$links[i, ]
      circos.link(sector.index1 = link$sender, point1 = 0.5,
                  sector.index2 = link$receiver, point2 = 0.5,
                  col = rgb(0,0,0,0.3), border = NA)
    }
  }
  
  # Each iteration creates a **new page** in the PDF
}

# Close PDF device
dev.off()

#### 11.2 NicheNet MS duraFibro #####
#convert Seurat to single cell experiment and define meta data
sc_combined$Fibro <- dplyr::case_when(
  sc_combined$snRNA_major_group_mix %in% c("duraFibro", "arach_bordFibro") ~ "duraFibro",
  TRUE ~ as.character(sc_combined$snRNA_major_group_mix)
)

sce <- Seurat::as.SingleCellExperiment(sc_combined, assay = "SCT")
sce <- alias_to_symbol_SCE(sce, "human") |> makenames_SCE()

sample_id <- "sample"
group_id <- "condition"
celltype_id <- "Fibro"
covariates <- NA
batches <- NA

exclude_celltypes <- c("mix", "unknown")
sce <- sce[, !colData(sce)[, celltype_id] %in% exclude_celltypes]


#define organism
organism = "human"

if(organism == "human"){
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  lr_network = lr_network |> dplyr::rename(ligand = from, receptor = to) |> distinct(ligand, receptor) |> mutate(ligand = make.names(ligand), receptor = make.names(receptor))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) |> make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) |> make.names()
} else if(organism == "mouse"){
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  lr_network = lr_network |> dplyr::rename(ligand = from, receptor = to) |> distinct(ligand, receptor) |> mutate(ligand = make.names(ligand), receptor = make.names(receptor))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  colnames(ligand_target_matrix) = colnames(ligand_target_matrix) |> make.names()
  rownames(ligand_target_matrix) = rownames(ligand_target_matrix) |> make.names()
}

#change names because _ and dots are not allowed
colData(sce)[, celltype_id] <- factor(gsub("-", "_", colData(sce)[, celltype_id]))
levels(colData(sce)[, celltype_id])

#Sender and receiver cell types need to be defined. Here all cell types in the dataset are used
senders_oi <- SummarizedExperiment::colData(sce)[,celltype_id] |> unique()
receivers_oi <- SummarizedExperiment::colData(sce)[,celltype_id] |> unique()

sce <- sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% c(senders_oi, receivers_oi)]

min_cells <- 10

abundance_info = get_abundance_info(
  sce = sce, 
  sample_id = sample_id, 
  group_id = group_id, 
  celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, 
  receivers_oi = receivers_oi, 
  batches = batches)

## check for condition-specific cell types
abundance_df_summarized = abundance_info$abundance_data %>% mutate(keep = as.logical(keep)) %>% group_by(group_id, celltype_id) %>% summarise(samples_present = sum((keep)))
celltypes_absent_one_condition = abundance_df_summarized %>% filter(samples_present == 0) %>% pull(celltype_id) %>% unique() # find truly condition-specific cell types by searching for cell types truely absent in at least one condition
celltypes_present_one_condition = abundance_df_summarized %>% filter(samples_present >= 2) %>% pull(celltype_id) %>% unique() # require presence in at least 2 samples of one group so it is really present in at least one condition
condition_specific_celltypes = intersect(celltypes_absent_one_condition, celltypes_present_one_condition)

condition_specific_celltypes_senders = condition_specific_celltypes %>% intersect(senders_oi)
condition_specific_celltypes_receivers = condition_specific_celltypes %>% intersect(receivers_oi)

fraction_cutoff <- 0.05
min_sample_prop <- 0.5
frq_list = get_frac_exprs(sce = sce, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, batches = batches, min_cells = min_cells, fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)

abundance_expression_info <- process_abundance_expression_info(
  sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, min_cells = min_cells,
  senders_oi = union(senders_oi, condition_specific_celltypes_senders),
  receivers_oi = union(receivers_oi, condition_specific_celltypes_receivers),
  lr_network = lr_network, batches = batches, frq_list = frq_list, abundance_info = abundance_info)

table(colData(sce)[,celltype_id], colData(sce)[,sample_id])

#If you compare against two groups, you should divide by 2, if you compare against three groups, you should divide by 3 etcetera
contrasts_oi <- c("'M-(S+A)/2','S-(M+A)/2','A-(S+M)/2'")
#if compare against two groups
contrasts_oi <- c("'MS-Ctrl','Ctrl-MS'") 
contrast_tbl <- tibble(contrast = c("MS-Ctrl","Ctrl-MS"), group = c("MS","Ctrl"))

assay_oi_pb <- "counts"
fun_oi_pb <- "sum"
de_method_oi <- "edgeR" 
findMarkers <- FALSE

DE_info = get_DE_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, batches = batches, covariates = covariates, contrasts_oi = contrasts_oi, min_cells = min_cells, 
                      assay_oi_pb = assay_oi_pb,
                      fun_oi_pb = fun_oi_pb,
                      de_method_oi = de_method_oi,
                      findMarkers = findMarkers, 
                      expressed_df = frq_list$expressed_df)

#check p-value distribution#should be quite straight horizontal line
DE_info$hist_pvals

#if it looks ok, run "empirical_pval == FALSE 
empirical_pval = FALSE
#if it looks suboptimal run empirical_pval == TRUE and check quality
empirical_pval = TRUE

if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
}
DE_info_emp$hist_pvals_emp

if(empirical_pval == FALSE){
  celltype_de = DE_info$celltype_de$de_output_tidy
} else {
  celltype_de = DE_info_emp$de_output_tidy_emp %>% dplyr::select(-p_val, -p_adj) %>% dplyr::rename(p_val = p_emp, p_adj = p_adj_emp)
}

sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)

logFC_threshold = 0.50
p_val_threshold = 0.05
fraction_cutoff = 0.05


# p_val_adj = TRUE #In case of more samples per group + sufficient high number of DE genes per group-celltype (> 20-50), use p_val_adj = TRUE
p_val_adj = T

top_n_target = 400

verbose = TRUE
cores_system = 8
n.cores = min(cores_system, sender_receiver_de$receiver %>% unique() %>% length()) # use one core per receiver cell type

ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(get_ligand_activities_targets_DEgenes(
  receiver_de = celltype_de,
  receivers_oi = receivers_oi,
  ligand_target_matrix = ligand_target_matrix,
  logFC_threshold = logFC_threshold,
  p_val_threshold = p_val_threshold,
  p_val_adj = T,
  top_n_target = top_n_target,
  verbose = T, 
  n.cores = n.cores
)))

prioritizing_weights_DE = c("de_ligand" = 1,
                            "de_receptor" = 1)
prioritizing_weights_activity = c("activity_scaled" = 2)

prioritizing_weights_expression_specificity = c("exprs_ligand" = 2,
                                                "exprs_receptor" = 2)

prioritizing_weights_expression_sufficiency = c("frac_exprs_ligand_receptor" = 1)

prioritizing_weights_relative_abundance = c( "abund_sender" = 0,
                                             "abund_receiver" = 0)

prioritizing_weights = c(prioritizing_weights_DE, 
                         prioritizing_weights_activity, 
                         prioritizing_weights_expression_specificity,
                         prioritizing_weights_expression_sufficiency, 
                         prioritizing_weights_relative_abundance)

sender_receiver_tbl = sender_receiver_de |> dplyr::distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) |> tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] |> tibble::as_tibble() |> dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] |> tibble::as_tibble() |> dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group")
}

scenario <- "regular" # all prioritization criteria will be weighted equally
ligand_activity_down <- FALSE # use only upregulatory ligand activities to prioritize

prioritization_tables = suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = scenario, 
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = ligand_activity_down 
))

lr_target_prior_cor = lr_target_prior_cor_inference(prioritization_tables$group_prioritization_tbl$receiver |> unique(), abundance_expression_info, celltype_de, grouping_tbl, prioritization_tables, ligand_target_matrix, logFC_threshold = logFC_threshold, p_val_threshold = p_val_threshold, p_val_adj = p_val_adj)

#extract MS-specific strucFibro interactions
interactions_MS <- prioritization_tables$group_prioritization_tbl %>%
  filter(group == "MS")

duraFibro_MS_interactions <- interactions_MS %>%
  filter(sender == "duraFibro" | receiver == "duraFibro") %>%
  arrange(desc(prioritization_score))

# Inspect partners
unique_senders_MS <- setdiff(unique(duraFibro_MS_interactions$sender), "duraFibro")
unique_receivers_MS <- setdiff(unique(duraFibro_MS_interactions$receiver), "duraFibro")

cat("Senders to duraFibro in MS:\n")
print(unique_senders_MS)

cat("Receivers from duraFibro in MS:\n")
print(unique_receivers_MS)

# Save table
write.csv(
  duraFibro_MS_interactions,
  file = file.path(result_dir, "11.2_MS_duraFibro_interactions.csv"),
  row.names = FALSE
)

#save previous calculations as one large list
multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
) 
multinichenet_output = make_lite_output(multinichenet_output)

saveRDS(multinichenet_output, file = file.path(result_dir, "11.2_multinichenet_MS_duraFibro.rds"))

###### Visualization Circos plot 
# Rank top 20 by prioritization_score
duraFibro_MS_top20 <- duraFibro_MS_interactions %>%
  dplyr::arrange(dplyr::desc(prioritization_score)) %>%
  dplyr::slice_head(n = 20)
senders_receivers <- union(duraFibro_MS_top20$sender %>% unique(), 
                           duraFibro_MS_top20$receiver %>% unique()) %>% sort()
colors_sender <- cluster_colors[senders_receivers]
colors_receiver <- cluster_colors[senders_receivers]

pdf(file.path(result_dir, "11.2_top20_circos_plot_MS_duraFibro.pdf"), width = 10, height = 10)
circos_list <- make_circos_group_comparison(duraFibro_MS_top20, colors_sender, colors_receiver)
dev.off()

#===== 12.0 ROI ======
roi_dir <- #roi information and coordinated from xenium explorer

roi_files <- list.files(
  roi_dir,
  pattern = "_ROI[0-9]*_cells_stats\\.csv$",
  full.names = TRUE
)

basename(roi_files)
length(roi_files)

#conversion because Cell.ID are slightly different in Xenium Explorer and Seurat
meta <- sc_combined@meta.data %>%
  as.data.frame()

meta$seurat_cell <- rownames(meta)
meta$xenium_cell <- sub("_[0-9]+$", "", meta$seurat_cell)

sum("fbipkacg-1" %in% meta$xenium_cell)


#put multiple samples of one sample into one seurat object --> collect ROI cells per sample
roi_cells_by_sample <- list()

for (f in roi_files) {
  
  file_base <- tools::file_path_sans_ext(basename(f))
  
  condition <- stringr::str_extract(file_base, "Ctrl|MS")
  
  sample_part <- file_base %>%
    stringr::str_remove("^S[0-9]+_") %>%
    stringr::str_remove("_Ctrl.*$") %>%
    stringr::str_remove("_MS.*$") %>%
    stringr::str_remove("^s")
  
  sample_part <- gsub("-", "_", sample_part)
  
  sample_key <- paste0(condition, "_", sample_part)
  roi_name <- stringr::str_extract(file_base, "ROI[0-9]*")
  
  roi_df <- read.csv2(f)
  
  roi_cells_xenium <- as.character(roi_df$Cell.ID) %>%
    trimws()
  
  roi_cells_seurat <- meta$seurat_cell[
    meta$xenium_cell %in% roi_cells_xenium &
      meta$sample == sample_key
  ]
  
  message(
    basename(f),
    " | ", sample_key,
    " | ", roi_name,
    ": ",
    length(roi_cells_seurat),
    " cells"
  )
  
  if (length(roi_cells_seurat) == 0) next
  
  roi_cells_by_sample[[sample_key]] <- unique(c(
    roi_cells_by_sample[[sample_key]],
    roi_cells_seurat
  ))
}

# Create one Seurat object per sample, with ROI1/ROI2 merged if present
roi_list <- list()

for (sample_key in names(roi_cells_by_sample)) {
  
  roi_cells <- unique(roi_cells_by_sample[[sample_key]])
  
  roi_obj <- subset(sc_combined, cells = roi_cells)
  
  roi_obj$ROI_sample_key <- sample_key
  
  roi_list[[sample_key]] <- roi_obj
  
  message(sample_key, ": ", length(roi_cells), " merged ROI cells")
}

names(roi_list)
sapply(roi_list, ncol)

# Merge all ROI samples into one object
if (length(roi_list) == 0) {
  stop("No ROI cells matched.")
} else if (length(roi_list) == 1) {
  sc_roi <- roi_list[[1]]
} else {
  sc_roi <- merge(
    roi_list[[1]],
    y = roi_list[-1]
  )
}

sc_roi
table(sc_roi$sample)
table(sc_roi$condition)

saveRDS(
  sc_roi,
  file = file.path(result_dir, "12.0_sc_combined_ROI_only.rds")
)

rm(meta, roi_cells_by_sample, roi_df, roi_list, roi_obj, condition, f, file_base, roi_cells, roi_cells_seurat, roi_cells_xenium, roi_dir, roi_files, roi_name, sample_id, sample_key, sample_part)


#### 12.1 correct saving w/ all info #####
# Combine all ROI cells from all samples
all_roi_cells <- unique(unlist(roi_cells_by_sample))

length(all_roi_cells)

# Subset the ORIGINAL Seurat object only once
# This preserves PCA / UMAP / other reductions
sc_roi <- subset(
  sc_combined,
  cells = all_roi_cells
)

sc_roi

# Check reductions
Reductions(sc_roi)

qsave(
  sc_roi,
  file = file.path(result_dir, "12.1_sc_ROI_only.rds")
)


#==== 13.0 ROI analysis no new clustering ====
#### 13.1 visualize ROI cells in whole datatset ####
# Mark ROI cells in full object
sc_combined$ROI_status <- ifelse(
  colnames(sc_combined) %in% colnames(sc_roi),
  "ROI",
  "non_ROI"
)

table(sc_combined$ROI_status)

p <- DimPlot(
  sc_combined,
  group.by = "ROI_status",
  reduction = "umap",
  cols = c("non_ROI" = "grey85", "ROI" = "darkred"),
  pt.size = 0.2
)
ggsave(file.path(result_dir, "13.1_DimPlot_ROI_seurat.pdf"), p, width = 7, height = 6)

p <- DimPlot(
  sc_combined,
  group.by = "ROI_status",
  reduction = "banksy_umap",
  cols = c("non_ROI" = "grey85", "ROI" = "darkred"),
  pt.size = 0.2
)
ggsave(file.path(result_dir, "13.1_DimPlot_ROI_banksy.pdf"), p, width = 7, height = 6)

split_objects <-  SplitObject(sc_combined, split.by = "sample")

cluster_maps <- lapply(names(split_objects), function(sample_name) {
  
  ImageDimPlot(
    split_objects[[sample_name]],
    group.by = "ROI_status",
    border.size = NA,
    dark.background = FALSE,
    size = 0.5, 
    cols = c("non_ROI" = "grey85", "ROI" = "darkred")
  ) +
    ggtitle(sample_name)
  
})

combined_plot <- wrap_plots(cluster_maps) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right"
  )

ggsave(
  file.path(result_dir, "13.1_ImageDimPlot_ROI_in_sc_combined.pdf"),
  combined_plot,
  width = 12,
  height = 8
)

#### 13.2 spatial distriburion ROI ####
split_objects <-  SplitObject(sc_roi, split.by = "sample")

#seurat_clusters
cluster_maps <- lapply(names(split_objects), function(sample_name) {
  
  ImageDimPlot(
    split_objects[[sample_name]],
    group.by = "seurat_clusters",
    border.size = NA,
    dark.background = FALSE,
    size = 0.5, 
    cols = cluster_colors
  ) +
    ggtitle(sample_name)
  
})

combined_plot <- wrap_plots(cluster_maps) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right"
  )

ggsave(
  file.path(result_dir, "13.2_ImageDimPlot_ROI_seurat_clusters.pdf"),
  combined_plot,
  width = 12,
  height = 8
)

#major_group_mix
cluster_maps <- lapply(names(split_objects), function(sample_name) {
  
  ImageDimPlot(
    split_objects[[sample_name]],
    group.by = "snRNA_major_group_mix",
    border.size = NA,
    dark.background = FALSE,
    size = 0.5, 
    cols = cluster_colors_major_groups_mix
  ) +
    ggtitle(sample_name)
  
})

combined_plot <- wrap_plots(cluster_maps) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right"
  )

ggsave(
  file.path(result_dir, "13.2_ImageDimPlot_ROI_major_group_mix.pdf"),
  combined_plot,
  width = 12,
  height = 8
)


#### 13.3 marker genes ####
sc_roi <- PrepSCTFindMarkers( sc_roi, assay = "SCT", verbose = TRUE)
qsave(sc_roi, file.path(result_dir, "13.3_sc_roi_downstream_analysis.qs"))
Idents(sc_roi) <- "seurat_clusters"
markers_seurat <- FindAllMarkers(object=sc_roi,
                                 test.use="wilcox",
                                 layer="data",
                                 only.pos=FALSE,
                                 max.cells.per.ident=1000)

markers_sorted <- markers_seurat %>%
  arrange(cluster, p_val_adj, desc(avg_log2FC))
marker_list <- split(markers_sorted, markers_sorted$cluster)
write.xlsx(marker_list, file.path(result_dir, "13.3_ROI_Cluster_Markers_seurat_res0.4.xlsx"))

Idents(sc_roi) <- "snRNA_major_group_mix"
markers_major_group_mix <- FindAllMarkers(object=sc_roi,
                                          test.use="wilcox",
                                          layer="data",
                                          only.pos=FALSE,
                                          max.cells.per.ident=1000)

markers_sorted <- markers_major_group_mix %>%
  arrange(cluster, p_val_adj, desc(avg_log2FC))
marker_list <- split(markers_sorted, markers_sorted$cluster)
write.xlsx(marker_list, file.path(result_dir, "13.3_ROI_Cluster_Markers_snRNA_major_group_mix.xlsx"))

#### 13.4 Fibroblasts ####
DefaultAssay(sc_roi) <- "SCT"
Idents(sc_roi) <- "snRNA_major_group_mix"
Fibro_roi_subset <- subset(sc_roi, idents = c("duraFibro", "ABC", "arach_bordFibro"))

#### 13.5 DE analysis #####
##### 13.5.1 DE seurat_clusters ####
DefaultAssay(sc_roi) <- "SCT"

sc_roi$annotated <- as.factor(sc_roi$seurat_clusters) #change later
Idents(sc_roi) <- "annotated"
de_nested_list <- list()

comparisons <- list(
  c("MS", "Ctrl")
)

for (comp in comparisons) {
  ident_1 <- comp[1]
  ident_2 <- comp[2]
  comp_name <- paste0(ident_1, "_vs_", ident_2)
  
  message("Running comparison: ", comp_name)
  de_nested_list[[comp_name]] <- list()
  
  for (i in levels(sc_roi$annotated)) {
    message("  Processing cluster: ", i)
    
    cells_in_cluster <- WhichCells(sc_roi, idents = i)
    cluster_meta <- sc_roi@meta.data[cells_in_cluster, ]
    
    group_counts <- table(cluster_meta$condition)
    
    cluster_name <- paste0("Cluster_", i)
    
    if (all(c(ident_1, ident_2) %in% names(group_counts)) &&
        all(group_counts[c(ident_1, ident_2)] >= 3)) {
      
      de <- tryCatch({
        FindMarkers(
          sc_roi,
          ident.1 = ident_1,
          ident.2 = ident_2,
          subset.ident = i,
          group.by = "condition",
          recorrect_umi = FALSE
        )
      }, error = function(e) {
        message("    ERROR in cluster ", i, ": ", conditionMessage(e))
        return(NULL)
      })
      
      if (!is.null(de) && nrow(de) > 0) {
        de$gene <- rownames(de)
        de <- de[order(de$p_val_adj, -de$avg_log2FC), ]
        de_nested_list[[comp_name]][[cluster_name]] <- de
      } else {
        message("    No DE genes found for cluster ", i)
        de_nested_list[[comp_name]][[cluster_name]] <- data.frame(
          Note = "No DE genes found despite sufficient cells"
        )
      }
      
    } else {
      message("    Skipping cluster ", i, " — too few cells in one or both groups")
      de_nested_list[[comp_name]][[cluster_name]] <- data.frame(
        Note = paste("Too few cells in", ident_1, "or", ident_2)
      )
    }
  }
}

saveRDS(de_nested_list, "13.5.1_ROI_de_nested_list_SCT_condition_seurat_clusters.rds")

# Save each comparison as its own Excel file
for (comp_name in names(de_nested_list)) {
  file_name <- file.path(
    result_dir,
    paste0("13.5.1_DE_SCT_ROI_seurat_clusters_", comp_name, ".xlsx")
  )
  write.xlsx(de_nested_list[[comp_name]], file = file_name)
  message("Saved: ", file_name)
}


#Volcano Plots
de_file <- file.path(result_dir, "13.5.1_DE_SCT_ROI_seurat_clusters_MS_vs_Ctrl.xlsx") 

all_clusters <- getSheetNames(de_file)

#Loop over clusters 
for (cluster in all_clusters) {
  
  message("Processing cluster: ", cluster)
  
  # Read cluster sheet
  de_data <- read.xlsx(de_file, sheet = cluster)
  volcano_genes <- de_data
  
  #remove MT, RPL, RPS genes and Lnc 
  volcano_genes <- volcano_genes %>% filter(!grepl("^MT|^RPL|^RPS|^ZNF|orf|^LINC|^ENSG|\\.", gene))
  #Define sig. up/down-reg genes
  volcano_genes$DE <- "NO"
  volcano_genes$DE[volcano_genes$avg_log2FC > 0.5 & volcano_genes$p_val_adj < 0.001] <- "Up"
  volcano_genes$DE[volcano_genes$avg_log2FC < -0.5 & volcano_genes$p_val_adj < 0.001] <- "Down" 
  #sign. up & down
  volcano_genes$sig <- ifelse(volcano_genes$DE == "NO", "", volcano_genes$gene)
  sig <- volcano_genes$sig
  volcano_genes$sig <- ifelse(volcano_genes$gene %in% sig, volcano_genes$gene, NA)
  volcano_genes$DE <- factor(volcano_genes$DE, levels = c("Down", "NO", "Up"))
  volcano_genes$logP <- -log10(volcano_genes$p_val_adj)
  volcano_genes$logP <- pmin(volcano_genes$logP, 50)
  
  # Count up- and down-regulated genes
  n_up <- sum(volcano_genes$DE == "Up", na.rm = TRUE)
  n_down <- sum(volcano_genes$DE == "Down", na.rm = TRUE)
  
  # Filter only significant DE genes
  sig_genes <- volcano_genes %>% 
    filter(DE != "NO") 
  # Rank first by adjusted p-value, then by absolute log2FC
  top_genes <- sig_genes %>%
    arrange(p_val_adj, desc(abs(avg_log2FC))) %>%  # first p-value, then log2FC
    slice_head(n = 200) %>%                         # top 200 genes
    pull(gene)
  
  # Update sig column to label only top 200 genes
  volcano_genes$sig <- ifelse(volcano_genes$gene %in% top_genes, volcano_genes$gene, NA)
  
  p <- ggplot(volcano_genes, aes(x = avg_log2FC, y = logP, color = DE, label = sig)) +
    geom_vline(xintercept = c(-0.5, 0.5), color = "gray40", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.001), color = "gray40", linetype = "dashed") + #changed for DE_seq2 0.1
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = c("Down" = "black", "NO" = "grey70", "Up" = "#3d70b5"),
                       labels = c("Down", "ns", "Up")) +
    geom_text_repel(max.overlaps = 50, size = 3) +
    theme_classic(base_size = 14) +
    labs(
      x = "log Fold Change",
      y = expression(-log[10]("adjusted p-value")),
      color = "Regulation"
    ) +
    annotate("text", x = max(volcano_genes$avg_log2FC), y = max(volcano_genes$logP),
             label = paste0("Up: ", n_up), hjust = 1, vjust = -0.5, color = "#3d70b5", size = 5) +
    annotate("text", x = min(volcano_genes$avg_log2FC), y = max(volcano_genes$logP),
             label = paste0("Down: ", n_down), hjust = 0, vjust = -0.5, color = "black", size = 5)
  
  # Save per cluster
  out_file <- paste0(file.path(result_dir, "13.5.1_VolcanoPlot_DE_sct_ROI_seurat_cluster_"), cluster, "_DE.pdf")
  ggsave(out_file, plot = p, width = 10, height = 7)
}



##### 13.5.2 DE snRNA_major_group_mix ####
DefaultAssay(sc_roi) <- "SCT"

sc_roi$annotated <- as.factor(sc_roi$snRNA_major_group_mix) #change later
Idents(sc_roi) <- "annotated"
de_nested_list <- list()

comparisons <- list(
  c("MS", "Ctrl")
)

for (comp in comparisons) {
  ident_1 <- comp[1]
  ident_2 <- comp[2]
  comp_name <- paste0(ident_1, "_vs_", ident_2)
  
  message("Running comparison: ", comp_name)
  de_nested_list[[comp_name]] <- list()
  
  for (i in levels(sc_roi$annotated)) {
    message("  Processing cluster: ", i)
    
    cells_in_cluster <- WhichCells(sc_roi, idents = i)
    cluster_meta <- sc_roi@meta.data[cells_in_cluster, ]
    
    group_counts <- table(cluster_meta$condition)
    
    cluster_name <- paste0("Cluster_", i)
    
    if (all(c(ident_1, ident_2) %in% names(group_counts)) &&
        all(group_counts[c(ident_1, ident_2)] >= 3)) {
      
      de <- tryCatch({
        FindMarkers(
          sc_roi,
          ident.1 = ident_1,
          ident.2 = ident_2,
          subset.ident = i,
          group.by = "condition",
          recorrect_umi = FALSE
        )
      }, error = function(e) {
        message("    ERROR in cluster ", i, ": ", conditionMessage(e))
        return(NULL)
      })
      
      if (!is.null(de) && nrow(de) > 0) {
        de$gene <- rownames(de)
        de <- de[order(de$p_val_adj, -de$avg_log2FC), ]
        de_nested_list[[comp_name]][[cluster_name]] <- de
      } else {
        message("    No DE genes found for cluster ", i)
        de_nested_list[[comp_name]][[cluster_name]] <- data.frame(
          Note = "No DE genes found despite sufficient cells"
        )
      }
      
    } else {
      message("    Skipping cluster ", i, " — too few cells in one or both groups")
      de_nested_list[[comp_name]][[cluster_name]] <- data.frame(
        Note = paste("Too few cells in", ident_1, "or", ident_2)
      )
    }
  }
}

saveRDS(de_nested_list, "13.5.2_ROI_de_nested_list_SCT_condition_major_groups_mix.rds")

# Save each comparison as its own Excel file
for (comp_name in names(de_nested_list)) {
  file_name <- file.path(
    result_dir,
    paste0("13.5.2_DE_SCT_ROI_major_groups_mix_", comp_name, ".xlsx")
  )
  write.xlsx(de_nested_list[[comp_name]], file = file_name)
  message("Saved: ", file_name)
}

#Volcano Plots
de_file <- file.path(result_dir, "13.5.2_DE_SCT_ROI_major_groups_mix_MS_vs_Ctrl.xlsx") 

all_clusters <- getSheetNames(de_file)

#Loop over clusters 
for (cluster in all_clusters) {
  
  message("Processing cluster: ", cluster)
  
  # Read cluster sheet
  de_data <- read.xlsx(de_file, sheet = cluster)
  volcano_genes <- de_data
  
  #remove MT, RPL, RPS genes and Lnc 
  volcano_genes <- volcano_genes %>% filter(!grepl("^MT|^RPL|^RPS|^ZNF|orf|^LINC|^ENSG|\\.", gene))
  #Define sig. up/down-reg genes
  volcano_genes$DE <- "NO"
  volcano_genes$DE[volcano_genes$avg_log2FC > 0.5 & volcano_genes$p_val_adj < 0.001] <- "Up"
  volcano_genes$DE[volcano_genes$avg_log2FC < -0.5 & volcano_genes$p_val_adj < 0.001] <- "Down" 
  #sign. up & down
  volcano_genes$sig <- ifelse(volcano_genes$DE == "NO", "", volcano_genes$gene)
  sig <- volcano_genes$sig
  volcano_genes$sig <- ifelse(volcano_genes$gene %in% sig, volcano_genes$gene, NA)
  volcano_genes$DE <- factor(volcano_genes$DE, levels = c("Down", "NO", "Up"))
  volcano_genes$logP <- -log10(volcano_genes$p_val_adj)
  volcano_genes$logP <- pmin(volcano_genes$logP, 50)
  
  # Count up- and down-regulated genes
  n_up <- sum(volcano_genes$DE == "Up", na.rm = TRUE)
  n_down <- sum(volcano_genes$DE == "Down", na.rm = TRUE)
  
  # Filter only significant DE genes
  sig_genes <- volcano_genes %>% 
    filter(DE != "NO") 
  # Rank first by adjusted p-value, then by absolute log2FC
  top_genes <- sig_genes %>%
    arrange(p_val_adj, desc(abs(avg_log2FC))) %>%  # first p-value, then log2FC
    slice_head(n = 200) %>%                         # top 200 genes
    pull(gene)
  
  # Update sig column to label only top 200 genes
  volcano_genes$sig <- ifelse(volcano_genes$gene %in% top_genes, volcano_genes$gene, NA)
  
  p <- ggplot(volcano_genes, aes(x = avg_log2FC, y = logP, color = DE, label = sig)) +
    geom_vline(xintercept = c(-0.5, 0.5), color = "gray40", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.001), color = "gray40", linetype = "dashed") + #changed for DE_seq2 0.1
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = c("Down" = "black", "NO" = "grey70", "Up" = "#3d70b5"),
                       labels = c("Down", "ns", "Up")) +
    geom_text_repel(max.overlaps = 50, size = 3) +
    theme_classic(base_size = 14) +
    labs(
      x = "log Fold Change",
      y = expression(-log[10]("adjusted p-value")),
      color = "Regulation"
    ) +
    annotate("text", x = max(volcano_genes$avg_log2FC), y = max(volcano_genes$logP),
             label = paste0("Up: ", n_up), hjust = 1, vjust = -0.5, color = "#3d70b5", size = 5) +
    annotate("text", x = min(volcano_genes$avg_log2FC), y = max(volcano_genes$logP),
             label = paste0("Down: ", n_down), hjust = 0, vjust = -0.5, color = "black", size = 5)
  
  # Save per cluster
  out_file <- paste0(file.path(result_dir, "13.5.2_VolcanoPlot_DE_sct_ROI_major_group_mix_"), cluster, "_DE.pdf")
  ggsave(out_file, plot = p, width = 10, height = 7)
}

##### 13.5.3 DE Fibroblasts merged ####
Idents(Fibro_roi_subset) <- "condition"

DefaultAssay(Fibro_roi_subset) <- "SCT"
DE_duraFibro_overall <- FindMarkers(
  Fibro_roi_subset,
  ident.1 = "MS",
  ident.2 = "Ctrl",
  group.by = "condition",
  recorrect_umi = FALSE
)

DE_duraFibro_overall <- DE_duraFibro_overall %>%
  rownames_to_column(var = "gene") %>%
  arrange(p_val_adj, desc(avg_log2FC))

write_xlsx(DE_duraFibro_overall, file.path(result_dir, "13.5.3_DE_duraFibro_roi_overall.xlsx"))

### volcano plot
de_file <- file.path(result_dir, "13.5.3_DE_duraFibro_roi_overall.xlsx") 

all_clusters <- getSheetNames(de_file)

#Loop over clusters 
for (cluster in all_clusters) {
  
  message("Processing cluster: ", cluster)
  
  # Read cluster sheet
  de_data <- read.xlsx(de_file, sheet = cluster)
  volcano_genes <- de_data
  
  #remove MT, RPL, RPS genes and Lnc 
  volcano_genes <- volcano_genes %>% filter(!grepl("^MT|^RPL|^RPS|^ZNF|orf|^LINC|^ENSG|\\.", gene))
  #Define sig. up/down-reg genes
  volcano_genes$DE <- "NO"
  volcano_genes$DE[volcano_genes$avg_log2FC > 0.5 & volcano_genes$p_val_adj < 0.001] <- "Up" 
  volcano_genes$DE[volcano_genes$avg_log2FC < -0.5 & volcano_genes$p_val_adj < 0.001] <- "Down" 
  #sign. up & down
  volcano_genes$sig <- ifelse(volcano_genes$DE == "NO", "", volcano_genes$gene)
  sig <- volcano_genes$sig
  volcano_genes$sig <- ifelse(volcano_genes$gene %in% sig, volcano_genes$gene, NA)
  volcano_genes$DE <- factor(volcano_genes$DE, levels = c("Down", "NO", "Up"))
  volcano_genes$logP <- -log10(volcano_genes$p_val_adj)
  volcano_genes$logP <- pmin(volcano_genes$logP, 50)
  
  # Count up- and down-regulated genes
  n_up <- sum(volcano_genes$DE == "Up", na.rm = TRUE)
  n_down <- sum(volcano_genes$DE == "Down", na.rm = TRUE)
  
  # Filter only significant DE genes
  sig_genes <- volcano_genes %>% 
    filter(DE != "NO") 
  # Rank first by adjusted p-value, then by absolute log2FC
  top_genes <- sig_genes %>%
    arrange(p_val_adj, desc(abs(avg_log2FC))) %>%  # first p-value, then log2FC
    slice_head(n = 400) %>%                         # top 500 genes
    pull(gene)
  
  # Update sig column to label only top 400 genes
  volcano_genes$interest <- ifelse(volcano_genes$gene %in% top_genes, volcano_genes$gene, NA)
  
  p <- ggplot(volcano_genes, aes(x = avg_log2FC, y = logP, color = DE, label = sig)) +
    geom_vline(xintercept = c(-0.5, 0.5), color = "gray40", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.001), color = "gray40", linetype = "dashed") + 
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = c("Down" = "black", "NO" = "grey70", "Up" = "#3d70b5"),
                       labels = c("Down", "ns", "Up")) +
    geom_text_repel(max.overlaps = 50, size = 3) +
    theme_classic(base_size = 14) +
    labs(
      x = "log Fold Change",
      y = expression(-log[10]("adjusted p-value")),
      color = "Regulation"
    ) +
    annotate("text", x = max(volcano_genes$avg_log2FC), y = max(volcano_genes$logP),
             label = paste0("Up: ", n_up), hjust = 1, vjust = -0.5, color = "#3d70b5", size = 5) +
    annotate("text", x = min(volcano_genes$avg_log2FC), y = max(volcano_genes$logP),
             label = paste0("Down: ", n_down), hjust = 0, vjust = -0.5, color = "black", size = 5)
  
  # Save per cluster
  out_file <- paste0(file.path(result_dir, "13.5.3_VolcanoPlot_duraFibro_roi_"), cluster, "_DE.pdf")
  ggsave(out_file, plot = p, width = 10, height = 7)
}

#### 13.6 fibrosis score ####
fibrosis <- readr::read_csv(
  file.path(result_dir, "fibrosis.txt"), #genes used in snRNA-seq
  col_names = FALSE,
  show_col_types = FALSE
) %>%
  dplyr::pull(X1) %>%
  unique() %>%
  na.omit()

##### 13.6.1 all clusters ####
#check which genes are present
fibrosis_present <- intersect(fibrosis, rownames(sc_roi))
fibrosis_missing <- setdiff(fibrosis, rownames(sc_roi))
fibrosis_present
fibrosis_missing
length(fibrosis_present)

#calculate fibrosis score 
DefaultAssay(sc_roi) <- "SCT"
fibrosis_score <- AddModuleScore(sc_roi, features = list(fibrosis_present), name = "fibrosisScore", nbin = 5)

df <- FetchData(fibrosis_score, vars = c("fibrosisScore1", "condition", "seurat_clusters", "sample"))
df_avg <- df %>%
  group_by(seurat_clusters, condition, sample) %>%
  summarise(fibrosisScore1 = mean(fibrosisScore1, na.rm = TRUE), .groups = "drop")
df_collapsed <- df_avg %>%
  group_by(sample, condition) %>%
  summarise(fibrosisScore1 = mean(fibrosisScore1), .groups = "drop")

pdf(file.path(result_dir, "13.6.1_Xenium_5k_dura_fibrosis_score_MS_vs_Ctrl_ROI_all_cluster.pdf"), width = 4, height = 6)
ggplot(df_collapsed, aes(x = condition, y = fibrosisScore1, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.6, size = 2, width = 0.2) +
  stat_compare_means(label = "p.signif", method = "wilcox.test") +
  scale_fill_manual(values = c("Ctrl" = "darkgrey", "MS" = "#3d70b5")) +
  theme_classic() +
  labs(title = "fibrosis score all cluster", x = "Condition", y = "average fibrosis score per sample")
dev.off()

##### 13.6.2 fibroblasts ####
#check which genes are present
fibrosis_present <- intersect(fibrosis, rownames(Fibro_roi_subset))
fibrosis_missing <- setdiff(fibrosis, rownames(Fibro_roi_subset))
fibrosis_present
fibrosis_missing
length(fibrosis_present)

#calculate fibrosis score fibroblast cluster
DefaultAssay(sc_roi) <- "SCT"
fibrosis_score <- AddModuleScore(Fibro_roi_subset, features = list(fibrosis_present), name = "fibrosisScore", nbin = 5)

df <- FetchData(fibrosis_score, vars = c("fibrosisScore1", "condition", "seurat_clusters", "sample"))
df_avg <- df %>%
  group_by(seurat_clusters, condition, sample) %>%
  summarise(fibrosisScore1 = mean(fibrosisScore1, na.rm = TRUE), .groups = "drop")
df_collapsed <- df_avg %>%
  group_by(sample, condition) %>%
  summarise(fibrosisScore1 = mean(fibrosisScore1), .groups = "drop")

pdf(file.path(result_dir, "13.6.2_Xenium_5k_dura_fibrosis_score_MS_vs_Ctrl_ROI_fibro_cluster.pdf"), width = 4, height = 6)
ggplot(df_collapsed, aes(x = condition, y = fibrosisScore1, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.6, size = 2, width = 0.2) +
  stat_compare_means(label = "p.signif", method = "wilcox.test") +
  scale_fill_manual(values = c("Ctrl" = "darkgrey", "MS" = "#3d70b5")) +
  theme_classic() +
  labs(title = "fibrosis score all cluster", x = "Condition", y = "average fibrosis score per sample")
dev.off()

#### 13.7 SquidPy distance analysis ####
https://github.com/mihem/dura_human_fibro


