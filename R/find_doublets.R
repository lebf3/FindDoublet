#' Find and Visualize Doublets in Subclusters of a Seurat Object
#'
#' This function performs sub-clustering of cell types within a Seurat object,
#' identifies potential doublets using doublet scores, and generates visualizations
#' for subclusters and doublet scores.
#'
#' @param so A Seurat object containing single-cell RNA-seq data.
#' @param cell_type_col A string specifying the column name in the Seurat object metadata
#'        that contains cell type labels. Default is "cell_types".
#' @param ndim An integer specifying the number of principal components to compute for
#'        clustering and UMAP. Default is 20.
#' @param var_regress A string specifying the variable to regress out using SCTransform.
#'        Default is "percent.mt".
#' @param harmony_vars A character vector specifying the variables to integrate using harmony.
#'        Default is c("sample").
#' @param clust_res A numeric value specifying the resolution for subclustering. Default is 0.1.
#' @param doublet_scores_col A string specifying the column name in the Seurat object metadata
#'        that contains doublet scores from scDblFinder. Default is "scDbt_score".
#' @param output A string specifying the directory to save output plots. Default is "subcluster_doublets".
#'
#' @return A character vector of subcluster identifiers.
#' @export

find_doublets <- function(so = seurat_object,
                          cell_type_col = "cell_types",
                          ndim = 20,
                          var_regress = "percent.mt",
                          harmony_vars = c("sample"),
                          clust_res = 0.1,
                          doublet_scores_col = "scDbt_score",
                          output = "subcluster_doublets") {

  # Validate input parameters
  if (!inherits(so, "Seurat")) stop("The 'so' parameter must be a Seurat object")
  if (!cell_type_col %in% colnames(so@meta.data)) stop("cell_type_col does not exist in the Seurat object metadata")
  if (!doublet_scores_col %in% colnames(so@meta.data)) stop("doublet_scores_col does not exist in the Seurat object metadata")

  set.seed(2024)

  markers.top1 <- presto:::wilcoxauc.Seurat(X = so,
                                            group_by = cell_type_col,
                                            assay = "data",
                                            seurat_assay = "RNA") %>%
    dplyr::filter(pct_out < 50) %>%
    dplyr::group_by(group) %>%
    dplyr::slice_max(auc, n = 1)

  Idents(so) <- so[[cell_type_col]]
  celltypes <- as.character(unique(Idents(so)))

  create_directory <- function(path) {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
    }
  }

  sub_clusters <- lapply(celltypes, function(ct_x){

    cat("Processing cell type:", ct_x, "\n")

    sub <- subset(so, idents = ct_x)

    DefaultAssay(sub) <- "RNA"
    sub <- SCTransform(sub,
                       verbose = FALSE,
                       vars.to.regress = var_regress) %>%
      RunPCA(npcs = ndim) %>%
      harmony::RunHarmony(group.by.vars = harmony_vars, assay.use = "SCT") %>%
      RunUMAP(dims = 1:ndim,
              reduction = "harmony",
              reduction.name = "umap.rna.harmony") %>%
      FindNeighbors(reduction = "harmony", dims = 1:ndim) %>%
      FindClusters(verbose = FALSE, resolution = clust_res)

    ## plot subclusters
    umap_clust <- DimPlot(object = sub,
                          reduction = "umap.rna.harmony",
                          group.by = "seurat_clusters",
                          label = T) +
      ylab(NULL) +
      xlab(NULL) +
      ggtitle("Clusters") +
      NoLegend()

    ## plot doublet scores
    DefaultAssay(sub) <- "RNA"
    umap_db_score <- FeaturePlot(object = sub,
                                 reduction = "umap.rna.harmony",
                                 features = doublet_scores_col) +
      ylab(NULL) +
      xlab(NULL) +
      ggtitle("doublet score RNA")

    vnl_db_score <- VlnPlot(object = sub,
                            features = doublet_scores_col) +
      ggtitle("doublet score RNA")
    patch <- umap_clust | umap_db_score | vnl_db_score

    ## save plots
    # dir.create(output)
    fig.path <- paste0(output, "/", ct_x, "/")
    create_directory(fig.path)
    # dir.create(fig.path)
    ggsave(plot = patch,
           paste0(fig.path, "scAF_Dbl_scores.png"),
           width = 18,
           height = 6)

    ## plot top RNA markers for cell-types on WNN UMAP
    p.features <- FeaturePlot(sub,
                              markers.top1$feature,
                              reduction = "umap.rna.harmony",
                              combine = F)

    p.features2 <- lapply(1:length(p.features), function(x) {
      p.features[[x]] +
        labs(subtitle = markers.top1$group[x]) +
        NoAxes()
    })

    ## save plots
    ggsave(plot = wrap_plots(p.features2),
           paste0(fig.path, "umap_scAF_RNA_markers.png"),
           width = 18,
           height = 12)

    ## return subclusters
    sub$sub.ct.clusters <- paste0(ct_x, "_", sub$seurat_clusters)
    sub$sub.ct.clusters
  })

  unlist(sub_clusters)
}
