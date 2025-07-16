#' Annotate cell types using SingleR and HumanPrimaryCellAtlasData
#'
#' This function uses SingleR to annotate cell types in a Seurat object based on
#' the HumanPrimaryCellAtlasData reference.
#' the function ASSUMES PCA, CLUSTERING and UMAP ALREADY BEEN RUN
#'
#' @param seurat_obj A Seurat object with PCA and clustering already performed.
#' @param use_genes Optional character vector of genes to use (default: VariableFeatures).
#'
#' @return The Seurat object with SingleR annotations added to metadata.
#' @export
#'
#' @examples
#' load libraries and run:
#' library(SingleR)
#' library(celldex)          # Reference datasets
#' library(Seurat)
#' library(SummarizedExperiment)  # Required by SingleR
#' seurat_obj <- annotate_cells_with_SingleR_basic(seurat_obj)
#' DimPlot(seurat_obj, group.by = "SingleR_label", label = TRUE) +
#' ggtitle("UMAP with predicted cell types (SingleR)")

annotate_cells_with_SingleR <- function(seurat_obj, use_genes = NULL) {
  if (!requireNamespace("SingleR", quietly = TRUE) ||
      !requireNamespace("celldex", quietly = TRUE)) {
    stop("packages SingleR and celldex to use this function.")
  }

  # Prepare count data
  counts <- Seurat::GetAssayData(seurat_obj, slot = "data")
  counts <- as.matrix(counts)

  # Create metadata: clusters
  clusters <- seurat_obj$seurat_clusters

  # Load reference dataset
  ref <- celldex::HumanPrimaryCellAtlasData()

  # Run SingleR on clusters
  pred <- SingleR::SingleR(
    test = counts,
    ref = ref,
    labels = ref$label.fine,
    clusters = clusters
  )

  # Map predicted labels back to each cell using cluster ID
  seurat_obj$SingleR_label <- pred$labels[match(seurat_obj$seurat_clusters, rownames(pred))]

  return(seurat_obj)
}
