#' Run UMAP and clustering on Seurat object
#'
#' This function performs dimensionality reduction and clustering on a Seurat object (with PCA already ran)
#' using the specified number of principal components (10, the most informative). method to visualize in 2D the PCA result.
#'
#' @param seurat_obj A Seurat object with PCA results already computed.
#' @param dims Number of principal components to use (10).
#'
#' @return The Seurat object with UMAP and clustering results added.
#' @export
#'
#' @examples
#' seurat_obj <- run_umap_and_cluster(seurat_obj, dims = 10)
#'DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
#'ggplot2::ggtitle("UMAP of clustered cells")

run_umap_and_cluster <- function(seurat_obj, dims = 10) {
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:dims)
  seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 0.8)
  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:dims)
  return(seurat_obj)
}
