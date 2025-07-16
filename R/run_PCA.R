#' Run basic PCA on a Seurat object
#'
#' This function normalizes the data, identifies variable genes, scales the data,
#' and performs Principal Component Analysis (PCA) on the dataset.
#'
#' @param seurat_obj A Seurat object containing raw count data
#'
#' @return A Seurat object updated with PCA results
#' @export
#'
#' @examples
#' tu use the function and visualize the PCA result:
#' seurat_obj <- run_basic_pca(seurat_obj)
#' print(seurat_obj[["pca"]])
#' ElbowPlot(seurat_obj)

run_basic_pca <- function(seurat_obj) {
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  return(seurat_obj)
}
