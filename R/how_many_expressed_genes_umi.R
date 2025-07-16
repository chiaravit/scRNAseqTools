#' Calculate number of genes expressed per cell (≥3 UMIs)
#'
#' For each cell in a Seurat object, counts how many genes have ≥3 UMIs.
#' create a matrix with TRUE/FALSE for genes with ≥3 UMIs and count how many are them
#'
#' @param seurat_obj A Seurat object.
#'
#' @return A vector with the number of expressed genes per cell.
#' @export
#'
calculate_genes_expressed <- function(seurat_obj) {
  counts_matrix <- seurat_obj[["RNA"]]@counts
  colSums(counts_matrix >= 3)
}
