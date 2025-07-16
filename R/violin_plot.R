#' Plot number of genes expressed per cell from previous function
#'
#' Plots a violin plot of the number of genes (≥3 UMIs) per cell to visualize result of previous function.
#' take the ggplot2 library, create a data frame for cells and genes expressed and plot the violin plot.
#'
#' @param seurat_obj A Seurat object with the metadata column `genes_expr_3umis`.
#'
#' @return A ggplot object showing the violin plot.
#' @export
#'
#' @examples
#' plot_genes_expressed_violin(seurat_obj)
plot_genes_expressed_violin <- function(seurat_obj) {
  df_violin <- data.frame(
    cell = colnames(seurat_obj),
    genes_expressed = seurat_obj$genes_expr_3umis
  )

  ggplot(df_violin, aes(x = "", y = genes_expressed)) +
    geom_violin(fill = "lightblue") +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    labs(
      title = "Numero di geni espressi (≥3 UMIs) per cellula",
      y = "Numero di geni",
      x = ""
    ) +
    theme_minimal()
}
