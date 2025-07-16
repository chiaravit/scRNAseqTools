#' Plot PCA histogram of first 20 PCs to assess variance
#'
#' This function creates a bar plot showing the proportion of variance explained
#' by the first 20 principal components (PCs) from a Seurat object.
#' Extract standard deviations of PCs, calculate variance explained.
#' Create a data frame of variance explained. (df) keeps only first 20 PCs.
#'
#' ggplot2 library creates the histogram.
#'
#' @param seurat_obj A Seurat object containing PCA results
#' @param num_pcs Integer. The number of top principal components to plot. Default is 20.
#'
#' @return A ggplot object showing the variance explained per PC
#' @export
#'
#' @examples
#' to use the function: load ggplot2 library and plot_pca_variance(seurat_obj)

plot_pca_variance <- function(seurat_obj, num_pcs = 20) {
  std_dev <- seurat_obj[["pca"]]@stdev
  var_explained <- std_dev^2 / sum(std_dev^2)

  df <- data.frame(
    PC = 1:length(var_explained),
    Variance = var_explained
  )
  df <- df[1:num_pcs, ]

  ggplot2::ggplot(df, ggplot2::aes(x = PC, y = Variance)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::labs(
      title = paste("Variance explained by first", num_pcs, "PCs"),
      x = "Principal Component",
      y = "Proportion of Variance Explained"
    ) +
    ggplot2::theme_minimal()
}
