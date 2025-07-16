#' Summarize number of unwanted genes by category
#'
#' Creates a summary table of how many mitochondrial genes, ribosomal genes, and pseudogenes
#' were identified in the GRanges object before filtering.gives amount og removed genes from each category.
#'
#' @param gtf_pc A GRanges object containing protein-coding genes (before filtering)
#'
#' @return A data.frame with the number of removed genes per category
#' @export
#'
#' @examples
#' to use the function:
#' summary_table <- summarize_removed_genes(gtf_pc)
#' print(summary_table)

summarize_removed_genes <- function(gtf_pc) {
  n_mito <- sum(grepl("^MT-", gtf_pc$gene_name))
  n_rps <- sum(grepl("^RPS", gtf_pc$gene_name))
  n_rpl <- sum(grepl("^RPL", gtf_pc$gene_name))
  n_ribo <- n_rps + n_rpl
  n_pseudo <- sum(grepl("pseudogene", gtf_pc$gene_biotype))

  data.frame(
    Category = c("mitochondrial", "Ribosomal", "Pseudogenes"),
    Removed_genes = c(n_mito, n_ribo, n_pseudo)
  )
}
