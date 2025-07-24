#' Filter protein-coding genes from GTF
#'
#' This function filters a GTF GRanges object to retain only entries that are
#' annotated as protein-coding genes.
#'
#' @param gtf A GRanges object imported from a GTF file.
#'
#' @return A filtered GRanges object containing only protein-coding genes.
#' @export
#'
filter_protein_coding_genes <- function(gtf) {
  # This line filters the GTF and implicitly returns it (or explicitly use 'return()')
  filtered_gtf <- gtf[gtf$type == "gene" & gtf$gene_biotype == "protein_coding"]
  return(filtered_gtf) # Explicitly return the filtered object
}

# The 'length(gtf_pc)' line should NOT be here.
# If you want to print the length, do it in run_all_analysis.R AFTER calling this function.
