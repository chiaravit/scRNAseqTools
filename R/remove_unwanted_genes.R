#' Remove unwanted genes: mitochondrial and ribosomal genes and ribosomal pseudogenes
#'
#' Filters out genes from a GRanges object (already filtered for protein-coding) based on their name or biotype.
#'
#' Specifically removes:
#' - Mitochondrial genes (names starting with "MT-")
#' - Ribosomal protein genes (names starting with "RPS" or "RPL")
#' - Pseudogenes ( "pseudogene" as biotype)
#'
#' @param gtf A GRanges object, typically already filtered for protein-coding genes
#'
#' @return A GRanges object with unwanted gene types removed
#' @export
#'
#' @examples
#' to use the function: gtf_clean <- remove_unwanted_genes(gtf_pc)

remove_unwanted_genes <- function(gtf) {
  keep <- !grepl("^MT-", gtf$gene_name) &           # mitochondrial genes
    !grepl("^RPS", gtf$gene_name) &           # ribosomal protein S
    !grepl("^RPL", gtf$gene_name) &           # ribosomal protein L
    !grepl("pseudogene", gtf$gene_biotype)    # pseudogenes

  gtf[keep]
}
