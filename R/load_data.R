#' Load 10X data and GTF annotation
#'
#' @param counts_dir Path to 10X directory (with barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
#' @param gtf_path Path to GTF file (e.g., Homo_sapiens.GRCh38.111.gtf.gz)
#'
#' @return A list with two elements: seurat_obj and gtf
#' @export
load_data_and_gtf <- function(counts_dir, gtf_path) {
  seurat_obj <- Seurat::CreateSeuratObject(counts = Seurat::Read10X(counts_dir))
  gtf <- rtracklayer::import(gtf_path)
  return(list(seurat_obj = seurat_obj, gtf = gtf))
}
