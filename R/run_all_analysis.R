# R/run_all_analysis.R

# --- 0. Impostazioni Iniziali e Caricamento Librerie ---

# Folder "output" generation to save results.
if (!dir.exists("output")) {
  dir.create("output")
}
message(" 'output/' created")

# Libraries
library(Seurat)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(celldex)
library(SingleR)
library(SummarizedExperiment) # Dipendenza di SingleR/celldex
library(Matrix) # Spesso usata con Seurat, specialmente con Read10X

#' @param counts_dir Path to 10X directory (with barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
#' @param gtf_path Path to GTF file (e.g., Homo_sapiens.GRCh38.111.gtf.gz)
#' @param output_dir Directory where all results will be saved.
#' @export
run_sc_analysis_pipeline <- function(counts_dir, gtf_path, output_dir = "output") {

  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # 1.1 Data and GTF load
  message("Step 1/10: Caricamento dati e GTF...")
  pipeline_data <- load_data_and_gtf(counts_dir = counts_dir, gtf_path = gtf_path)
  seurat_obj <- pipeline_data$seurat_obj
  gtf_obj <- pipeline_data$gtf #gtf_obj for GRanges object

  # 1.2 Protein-coding genes filtering
  message("Step 2/10: Filtering protein-coding genes...")
  gtf_pc <- filter_protein_coding_genes(gtf_obj)
  message(paste("   - protein-coding genes after filtering:", length(gtf_pc)))
  saveRDS(gtf_pc, file.path(output_dir, "gtf_protein_coding.rds"))

  # 1.3 Remove unwanted genes
  message("Step 3/10: Removal unwanted genes...")
  gtf_clean <- remove_unwanted_genes(gtf_pc)
  genes_to_keep <- gtf_clean$gene_name
  seurat_obj <- subset(seurat_obj, features = intersect(rownames(seurat_obj), genes_to_keep))
  message(paste("   - Genes after filtering:", nrow(seurat_obj)))

  # 1.4 Removed genes table
  message("Step 4/10: Table of removed genes generation...")
  summary_removed_genes_df <- summarize_removed_genes(gtf_pc)
  print(summary_removed_genes_df)
  write.csv(summary_removed_genes_df, file.path(output_dir, "summary_removed_genes.csv"), row.names = FALSE)

  # 1.5 how many genes are expressed per cell (>=3 UMIs)
  message("Step 5/10: Calculating expressed genes...")
  genes_expr_3umis <- calculate_genes_expressed(seurat_obj)
  seurat_obj$genes_expr_3umis <- genes_expr_3umis

  # 1.6 Violin plot generation
  message("Step 6/10: Violin plot for expressed genes generation...")
  violin_plot_obj <- plot_genes_expressed_violin(seurat_obj)
  ggsave(file.path(output_dir, "violin_plot_genes_expressed.png"), plot = violin_plot_obj, width = 8, height = 6, dpi = 300)

  # 1.7 PCA analysis
  message("Step 7/10: PCA execution...")
  seurat_obj <- run_basic_pca(seurat_obj)
  png(file.path(output_dir, "pca_elbow_plot.png"), width = 800, height = 600, res = 100)
  ElbowPlot(seurat_obj, ndims = 30)
  dev.off()

  # 1.8 Histogram variance first 20 PCs
  message("Step 8/10: Generating histogram of variance")
  pca_variance_plot_obj <- plot_pca_variance(seurat_obj, num_pcs = 20)
  ggsave(file.path(output_dir, "histogram_variance_20PCs.png"), plot = pca_variance_plot_obj, width = 9, height = 6, dpi = 300)

  # 1.9 UMAP and clustering
  message("Step 9/10: UMAP and Clustering...")
  seurat_obj <- run_umap_and_cluster(seurat_obj, dims = 10)
  umap_cluster_plot_obj <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
    ggplot2::ggtitle("UMAP of Clustered Cells")
  ggsave(file.path(output_dir, "umap_clustering.png"), plot = umap_cluster_plot_obj, width = 8, height = 7, dpi = 300)

  # 1.10 Cell annotation
  message("Step 10/10: Cell annotation on UMAP plot...")
  seurat_obj <- annotate_cells_with_SingleR(seurat_obj)
  umap_annotated_plot_obj <- DimPlot(seurat_obj, group.by = "SingleR_label", label = TRUE) +
    ggplot2::ggtitle("UMAP with Predicted Cell Types (SingleR)")
  ggsave(file.path(output_dir, "umap_cell_annotation_final.png"), plot = umap_annotated_plot_obj, width = 9, height = 8, dpi = 300)

  # --- 2. Save final Seurat object ---
  message("Analysis completed. Saving final Seurat object")
  saveRDS(seurat_obj, file.path(output_dir, "final_seurat_object.rds"))

  message(paste("All results are saved in the Directory '", output_dir, "/'.", sep=""))
}
