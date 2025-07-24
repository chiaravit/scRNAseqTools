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

# --- 1. DATA ---
counts_dir <- "filtered_feature_bc_matrix"
gtf_path <- "Homo_sapiens.GRCh38.111.gtf.gz"

# --- 2.
# 2.1 Data and GTF load
source("R/load_data.R")
message("Step 1/9: Caricamento dati e GTF...")
pipeline_data <- load_data_and_gtf(counts_dir = counts_dir, gtf_path = gtf_path)
seurat_obj <- pipeline_data$seurat_obj
gtf_obj <- pipeline_data$gtf #gtf_obj for GRanges object

# 2.2 Protein-coding genes filtering
source("R/filter_genes.R") # Carica la definizione della funzione filter_protein_coding_genes
message("Step 2/10: Filtraggio geni protein-coding...")
# Call the function and assign its return value to gtf_pc
gtf_pc <- filter_protein_coding_genes(gtf_obj)
message(paste("   - Numero di geni protein-coding dopo il filtro:", length(gtf_pc)))
saveRDS(gtf_pc, "output/gtf_protein_coding.rds")

# 2.3 Remove unwanted genes
source("R/remove_unwanted_genes.R")
message("Step 3/9: Removal unwanted genes...")
gtf_clean <- remove_unwanted_genes(gtf_pc)
genes_to_keep <- gtf_clean$gene_name
seurat_obj <- subset(seurat_obj, features = intersect(rownames(seurat_obj), genes_to_keep))
message(paste("   - Removed genes after filtering:",  nrow(seurat_obj@assays$RNA@counts)))

# 2.4 Removed genes table
source("R/table_removed_genes.R")
message("Step 4/9: Table of removed genes generation...")
summary_removed_genes_df <- summarize_removed_genes(gtf_pc)
print(summary_removed_genes_df)
# Saving table:
write.csv(summary_removed_genes_df, "output/summary_removed_genes.csv", row.names = FALSE)

# 2.5 how many genes are expressed per cell (>=3 UMIs)
source("R/how_many_expressed_genes_umi.R")
message("Step 5/9: Calculating expressed genes...")
genes_expr_3umis <- calculate_genes_expressed(seurat_obj)
# Add information to Seurat object for violin plot generation
seurat_obj$genes_expr_3umis <- genes_expr_3umis

# 2.6 Violin plot generation
source("R/violin_plot.R")
message("Step 6/9:  Violin plot for expressed genes generation...")
violin_plot_obj <- plot_genes_expressed_violin(seurat_obj)
#  Save plot:
ggsave("output/violin_plot_genes_expressed.png", plot = violin_plot_obj, width = 8, height = 6, dpi = 300)

# 2.7 PCA analysis
source("R/run_PCA.R")
message("Step 7/9: PCA execution...")
seurat_obj <- run_basic_pca(seurat_obj)
png("output/pca_elbow_plot.png", width = 800, height = 600, res = 100)
ElbowPlot(seurat_obj, ndims = 30)
dev.off()

# 2.8 Histogram variance first 20 PCs
source("R/histogram_variance_first_20PCs.R")
message("Step 8/9: Generating histogram of variance")
pca_variance_plot_obj <- plot_pca_variance(seurat_obj, num_pcs = 20)
# save plot
ggsave("output/histogram_variance_20PCs.png", plot = pca_variance_plot_obj, width = 9, height = 6, dpi = 300)

# 2.9 UMAP and clustering
source("R/run_UMAP_and_clustering.R")
message("Step 9/9:  UMAP and Clustering...")
seurat_obj <- run_umap_and_cluster(seurat_obj, dims = 10)
# save plot
umap_cluster_plot_obj <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
  ggtitle("UMAP of Clustered Cells")
ggsave("output/umap_clustering.png", plot = umap_cluster_plot_obj, width = 8, height = 7, dpi = 300)

# 2.10 Cell annotation
source("R/UMAP_with_cell_annotation.R")
message("Step 10/9 (continua): Cell annotation on UMAP plot...")
seurat_obj <- annotate_cells_with_SingleR(seurat_obj)
# save plot
umap_annotated_plot_obj <- DimPlot(seurat_obj, group.by = "SingleR_label", label = TRUE) +
  ggtitle("UMAP with Predicted Cell Types (SingleR)")
ggsave("output/umap_cell_annotation_final.png", plot = umap_annotated_plot_obj, width = 9, height = 8, dpi = 300)

# --- 3. Save finale Seurat object  ---
message("Analysis completed. Saving final Seurat object")
saveRDS(seurat_obj, "output/final_seurat_object.rds")

message("All results are saved in the Directory 'output/'.")
