% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_data.R
\name{load_data_and_gtf}
\alias{load_data_and_gtf}
\title{Load 10X data and GTF annotation}
\usage{
load_data_and_gtf(counts_dir, gtf_path)
}
\arguments{
\item{counts_dir}{Path to 10X directory (with barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)}

\item{gtf_path}{Path to GTF file (e.g., Homo_sapiens.GRCh38.111.gtf.gz)}
}
\value{
A list with two elements: seurat_obj and gtf
}
\description{
Load 10X data and GTF annotation
}
