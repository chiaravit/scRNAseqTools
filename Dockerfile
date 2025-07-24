FROM rocker/r-ver:4.3.1

# STEP 2: Install dependencies for R
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libhdf5-dev \
    libzstd-dev \
    libbz2-dev \
    liblzma-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    build-essential \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# STEP 3: working directory within the container
WORKDIR /app

# STEP 4: copy fikes inside the working directory of the container
COPY . /app/

# STEP 5: Install R packages
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

# Install CRAN packages
RUN R -e "install.packages(c('ggplot2', 'dplyr', 'Matrix', 'patchwork', 'readr', 'uwot', 'Rtsne', 'tibble', 'stringr', 'scales', 'viridisLite'), repos='https://cloud.r-project.org')"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c('Seurat', 'SingleR', 'celldex', 'SummarizedExperiment', 'rtracklayer', 'org.Hs.eg.db', 'AnnotationDbi'))"

# STEP 6: Starting contrainer
CMD ["Rscript", "R/run_all_analysis.R"]
