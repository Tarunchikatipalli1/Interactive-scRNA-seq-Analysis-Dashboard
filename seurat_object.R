

# Install and load required libraries if not already installed
if (!requireNamespace("shiny", quietly = TRUE)) {
  install.packages("shiny")
}
if (!requireNamespace("shinydashboard", quietly = TRUE)) {
  install.packages("shinydashboard")
}
if (!requireNamespace("shinyjs", quietly = TRUE)) {
  install.packages("shinyjs")
}
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}
if (!requireNamespace("shinydashboardPlus", quietly = TRUE)) {
  install.packages("shinydashboardPlus")
}
if (!requireNamespace("shinyWidgets", quietly = TRUE)) {
  install.packages("shinyWidgets")
}

# Load libraries
library(shiny)
library(shinydashboard)
library(shinyjs)
library(Seurat)
library(shinydashboardPlus)
library(shinyWidgets)
library(tidyverse)

# Define correct data directory
data_dir <- "/Users/tarun/Desktop/scRNAseq Shiny/filtered_gene_bc_matrices/hg19"

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = data_dir)

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)


#Quality Control
# Calculate mitochondrial gene percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


# Normalize the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(pbmc), 10)

# Scale the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Perform PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:12)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Perform UMAP
pbmc <- RunUMAP(pbmc, dims = 1:12)

# Assign cluster identities
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# Assign cell types
pbmc[["celltype"]] <- Idents(pbmc)

# Randomly assign metadata: Age, Sex, SampleID
set.seed(123) # for reproducibility
pbmc[["sex"]] <- sample(c("male", "female"), ncol(pbmc), replace = TRUE)
pbmc[["age"]] <- round(runif(ncol(pbmc), 18, 65))
pbmc[["sampleID"]] <- paste0("sample", sample(1:10, ncol(pbmc), replace = TRUE))

# Define the save path
save_path <- "/Users/tarun/Desktop/scRNAseq Shiny/seurat_object.rds"
# Ensure directory exists before saving
dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)

# Save the Seurat object
saveRDS(pbmc, file = save_path)
