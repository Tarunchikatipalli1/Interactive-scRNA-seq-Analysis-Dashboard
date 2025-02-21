## Interactive scRNA-seq Analysis Dashboard

### Overview

This Shiny app provides an interactive dashboard for analyzing single-cell RNA sequencing (scRNA-seq) data using the Seurat package in R. It allows users to upload a Seurat object, generate UMAP visualizations, explore gene expression, and download analysis results.

### Features

- Upload a Seurat object (.rds file)

- Generate UMAP plots colored by metadata columns

- Visualize expression of selected genes

- Download UMAP and gene expression plots

- Reset the app to start a new analysis

### Installation

To run the app locally, follow these steps:

1. Clone the repository:

```bash
git clone https://github.com/Tarunchikatipalli1/Interactive-scRNA-seq-Analysis-Dashboard.git
cd Interactive-scRNA-seq-Analysis-Dashboard
```

2. Open R and install required dependencies if not already installed:

```bash
install.packages(c("shiny", "shinydashboard", "shinyjs", "Seurat", "shinydashboardPlus", "shinyWidgets", "dplyr"))
```

3. Run the Shiny app:

```bash
library(shiny)
runApp(".")
```

### How to Use This App

Upload Seurat Object:

- Use the "Upload File" button to upload a Seurat object (.rds format).

Run Analysis:

- Press "Run" to start processing the data.

- View UMAP plots and gene expression features.

- Select a metadata column to color the UMAP plot.

- Select a gene to visualize its expression.

Download Results:

- Use the download buttons to save UMAP and gene expression plots.


### File Structure

Shinyapp.R: Defines the UI and server logic for the app.

config.R: Contains global configurations and functions for loading and processing Seurat objects.

seurat_object.R: Sample script for creating a Seurat object from raw data.

### Acknowledgments

This project is built using the Seurat package and the Shiny framework for interactive web applications in R.

License

This project is open-source and available under the MIT License.
