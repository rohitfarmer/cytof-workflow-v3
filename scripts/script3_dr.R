#!/usr/bin/env Rscript --vanilla

# PURPOSE: To run CyTOF workflow (Robinson's) version 3. This script carries
# out dimensionality reduction.
# PUBLICATION: https://f1000research.com/articles/6-748

# Read command line arguments passed to the script. 
cmd_args = commandArgs(trailingOnly=TRUE)

# Temp. If using RStudio.
# setwd("/Users/farmerr2/locus/sandbox/projects/mini")

# If using example PBMC8 data.
RNGversion("3.5.3")

# Read yaml file. 
suppressMessages(library(yaml))
if(interactive()){
        cat("Running in interactive mode.\n")
        yam <- read_yaml("meta/PBMC8.yaml", fileEncoding = "UTF-8") # Change yaml file for interactive execution.
}else{
        cat("Running in Rscript mode.\n")
        if (length(cmd_args) < 1){
        cat("Missing command line argument(s).\n")
        cat("Usage: script.R name.yaml\n")
        stopifnot(length(cmd_args) > 1)
        }else{
                yam <- read_yaml(cmd_args[1], fileEncoding = "UTF-8") # Pass yaml file as a command line argument.
        }
}

analysis_name <- yam$analysis_name
data_location <- yam$data_location
args_file <- yam$args_file
panel_file <- yam$panel_file
condition_levels <- yam$condition_levels
no_of_clusters <- yam$no_of_clusters
tsne_no_cells <- yam$tsne_no_cells
umap_no_cells <- yam$umap_no_cells
meta_string <- yam$meta_string
# color_conditions <- yam$color_conditions

# Load external libraries.
cat("Loading external libraries.\n")
suppressMessages(library(cytofWorkflow))
suppressMessages(library(logging))
suppressMessages(library(readr))
suppressMessages(library(ggplot2))

# LOGGING
log_file <- file.path("logs", paste(analysis_name, ".log", sep = ""))
basicConfig(level='FINEST')
addHandler(writeToFile, file=log_file, level='DEBUG')

# Assigning results folders (make sure that the script1 is executed before.)
results_folder <- file.path("results", analysis_name)
figures_folder <- file.path("figures", analysis_name)

# Load daFrame with clustering results saved in script2.
loginfo("Loading daFrame with clustering results saved in script2.")
daf <- readRDS(file.path(results_folder, "daf_clust.rds"))

# run t-SNE & UMAP                           
set.seed(1234)          
loginfo("Calculating t-SNE and UMAP.")
daf <- runDR(daf, "TSNE", rows_to_use = tsne_no_cells) 
daf <- runDR(daf, "UMAP", rows_to_use = umap_no_cells)

# Save daFrame 
loginfo("Saving daFrame with DR data.")
saveRDS(daf, file.path(results_folder, "daf_dr.rds"))

loginfo("Done")
