#!/usr/bin/env Rscript --vanilla

# PURPOSE: To run CyTOF workflow (Robinson's) version 3. This script plots
# figures related to quality control.
# PUBLICATION: https://f1000research.com/articles/6-748

# Read yaml file. 
suppressMessages(library(yaml))
if(interactive()){
        cat("Running in interactive mode.\n")
        yam <- read_yaml("meta/all_109_20.yaml", fileEncoding = "UTF-8") # Change yaml file for interactive execution.
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
merging1_file  <- yam$merging1_file
# color_conditions <- yam$color_conditions

# Load external libraries.
cat("Loading external libraries.\n")
suppressMessages(library(cytofWorkflow))
suppressMessages(library(logging))
suppressMessages(library(tidyverse))

# LOGGING
log_file <- file.path("logs", paste(analysis_name, ".log", sep = ""))
basicConfig(level='FINEST')
addHandler(writeToFile, file=log_file, level='DEBUG')

# Assigning results folders (make sure that the script1 is executed before.)
results_folder <- file.path("results", analysis_name)
figures_folder <- file.path("figures", analysis_name)

# Load daFrame with dimensionality reduction results.
loginfo("Loading daFrame with dimensionality reduction results.")
daf <- readRDS(file.path(results_folder, "daf_dr_script4.rds"))

# Convert rowData from daFrame to a tibble.
rd <- rowData(daf) %>% as_tibble()

# Fetch consensus cluster codes.
clust_codes <- metadata(daf)$cluster_codes
clust_codes <- select(clust_codes, "som100", meta_string) %>% as_tibble()

# Convert som cluster IDs to meta cluster IDs.
rd$cluster_id <- clust_codes[[meta_string]][match(rd$cluster_id, clust_codes$som100)]

# Summarize data by grouping sample_id and cluster_id and counting 
# the number of cells clustered per cluster per sample.
rd_summary <- summarise(group_by(rd, sample_id, cluster_id), count =n())

# Two ways of creating data frame from the summary table created above.
# 1. Spread tibble with samples per row, cluster id per column and frequency in cells.
sam_som_1 <- rd_summary %>% spread(key = 'cluster_id', value = 'count')  %>% ungroup %>% select(-c(sample_id))
# Transpose data.
# sam_som_t <- t(sam_som)

# 2. Spread tibble with cluster id per row, sample id per column and frquency in cells. 
sam_som_2 <- rd_summary %>% spread(key = 'sample_id', value = 'count') %>% ungroup %>% select(-c(cluster_id))

# Plot a boxplot with samples of the x-axis and frequncy per cluster on the y-axis.
loginfo("Generating sample vs frequency box plot.")
tiff(filename = file.path(figures_folder, "figure_qc_boxplot.tiff"), width = 11, height = 8.5, units = "in", res = 300)
fig_qc_box <- boxplot(sam_som_2, use.cols = TRUE, xlab = "Samples", ylab = "Frequency", outline=FALSE)
dev.off()
       
