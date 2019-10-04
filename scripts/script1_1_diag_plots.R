#!/usr/bin/env Rscript --vanilla

# PURPOSE: To run CyTOF workflow (Robinson's) version 3. This script plots
# diagnostic plots.
# PUBLICATION: https://f1000research.com/articles/6-748

# Read command line arguments passed to the script. 
cmd_args = commandArgs(trailingOnly=TRUE)

# Temp. If using RStudio.
# setwd("/Users/farmerr2/locus/sandbox/projects/mini")

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
data_type  <- yam$data_type
# color_conditions <- yam$color_conditions

# Load external libraries.
cat("Loading external libraries.\n")
suppressMessages(library(cytofWorkflow))
suppressMessages(library(logging))
suppressMessages(library(readr))
suppressMessages(library(ggplot2))
suppressMessages(library(Cairo))
suppressMessages(library(ComplexHeatmap))

# LOGGING
log_file <- file.path("logs", paste(analysis_name, ".log", sep = ""))
basicConfig(level='FINEST')
addHandler(writeToFile, file=log_file, level='DEBUG')

# Assigning results folders (make sure that the script1 is executed before.)
results_folder <- file.path("results", analysis_name)
figures_folder <- file.path("figures", analysis_name)

# Load daFrame with arcsinh transformed data.
loginfo("Loading daFrame with arcsinh transformed data.")
daf <- readRDS(file.path(results_folder, "daf_arcsinh.rds"))

# DIAGNOSTIC PLOTS
loginfo("DIAGNOSTIC PLOTS")

# Figure 1. Per-sample smoothed densities of marker expression (arcsinh-transformed).
loginfo("Generating figure1.")
# pdf(file = file.path(figures_folder, "figure1.pdf"), paper = "letter")
# Cairo(width = 8.5, height = 11, file = file.path(figures_folder, "figure1.pdf"), type="pdf", pointsize=12, bg = "transparent", canvas = "white", units = "in", dpi = "auto")
fig1  <- plotExprs(daf, color_by = "condition")
logdebug("plotExprs done.")
fig1$facet$params$ncol <- 5
ggsave("figure1.pdf", plot = fig1, device = "pdf", path = figures_folder, width = 12, height = 10, units = "in")
# dev.off()
rm(fig1)
gc()

# Figure 2. Barplot showing the number of cells measured for each sam ple in the PBMC dataset.
loginfo("Generating figure2.")
loginfo("Number of cells per sample.")
loginfo(metadata(daf)$experiment_info$n_cells)
fig2 <- plotCounts(daf, color_by = "condition")
ggsave("figure2.pdf", plot = fig2, device = "pdf", path = figures_folder, width = 12, height = 10, units = "in")
rm(fig2)
gc()

# Figure 3. MDS plot.
loginfo("Generating figure3.")
fig3 <- plotMDS(daf, color_by = "condition")
ggsave("figure3.pdf", plot = fig3, device = "pdf", path = figures_folder, width = 12, height = 12, units = "in")
rm(fig3)
gc()

# Figure 4. Heatmap of the median (archsinh-transformed) marker expression 
# of lineage markers and functional markers across all cells measured for 
# each sample in the PBMC dataset.
# For phenotyping data it will not produce any figure because there are no
# functional markers. 
if(data_type == "stimulation"){
        loginfo("Stimulation data. Generating figure4.")
        pdf(file =  file.path(figures_folder, "figure4.pdf"), width = 11, height = 8.5)
        fig4 <- plotExprHeatmap(daf, bin_anno = TRUE, row_anno = TRUE)
        draw(fig4)
        dev.off()
        rm(fig4)
        gc()
}else if(data_type == "phenotyping"){
        loginfo("Phenotyping data skipping figure 4 generation.")
}

# Figure 5. Non-redundancy scores for each of the markers and all samples in the PBMC dataset.
loginfo("Generating figure5.")
fig5 <- plotNRS(daf, markers = type_markers(daf), color_by = "condition")
ggsave("figure5.pdf", plot = fig5, device = "pdf", path = figures_folder, width = 12, height = 10, units = "in")
rm(fig5)
gc()

loginfo("Done")
