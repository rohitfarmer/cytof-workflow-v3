#!/usr/bin/env Rscript --vanilla

# PURPOSE: To run CyTOF workflow (Robinson's) version 3. This script plots
# figures differential cell population abundance analysis.
# PUBLICATION: https://f1000research.com/articles/6-748

# Read command line arguments passed to the script. 
cmd_args = commandArgs(trailingOnly=TRUE)

# Temp. If using RStudio.
# setwd("/Users/farmerr2/locus/sandbox/projects/mini")


# Read yaml file. 
suppressMessages(library(yaml))
if(interactive()){
        cat("Running in interactive mode.\n")
        yam <- read_yaml("meta/all_phospho_20.yaml", fileEncoding = "UTF-8") # Change yaml file for interactive execution.
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
suppressMessages(library(readr))
suppressMessages(library(ggplot2))

# LOGGING
log_file <- file.path("logs", paste(analysis_name, ".log", sep = ""))
basicConfig(level='FINEST')
addHandler(writeToFile, file=log_file, level='DEBUG')

# Assigning results folders (make sure that the script1 is executed before.)
results_folder <- file.path("results", analysis_name)
figures_folder <- file.path("figures", analysis_name)

# Load daFrame with merging1 data.
loginfo("Loading daFrame with merging1 data.")
daf <- readRDS(file.path(results_folder, "daf_dr_script4.rds"))

# Differential cell population abundance.
FDR_cutoff <- 0.5

# Figure 20. Relative abundance of the PBMC populations in each sample 
# (x-axis), in the PBMC dataset, represented with a barplot.
loginfo("Generating figure 20.")
# metadata(daf)$experiment_info$sample_id[1:33]
fig20 <- plotAbundances(daf, k = "meta20", by = "sample_id")
ggsave("figure20.pdf", plot = fig20, device = "pdf", path = figures_folder, width = 36, height = 8.5, units = "in", scale = 1, limitsize = F)
rm(fig20)
gc()

# Figure 21. Relative abundance of the PBMC populations in each sample, 
# in the PBMC dataset, represented with boxplots.
# loginfo("Generating figure 21.")
# fig21 <- plotAbundances(daf, k = "merging1", by = "cluster_id", shape = "patient_id")
# ggsave("figure21.pdf", plot = fig21, device = "pdf", path = figures_folder, width = 11, height = 8.5, units = "in", scale = 1)
# rm(fig21)
# gc()

# Figure 22. DA test results and normalized proportions for PBMC 
# cell populations in BCR/FcR-XL stimulated and unstimulated conditions.
# loginfo("Generating figure 22.")
# ei <- metadata(daf)$experiment_info 
# da_formula1 <- createFormula(ei, cols_fixed = "condition", cols_random = "sample_id")
# da_formula2 <- createFormula(ei, cols_fixed = "condition", cols_random = c("sample_id", "patient_id"))
# 
# contrast <- createContrast(c(0, 1))
# 
# da_res1 <- diffcyt(daf,                                            
#     formula = da_formula1, contrast = contrast,                    
#     analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
#     clustering_to_use = "merging1", verbose = FALSE)               
# da_res2 <- diffcyt(daf,                                            
#     formula = da_formula2, contrast = contrast,
#     analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
#     clustering_to_use = "merging1", verbose = FALSE)   
# 
# rowData(da_res1$res)$p_adj < FDR_cutoff
# rowData(da_res2$res)$p_adj < FDR_cutoff
# 
# pdf(file = file.path(figures_folder, "figure22.pdf"), width = 28, height = 26)
# fig22 <- plotDiffHeatmap(daf, da_res2, th = FDR_cutoff, normalize = TRUE, hm1 = FALSE)
# dev.off()
# rm(fig22)
# gc()

# Differential analysis of marker expression stratified by cell population.

# Figure 23. Median (arcsinh-transformed) expression of 14 signaling markers (x-axis) 
# across the 8 identified PBMC cell populations (individual panels).
# fig23 <- plotMedExprs(daf, k = "merging1", facet = "cluster_id",
#     shape_by = "patient_id")                                
# fig23$facet$params$ncol <- 2                                    
# ggsave("figure23.pdf", plot = fig23, device = "pdf", path = figures_folder, width = 11, height = 8.5, units = "in", scale = 1)
# rm(fig23)
# gc()
