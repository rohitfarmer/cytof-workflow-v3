#!/usr/bin/env Rscript --vanilla

# PURPOSE: To run CyTOF workflow (Robinson's) version 3. This script plots
# figures related to differential cell population abundance analysis.
# PUBLICATION: https://f1000research.com/articles/6-748

# Read command line arguments passed to the script. 
cmd_args = commandArgs(trailingOnly=TRUE)

# Temp. If using RStudio.
# setwd("/Users/farmerr2/locus/sandbox/projects/mini")

# Read yaml file. 
suppressMessages(library(yaml))
if(interactive()){
        cat("Running in interactive mode.\n")
        yam <- read_yaml("meta/detector_qc_20.yaml", fileEncoding = "UTF-8") # Change yaml file for interactive execution.
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

# Load external libraries.
cat("Loading external libraries.\n")
suppressMessages(library(cytofWorkflow))
suppressMessages(library(logging))
suppressMessages(library(readr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))

# LOGGING
log_file <- file.path("logs", paste(analysis_name, ".log", sep = ""))
basicConfig(level='FINEST')
addHandler(writeToFile, file=log_file, level='DEBUG')

# Assigning results folders (make sure that the script1 is executed before.)
results_folder <- file.path("results", analysis_name)
figures_folder <- file.path("figures", analysis_name)

# Load daFrame with merging1 data.
loginfo("Loading daFrame with merging1 data.")
daf <- readRDS(file.path(results_folder, "daf_clust.rds"))

# Differential cell population abundance.
# FDR_cutoff <- 0.05

# Figure 20. Relative abundance of the PBMC populations in each sample 
# (x-axis), in the PBMC dataset, represented with a barplot.
loginfo("Generating figure 20.")
# metadata(daf)$experiment_info$sample_id[1:33]
fig20 <- plotAbundances(daf, k = "meta20", by = "sample_id")
ggsave("figure20.pdf", plot = fig20, device = "pdf", path = figures_folder, width = 11, height = 8.5, units = "in", scale = 1, limitsize = F)
fig20_data  <- ggplot_build(fig20)$plot$data %>% 
        select(c(sample_id, cluster_id, patient_id, condition,freq)) %>%
        spread(cluster_id, freq)
write.table(fig20_data, file.path("results", analysis_name, "fig20_data.tsv"), sep = "\t", row.names = F)
rm(fig20)
gc()

# Figure 21. Relative abundance of the PBMC populations in each sample, 
# in the PBMC dataset, represented with boxplots.
loginfo("Generating figure 21.")
fig21 <- plotAbundances(daf, k = "meta20", by = "cluster_id", shape = "patient_id")
ggsave("figure21.pdf", plot = fig21, device = "pdf", path = figures_folder, width = 11, height = 8.5, units = "in", scale = 1)
rm(fig21)
gc()

# Figure 22. DA test results and normalized proportions for PBMC 
# cell populations in BCR/FcR-XL stimulated and unstimulated conditions.
# loginfo("Generating figure 22.")
# 
# Fetch experiment info.
# ei <- metadata(daf)$experiment_info 
# 
# Create formula for DA using GLMM.
# da_formula0  <- createFormula(ei, cols_fixed = "condition", cols_random = "patient_id")
# da_formula1 <- createFormula(ei, cols_fixed = "condition", cols_random = "sample_id")
# da_formula2 <- createFormula(ei, cols_fixed = "condition", cols_random = c("sample_id", "patient_id"))
# 
# Create contrast matrix for GLMM.
# contrast <- createContrast(c(0,1))
# 
# Run DA for GLMM.
# da_res0 <- diffcyt(daf,                                            
#     formula = da_formula0, contrast = contrast,                    
#     analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
#     clustering_to_use = "merging1", verbose = T)               
# 
# da_res1 <- diffcyt(daf,                                            
#     formula = da_formula1, contrast = contrast,                    
#     analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
#     clustering_to_use = "merging1", verbose = T)               
# 
# da_res2 <- diffcyt(daf,                                            
#     formula = da_formula2, contrast = contrast,
#     analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
#     clustering_to_use = "merging1", verbose = T)   
# 
# Check for clusters below the FDR cuttoff. 
# table(rowData(da_res0$res)$p_adj < FDR_cutoff)
# table(rowData(da_res1$res)$p_adj < FDR_cutoff)
# table(rowData(da_res2$res)$p_adj < FDR_cutoff)
# 
# Save p and adj-p values to a TSV file. 
# p_glmm0 <- rowData(da_res0$res)
# write.table(p_glmm0, file = file.path("results", analysis_name, "p_glmm0.tsv"), sep = "\t", row.names = F)
# p_glmm1 <- rowData(da_res1$res)
# write.table(p_glmm1, file = file.path("results", analysis_name, "p_glmm1.tsv"), sep = "\t", row.names = F)
# p_glmm2 <- rowData(da_res2$res)
# write.table(p_glmm2, file = file.path("results", analysis_name, "p_glmm2.tsv"), sep = "\t", row.names = F)
# 
# Summary table displaying the results (raw and adjusted p-values) 
# together with the observed cell population proportions by sample.
# sum_glmm0 <- topTable(da_res0, show_props = TRUE, format_vals = TRUE, digits = 2)
# write.table(sum_glmm0, file = file.path("results", analysis_name, "summary_glmm0.tsv"), sep = "\t", row.names = F)
# sum_glmm1 <- topTable(da_res1, show_props = TRUE, format_vals = TRUE, digits = 2)
# write.table(sum_glmm1, file = file.path("results", analysis_name, "summary_glmm1.tsv"), sep = "\t", row.names = F)
# sum_glmm2 <- topTable(da_res2, show_props = TRUE, format_vals = TRUE, digits = 2)
# write.table(sum_glmm2, file = file.path("results", analysis_name, "summary_glmm2.tsv"), sep = "\t", row.names = F)
# 
# Create design matrix for edgeR method.
# design <- createDesignMatrix(ei, cols_design = c("condition"))
# 
# Create contrast matrix for edgeR.
# design_col  <- ncol(design) - 2
# contrast_edger <- createContrast(c(0,1))
# 
# Check if the nuber of rows in the contrast matrix is equal to the number of
# columns in the desing matrix. 
# nrow(contrast_edger) == ncol(design)
# data.frame(parameters = colnames(design), contrast_edger)
# 
# Run DA for edgeR.
# da_res_edger <- diffcyt(daf, design = design, analysis_type = "DA", contrast = contrast_edger, 
#                         method_DA = "diffcyt-DA-edgeR", clustering_to_use = "merging1",
#                         verbose = T)
# 
# Check for clusters below the FDR cuttoff. 
# table(rowData(da_res_edger$res)$p_adj < FDR_cutoff)
# 
# Save p and adj-p values to a TSV file. 
# p_edger <- rowData(da_res_edger$res)
# write.table(p_edger, file = file.path("results", analysis_name, "p_edger.tsv"), sep = "\t", row.names = F)
# 
# Summary table displaying the results (raw and adjusted p-values) 
# together with the observed cell population proportions by sample.
# sum_edger <- topTable(da_res_edger, show_props = TRUE, format_vals = TRUE, digits = 2)
# write.table(sum_edger, file = file.path("results", analysis_name, "summary_edger.tsv"), sep = "\t", row.names = F)
# 
# Save figure 22 heatmaps. Name them differently depending up the set of data
# used. "all" means all the conditons used. "cc" means only Case and Control
# samples used. 
# 
# Figure 22 for GLMM.
# pdf(file = file.path(figures_folder, "figure22_glmm_0_cc.pdf"), width = 28, height = 26)
# fig22_0 <- plotDiffHeatmap(daf, da_res0, th = FDR_cutoff, normalize = TRUE, hm1 = FALSE)
# dev.off()
# rm(fig22_0)
# gc()
# 
# pdf(file = file.path(figures_folder, "figure22_glmm_1_cc.pdf"), width = 28, height = 26)
# fig22_1 <- plotDiffHeatmap(daf, da_res1, th = FDR_cutoff, normalize = TRUE, hm1 = FALSE)
# dev.off()
# rm(fig22_1)
# gc()
# 
# pdf(file = file.path(figures_folder, "figure22_glmm_2_cc.pdf"), width = 28, height = 26)
# fig22_2 <- plotDiffHeatmap(daf, da_res2, th = FDR_cutoff, normalize = TRUE, hm1 = FALSE)
# dev.off()
# rm(fig22_2)
# gc()
# 
# Figure 22 for edgeR.
# pdf(file = file.path(figures_folder, "figure22_edgeR_cc.pdf"), width = 28, height = 26)
# fig22_edger <- plotDiffHeatmap(daf, da_res_edger, th = FDR_cutoff, normalize = TRUE, hm1 = FALSE)
# dev.off()
# rm(fig22_edger)
# gc()
# 
