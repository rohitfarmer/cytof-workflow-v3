#!/usr/bin/env Rscript --vanilla

# PURPOSE: To run CyTOF workflow (Robinson's) version 3. This script plots
# figures for differential state (DS) analysis.
# PUBLICATION: https://f1000research.com/articles/6-748

# Read command line arguments passed to the script. 
cmd_args = commandArgs(trailingOnly=TRUE)

# Temp. If using RStudio.
# setwd("/Users/farmerr2/locus/sandbox/projects/mini")

# Read yaml file. 
suppressMessages(library(yaml))
if(interactive()){
        cat("Running in interactive mode.\n")
        yam <- read_yaml("meta/phospho_noqc_20_10k.yaml", fileEncoding = "UTF-8") # Change yaml file for interactive execution.
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
daf <- readRDS(file.path(results_folder, "daf_merging1.rds"))

# Differential cell population abundance.
FDR_cutoff <- 0.05

# Differential analysis of marker expression stratified by cell population.

# Figure 23. Median (arcsinh-transformed) expression of 14 signaling markers (x-axis) 
# across the 8 identified PBMC cell populations (individual panels).
fig23 <- plotMedExprs(daf, k = "merging1", facet = "cluster_id",
    shape_by = "patient_id")                                
fig23$facet$params$ncol <- 2                                    
ggsave("figure23.pdf", plot = fig23, device = "pdf", path = figures_folder, width = 11, height = 8.5, units = "in", scale = 1)
rm(fig23)
gc()

# Fetch experiment info.
ei <- metadata(daf)$experiment_info 

# Create formula for LMM.
ds_formula1 <- createFormula(ei, cols_fixed = "condition")
ds_formula2 <- createFormula(ei, cols_fixed = "condition", cols_random = "patient_id") 

# Create contrast matrix for GLMM.
contrast <- createContrast(c(0,1))

# Run LMM.
ds_res1 <- diffcyt(daf,                                            
    formula = ds_formula1, contrast = contrast,                    
    analysis_type = "DS", method_DS = "diffcyt-DS-LMM",            
    clustering_to_use = "merging1", verbose = T)    

ds_res2 <- diffcyt(daf,                                            
    formula = ds_formula2, contrast = contrast,                    
    analysis_type = "DS", method_DS = "diffcyt-DS-LMM",            
    clustering_to_use = "merging1", verbose = T)

# Check for clusters below the FDR cuttoff.
table(rowData(ds_res1$res)$p_adj < FDR_cutoff)
table(rowData(ds_res2$res)$p_adj < FDR_cutoff)

# Save p and adj-p values to a TSV file. 
p_lmm1 <- rowData(ds_res1$res)
write.table(p_lmm1, file = file.path("results", analysis_name, "p_ds_lmm1.tsv"), sep = "\t", row.names = F)
p_lmm2 <- rowData(ds_res2$res)
write.table(p_lmm2, file = file.path("results", analysis_name, "p_ds_lmm2.tsv"), sep = "\t", row.names = F)

# Summary table displaying the results (raw and adjusted p-values) 
# together with the observed cell population proportions by sample.
sum_lmm1 <- topTable(ds_res1, order_by = "p_adj", show_meds = TRUE, format_vals = TRUE, digits = 3)
write.table(sum_lmm1, file = file.path("results", analysis_name, "summary_ds_lmm1.tsv"), sep = "\t", row.names = F)
sum_lmm2 <- topTable(ds_res2, order_by = "p_adj", show_meds = TRUE, format_vals = TRUE, digits = 3)
write.table(sum_lmm2, file = file.path("results", analysis_name, "summary_ds_lmm2.tsv"), sep = "\t", row.names = F)

# Figure 24 DS test results and normalized expression of signaling markers in
# PBMC populations in stimulated and unstimulated conditions.
pdf(file = file.path(figures_folder, "figure24_lmm_1.pdf"), width = 28, height = 14)
fig24_1 <- plotDiffHeatmap(daf, ds_res1, top_n = 10, order = TRUE, 
                         th = FDR_cutoff, normalize = TRUE, hm1 = FALSE)
dev.off()
rm(fig24_1)
gc()

pdf(file = file.path(figures_folder, "figure24_lmm_2.pdf"), width = 28, height = 14)
fig24_2 <- plotDiffHeatmap(daf, ds_res2, top_n = 10, order = TRUE, 
                         th = FDR_cutoff, normalize = TRUE, hm1 = FALSE)
dev.off()
rm(fig24_2)
gc()

# Differential analysis of overall marker expression
daf <- mergeClusters(daf, k = "meta20", id = "merging_all",
                     table = data.frame(old_cluster = seq_len(20), new_cluster = "all"))

# Figure 25. Median (arcsinh-transformed) expression of state markers
# calculated from all the cells in a given sample in the PBMC dataset.
fig25 <- plotMedExprs(daf[, state_markers(daf)], shape_by = "patient_id")
fig25$facet$params$ncol <- 3  
ggsave("figure25.pdf", plot = fig25, device = "pdf", path = figures_folder, width = 11, height = 8.5, units = "in", scale = 1)
rm(fig25)
gc()

# Fit linear model.
ds_res3 <- diffcyt(daf, formula = ds_formula1, contrast = contrast,               
                   analysis_type = "DS", method_DS = "diffcyt-DS-LMM",            
                   clustering_to_use = "merging_all", verbose = T)

# fit linear mixed model with patient ID as random effect          
ds_res4 <- diffcyt(daf, formula = ds_formula2, contrast = contrast,                    
                   analysis_type = "DS", method_DS = "diffcyt-DS-LMM",            
                   clustering_to_use = "merging_all", verbose = T)

# Check adjusted p values that are below the FDR cutoff.
table(rowData(ds_res3$res)$p_adj < FDR_cutoff)                     
table(rowData(ds_res4$res)$p_adj < FDR_cutoff)                     

# Save p and adj-p values to a TSV file. 
p_lmm3 <- rowData(ds_res3$res)
write.table(p_lmm3, file = file.path("results", analysis_name, "p_ds_lmm3.tsv"), sep = "\t", row.names = F)
p_lmm4 <- rowData(ds_res4$res)
write.table(p_lmm4, file = file.path("results", analysis_name, "p_ds_lmm4.tsv"), sep = "\t", row.names = F)

# Summary table displaying the results (raw and adjusted p-values) 
# together with the observed cell population proportions by sample.
sum_lmm3 <- topTable(ds_res3, order_by = "p_adj", show_meds = TRUE, format_vals = TRUE, digits = 3)
write.table(sum_lmm3, file = file.path("results", analysis_name, "summary_ds_lmm3.tsv"), sep = "\t", row.names = F)
sum_lmm4 <- topTable(ds_res4, order_by = "p_adj", show_meds = TRUE, format_vals = TRUE, digits = 3)
write.table(sum_lmm4, file = file.path("results", analysis_name, "summary_ds_lmm4.tsv"), sep = "\t", row.names = F)

# Figure 26. DS test results and normalized expression of signaling markers
# calculated over all cells in PBMC populations in BCR/FcR-XL stimulated and
# unstimulated conditions.
pdf(file = file.path(figures_folder, "figure26_lmm_3.pdf"), width = 28, height = 14)
fig26_1 <- plotDiffHeatmap(daf, ds_res3, all = TRUE, hm1 = FALSE)
dev.off()
rm(fig26_1)
gc()

pdf(file = file.path(figures_folder, "figure26_lmm_4.pdf"), width = 28, height = 14)
fig26_2 <- plotDiffHeatmap(daf, ds_res4, all = TRUE, hm1 = FALSE)
dev.off()
rm(fig26_2)
gc()

