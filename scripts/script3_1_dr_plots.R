#!/usr/bin/env Rscript --vanilla

# PURPOSE: To run CyTOF workflow (Robinson's) version 3. This script plots
# figures related to dimenionality reduction.
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

# Load daFrame with dimensionality reduction results saved in script4.
loginfo("Loading daFrame with dimensionality reduction results saved in script4.")
daf <- readRDS(file.path(results_folder, "daf_dr.rds"))

# Extract marker names. 
marker_names  <- colnames(daf@metadata$SOM_codes)

# VISUAL REPRESENTATION WITH UMAP
loginfo("VISUAL REPRESENTATION WITH UMAP")

# Figure 9. UMAP based on the arcsinh-transformed expression of the lineage markers 
# in the cells from the PBMC dataset.
loginfo("Generating figure9.")
for (i in marker_names){
        logdebug("Generating figure9 for %s", i)
        fig9 <- plotDR(daf, "UMAP", color_by = i)
        ggsave(paste("figure9_", i, ".pdf", sep = ""), plot = fig9, device = "pdf", path = figures_folder, width = 8.5, height = 11, units = "in", scale = 1)
        rm(fig9)
        gc()
}

# Figure 10. t-SNE and UMAP based on the arcsinh-transformed expression of lineage markers 
# in the cells from the PBMC dataset.
loginfo("Generating figure10.")
p1 <- plotDR(daf, "TSNE", color_by = meta_string) +             
  theme(legend.position = "none")                          
p2 <- plotDR(daf, "UMAP", color_by = meta_string)               
lgd <- get_legend(p2)                                        
p2 <- p2 + theme(legend.position = "none")                   
fig10 <- plot_grid(p1, p2, lgd, nrow = 1, rel_widths = c(5, 5, 2))  
ggsave("figure10.pdf", plot = fig10, device = "pdf", path = figures_folder, width = 11, height = 8.5, units = "in", scale = 1)
rm(fig10)
gc()

# Figure 11. UMAP as in Figure 10, but stratified by sample.
## Facet per sample
loginfo("Generating figure11.")
fig11 <- plotDR(daf, "UMAP", color_by = meta_string, facet = "sample_id")
ggsave("figure11.pdf", plot = fig11, device = "pdf", path = figures_folder, width = 40, height = 48, units = "in", scale = 1)
rm(fig11)
gc()

# Figure 12. UMAP as in Figure 10, but stratified by condition.
## Facet per condition
loginfo("Generating figure12.")
fig12 <- plotDR(daf, "UMAP", color_by = meta_string, facet = "condition")
ggsave("figure12.pdf", plot = fig12, device = "pdf", path = figures_folder, width = 18, height = 16, units = "in", scale = 1)
rm(fig12)
gc()

# Figure 13. The 100 SOM codes in the PBMC dataset colored according to the 
# metaclustering with Consen-susClusterPlus into 20 cell populations presented 
# after the dimension reduction with (A) t-SNE and (B) PCA. The SOM codes 
# represent characteristics of the 100 (by default) clusters generated in 
# the first step of the FlowSOM pipeline. The size of the points corresponds 
# to the number of cells that were assigned to a given code.
loginfo("Generating figure13.")
fig13 <- plotCodes(daf, k = meta_string)
ggsave("figure13.pdf", plot = fig13, device = "pdf", path = figures_folder, width = 16, height = 14, units = "in", scale = 1)
rm(fig13)
gc()

# Figure 14. Heatmap of the median marker intensities of lineage markers 
# (left panel) and one signaling marker pS6 (right panel) across the 
# 100 SOM codes in the PBMC dataset.
loginfo("Generating figure14.")
for (i in marker_names){
# tiff(filename = file.path(figures_folder, "figure14.tiff"), width = 8.5, height = 11, units = "in", res = 300)
        logdebug("Generating figure14 for %s", i)
        pdf(file = file.path(figures_folder, paste("figure14_", i, ".pdf", sep = "")), width = 20, height = 14)
        fig14 <- plotClusterHeatmap(daf, hm2 = i, k = "som100", m = meta_string, cluster_anno = FALSE, draw_freqs = TRUE)
        dev.off()
        rm(fig14)
        gc()
}

loginfo("Done")
