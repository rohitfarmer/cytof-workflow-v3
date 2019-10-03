#!/usr/bin/env Rscript --vanilla

# PURPOSE: To run CyTOF workflow (Robinson's) version 3. This script plots
# figures related to clustering.
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
# color_conditions <- yam$color_conditions

# Load external libraries.
cat("Loading external libraries.\n")
suppressMessages(library(cytofWorkflow))
suppressMessages(library(logging))
suppressMessages(library(readr))
suppressMessages(library(ggplot2))
suppressMessages(library(foreach))
suppressMessages(library(doMC))
registerDoMC(cores=2)

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

# Extract marker names from daf. There are multiple places from where it can be
# extracted. I am extracting it from colnames of SOM_codes matrix.
marker_names  <- colnames(daf@metadata$SOM_codes)

# Figure 6. Heatmap of the meadian marker intensities across the 20 cell populations.
# obtained with FlowSOM after the metaclustering step with ConsensusClusterPlus.
loginfo("Generating figure6.")
pdf(file =  file.path(figures_folder, "figure6.pdf"), width = 11, height = 8.5)
fig6 <- plotClusterHeatmap(daf, hm2 = NULL, k = meta_string, m = NULL, cluster_anno = TRUE, draw_freqs = TRUE) 
dev.off()
rm(fig6)
gc()

# Figure 7. Distributions of marker intensities (arcsinh-transformed) 
# in the 20 cell populations obtained with FlowSOM after the metaclustering step 
# with ConsensusClusterPlus.
loginfo("Generating figure7.")
pdf(file =  file.path(figures_folder, "figure7.pdf"), width = 11, height = 8.5)
fig7 <- plotClusterExprs(daf, k = meta_string, markers = "type")  
ggsave("figure7.pdf", plot = fig7, device = "pdf", path = figures_folder, width = 10, height = 12, units = "in")
dev.off()
rm(fig7)
gc()

# Figure 8. Heatmap of the median marker intensities of the 10 lineage markers
# and one signaling marker (pS6) across the 20 cell populations obtained with 
# FlowSOM after the metaclustering step with ConsensusClusterPlus (PBMC data).
loginfo("Generating figure8.")
for (i in marker_names){
        logdebug("Generating plot for %s", i)
        pdf(file = file.path(figures_folder, paste("figure8_", i, ".pdf", sep = "")), width = 20, height = 8.5)
        fig8 <- plotClusterHeatmap(daf, hm2 = i, k = meta_string, draw_freqs = TRUE)
        logdebug("Done plotClusterHeatmap for %s", i)
        dev.off()
        logdebug("Done dev off for %s", i)
        rm(fig8)
        gc()
}

# z = length(marker_names)
# foreach(n = 1:z) %dopar% {
#         i = marker_names[n]
#         logdebug("Generating plot for %s", i)
#         pdf(file = file.path(figures_folder, paste("figure8_", i, ".pdf", sep = "")), width = 20, height = 8.5)
#         fig8 <- plotClusterHeatmap(daf, hm2 = i, k = meta_string, draw_freqs = TRUE)
#         logdebug("Done plotClusterHeatmap for %s", i)
#         dev.off()
#         logdebug("Done dev off for %s", i)
#         rm(fig8)
#         gc()
# }

        loginfo("Done")
