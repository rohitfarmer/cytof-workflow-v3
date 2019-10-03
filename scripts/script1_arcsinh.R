#!/usr/bin/env Rscript --vanilla

# PURPOSE: To run CyTOF workflow (Robinson's) version 3. This script produces
# arcsinh transformed and scaled data for further analysis.
# PUBLICATION: https://f1000research.com/articles/6-748

# Read command line arguments passed to the script. 
cmd_args = commandArgs(trailingOnly=TRUE)

# Temp. If using RStudio.
# setwd("/Users/farmerr2/locus/sandbox/projects/mini")

# Read yaml file. 
suppressMessages(library(yaml))
cat("Running script1_archsinh.R\n")

if(interactive()){
        cat("Running in interactive mode.\n")
        yam <- read_yaml("meta/phospho_qc.yaml", fileEncoding = "UTF-8") # Change yaml file for interactive execution.
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

# CREATE FOLDERS
loginfo("Checking if the required folders need to be created.")
results_folder <- file.path("results", analysis_name)
figures_folder <- file.path("figures", analysis_name)

if(!dir.exists(results_folder)){
  loginfo("Creating %s folder", results_folder)
  dir.create(results_folder)
} else {
  loginfo("%s folder already exists. Output will be over written.", results_folder)
}

if(!dir.exists(figures_folder)){
  loginfo("Creating %s folder", figures_folder)
  dir.create(figures_folder)
} else {
  loginfo("%s folder already exists. Output will be over written.", figures_folder)
}

# LOAD DATA
# Read experiment metadata.
loginfo("Loading experiment metadata: %s", args_file)
md<-suppressMessages(read_tsv(file.path("meta",args_file)))

# Specify levels for conditions & sample IDs to assure desired ordering.
md$condition <- factor(md$condition, levels = condition_levels)
md$sample_id  <- factor(md$sample_id, levels = md$sample_id[order(md$condition)])
md$file_name <- file.path("data", data_location, md$file_name)

# Define colors for conditions (may not be used in version 3; keep it for now).
# names(color_conditions) <- levels(md$condition)

# Read panel data.
loginfo("Loading panel information: %s", panel_file)
panel <- suppressMessages(read_tsv(file.path("meta", panel_file)))

# Load .fcs files into a flowSet object (not needed in version 3; only do it for debugging).
# library(flowCore)
# data_prefix = data_location
loginfo("Loading fcs files mentioned in the experiment metadata into a flowSet: %s", data_location)
fcs_raw <- read.flowSet(file = md$file_name, transformation = FALSE, truncate_max_range = FALSE)

# Check for file names. They should match to what is in the md$file_name.
ids <- c(keyword(fcs_raw, "FILENAME"))
loginfo("Checking .fcs filenames in flowSet.")
loginfo(ids)

# Spot check that all panel columns are in the flowSet object
loginfo("Spot check that all panel columns are in the flowSet object %s", all(panel$fcs_colname %in% colnames(fcs_raw)))

# Construct daFrame
# daFrame function from CATALYST will read the fcs files from the specified location.
# arcsinh transform and scale the data for subsequent analysis.
# It's source code: https://rdrr.io/bioc/CATALYST/src/R/daFrame.R

# library(CATALYST)
loginfo("Constructing daFrame")
daf  <- daFrame(fcs_raw, panel, md, cols_to_use = panel$fcs_colname)

# Remove fcs raw to save memory.
rm(fcs_raw)

# Save daFrame
loginfo("Saving daFrame with arcsinh transformed data.")
saveRDS(daf, file.path(results_folder, "daf_arcsinh.rds"))

loginfo("Done")

