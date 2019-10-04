# CyTOF Workflow Version 3
This repository is ready to use code for implementing CyTOF workflow V3 from Mark D. Robinson's group. The code is mentioned in their paper [ref](https://f1000research.com/articles/6-748/v3). I have tried to break the pipeline down into logical chunks to make it easier to run it analysis wise. It also makes it easier to troubleshoot errors and resume analysis from a crash point. The code mentioned in the paper doesn't have commands to save the figures. I am presuming if you are using Rstudio, then the figures will pop up in the figures panel, and from there,you can manually save them. However, if you are running the code in batches, then it would be ideal to have code to save them. The majority of the figures are ggplot objects; however, there few figures that are list type. I have included the code to save the figures accordingly. All the figures are saved through a variable; therefore, no figure panel will popup in Rstudio or an Xserver if running remotely. R 3.6.0 compiled in the Singularity container mentioned below also lacks support for x server. Therefore, writing figures to a file in the only option. I also run R interactively in Vim editor using Nvim-R plugin therefore it's easier for me to switch the plot pop up option off. Figure names are according to the figure numbers in the paper. I would recommend running the figure generation scripts interactively. You may have a different number of samples; therefore, you may want to change the figure dimensions accordingly. Some of the figures, especially the ones that generate maker intensity distributions, require a lot of memory. I run them on a cluster with memory ranging from 64Gb to 500Gb. 

In addition to the figures, the code also saves daFrame object (that carries all the input data and results) in an R data structure file after each major calculation step. To ensure if the calculation breaks, then you don't have to run the entire pipeline from the beginning. You can resume from where the last daFrame is saved. I have also implemented code to extract information from daFrame object and store them in TSV files.

# Requirements for the CyTOF Workflow V3
This implementation of the pipeline requires two meta files (args and panel) describing the experimental design; a YAML file that defines the run specific arguments to the scripts and a manually curated cluster merging file.

1. **A panel file.** It's a tab delimited file with a column for all the metals used,
   associated antigens per metal, and to define marke class; type for lineage
   and state for functional.
  1. An example phenotyping panel is available at meta/panel/phenotyping_panel_v3.txt
  2. An example stimulation panel availabe at meta/panel/phospho_panel_v3.txt
2. **An args file.** It's a tab delimited file with a column for FCS filenames to be
   used in a specific analysis, a number given to each sample as an ID,
   conditions to study, and a patient id. Sample ID is unique to each FCS file
   whereas a patient id can be redundant if more than two samples belong to the
   same subject/patient. 
   * A separate args file need to be created for every distinct analysis.
3. In addition to the panel and args file CHI's implementaion of the workflow
   requires a **YAML file** that specifies the location of data, panel and args
   file, and parameters to be passed to the workflow scripts.
   * A separate YAML file need to be created for every distinct analysis.
4. **Cluster Merging 1 file.** After the meta clustering step a manually curated merging file is required that gives a cell population name to each of the meta clusters.

# Setup an Example Project Environment 
## Clone this Repository
`git clone https://github.com/rohitfarmer/cytof-workflow-v3.git <project_name>`

## Setup Working Directory
After cloning the repository:
```
cd <project_name>
sh setup.sh
```
## Download Singularity Container
Download the Singularity container with the CyTOF workflow version 3 in the project directory.

[https://cloud.sylabs.io/library/_container/5d76b50b2c3454e3496d88c9](https://cloud.sylabs.io/library/_container/5d76b50b2c3454e3496d88c9)  
(Please see the description associated with each container for more details.) 

or

`singularity pull library://rohitfarmer/default/cytof_workflow_v3:<container_unique_id>`  

## Download Example Data
* FCS files: [PBMC8_fcs_files.zip](http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow/PBMC8_fcs_files.zip)
```
cd <project_name>
cd data
wget http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow/PBMC8_fcs_files.zip
unzip PBMC8_fcs_files.zip -d PBMC8_fcs_files
```
**Meta data, panel, and merging file are already in the meta folder in this repository; however, below are the instructions to download them. After download convert all the MS Excel files to TSV text files. I have modified the code to read TSVs than MS Excel files. I find it easier to work with TSV text files in a cross platform environment.** 
* Metadata: [PBMC8_metadata.xlsx](http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow/PBMC8_metadata.xlsx)
```
cd <project_name>
cd meta
wget http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow/PBMC8_metadata.xlsx
```

* Panel (v3): [PBMC8_panel_v3.xlsx](http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow/PBMC8_panel_v3.xlsx)
```
cd <project_name>
cd meta/panel
wget http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow/PBMC8_panel_v3.xlsx
```

* Merging of 20 metaclusters: [PBMC8_cluster_merging1.xlsx](http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow/PBMC8_cluster_merging1.xlsx)
```
cd <project_name>
cd meta
wget http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow/PBMC8_cluster_merging1.xlsx
```

* Merging of 12 metaclusters (v3): [PBMC8_cluster_merging2_v3.xlsx](http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow/PBMC8_cluster_merging2_v3.xlsx)
```
cd <project_name>
cd meta
wget http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow/PBMC8_cluster_merging2_v3.xlsx
```

## Prepare YAML file
Prepare a YAML file per analysis in the meta folder.
```
analysis_name: PBMC8
data_type: stimulation # or phenotyping
data_location: PBMC8_fcs_files # Within data folder.
args_file: PBMC8_metadata.txt # Within meta folder.
panel_file: panel/PBMC8_panel_v3.txt # Within meta folder.
condition_levels:
  - Ref
  - BCRXL
no_of_clusters: 20
meta_string: meta20
tsne_no_cells: 500
umap_no_cells: 1000
merging1_file: PBMC8_cluster_merging2_v3.txt
```

# Scripts to Run in Order
All the scripts can be run interactively by specifying the path to the YAML file at the very beginning of the script. Scripts can also be run with Rscript and passing the YAML file as a command line argument. In the workflow below the script are shown to run through Rscript within the singularity container mentioned above.  

## Data Import, Transformation and Diagnostic Plots (Figures 1 to 5)
```
singularity exec cytof_workflow_v3.sif Rscript --vanilla scripts/script1_arcsinh.R meta/PBMC8.yaml
singularity exec cytof_workflow_v3.sif Rscript --vanilla scripts/script1_1_diag_plots.R meta/PBMC8.yaml
```
**Output Data**
* R data structure with arcsinh transformed data: `results/PBMC8/daf_arcsinh.rds`

## Clustering and Meta Clustering (Figures 6 to 8)
```
singularity exec cytof_workflow_v3.sif Rscript --vanilla scripts/script2_clust.R meta/PBMC8.yaml
singularity exec cytof_workflow_v3.sif Rscript --vanilla scripts/script2_1_clust_plots.R meta/PBMC8.yaml 
```
**Output Data**
* R data structure with added clustering data: `results/PBMC8/daf_clust.rds`

## Dimensionality Reduction (Figures 9 to 14)
```
singularity exec cytof_workflow_v3.sif Rscript --vanilla scripts/script3_dr.R meta/PBMC8.yaml
singularity exec cytof_workflow_v3.sif Rscript --vanilla scripts/script3_1_dr_plots.R meta/PBMC8.yaml
```
**Output Data**
* R data structure with added dimensionality reduction data: `results/PBMC8/daf_dr.rds`

## Cluster Merging 1 (Figures 15 to 18)
```
singularity exec cytof_workflow_v3.sif Rscript --vanilla scripts/script4_merging1.R meta/PBMC8.yaml
```
**Output Data**
* R data structure with added merging 1 data: `results/PBMC8/daf_merging1.rds`

## Differential Abundance (DA) Analysis using GLMM and edgeR (Figures 20 to 22)
```
singularity exec cytof_workflow_v3.sif Rscript --vanilla scripts/script5_da.R meta/PBMC8.yaml
```
**Output Data**
* Figure 20 data: `results/PBMC8/fig20_data.tsv`
* TSV files with p and adjusted p-values from GLMM (Figure 22): `results/PBMC8/p_glmm0.tsv, p_glmm1.tsv, p_glmm2.tsv`
* TSV file with p and adjusted p-values from edger (Figure 22): `results/PBMC8/p_edger.tsv`              
* Summary statistics from GLMM and edge runs (Figure 22): `results/PBMC8/summary_glmm0.tsv, summary_glmm1.tsv, summary_glmm2.tsv, summary_edger.tsv`

## Differential State (DS) Analysis using LMM (Figures 23 to 26)
```
singularity exec cytof_workflow_v3.sif Rscript --vanilla scripts/script6_ds.R meta/PBMC8.yaml
```
* TSV files with p and adjusted p-values from LMM (Figure 24 and 26): `results/PBMC8/p_ds_lmm1.tsv, p_ds_lmm2.tsv, p_ds_lmm3.tsv, p_ds_lmm4.tsv`
* TSV files with summary statistics from LMM (Figure 24 and 26): `results/PBMC8/summary_ds_lmm1.tsv, summary_ds_lmm2.tsv, summary_ds_lmm3.tsv, summary_ds_lmm4.tsv`