# CyTOF Workflow Version 3
This repository is ready to use code for implementing CyTOF workflow V3 from Mark D. Robinson's group. The code is mentioned in their paper [ref](https://f1000research.com/articles/6-748/v3). I have tried to break the pipeline down into logical chunks to make it easier to run it analysis wise. It also makes it easier to troubleshoot errors and resume analysis from a crash point. The code mentioned in the paper doesn't have commands to save the figures. I am presuming if you are using Rstudio, then the figures will pop up in the figures panel, and from there,you can manually save them. However, if you are running the code in batches, then it would be ideal to have code to save them. The majority of the figures are ggplot objects; however, there few figures that are list type. I have included the code to save the figures accordingly. All the figures are saved through a variable; therefore, no figure panel will popup in Rstudio or an Xserver if running remotely. R 3.6.0 compiled in the Singularity container mentioned below also lacks support for x server. Therefore, writing figures to a file in the only option. I also run R interactively in Vim editor using Nvim-R plugin therefore it's easier for me to switch the plot pop up option off. Figure names are according to the figure numbers in the paper. I would recommend running the figure generation scripts interactively. You may have a different number of samples; therefore, you may want to change the figure dimensions accordingly. Some of the figures, especially the ones that generate maker intensity distributions, require a lot of memory. I run them on a cluster with memory ranging from 64Gb to 500Gb. 

In addition to the figures, the code also saves daFrame object (that carries all the input data and results) in an R data structure file after each major calculation step. To ensure if the calculation breaks, then you don't have to run the entire pipeline from the beginning. You can resume from where the last daFrame is saved. I have also implemented code to extract information from daFrame object and store them in TSV files.

# Clone this Repository
`git clone https://github.com/rohitfarmer/cytof-workflow-v3.git <project_name>`

# Setup Working Directory
After cloning the repository:
```
cd <project_name>
sh setup.sh
```
# Singularity Container
Singularity container with the CyTOF workflow version 3 can be downloaded from:

[https://cloud.sylabs.io/library/_container/5d76b50b2c3454e3496d88c9](https://cloud.sylabs.io/library/_container/5d76b50b2c3454e3496d88c9)  
(Please see the description associated with each container for more details.) 

or

`singularity pull library://rohitfarmer/default/cytof_workflow_v3`  

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

# Scripts to Run in Order
All the scripts can be run interactively by specifying the path to the YAML file at the very beginning of the script. Scripts can also be run with Rscript and passing the YAML file as a command line argument. In the workflow below the script are shown to run through Rscript within the singularity container mentioned above.  

## Data Import, Transformation and Diagnostic Plots (Figures 1 to 5)
```
singularity exec cytof_workflow_v3.sif Rscript --vanill scripts/script1_arcsinh.R meta/pheno_noqc_20_10k.yaml
singularity exec cytof_workflow_v3.sif Rscript --vanilla scripts/script1_1_diag_plots.R meta/pheno_noqc_20_10k.yaml
```

## Clustering and Meta Clustering (Figures 6 to 8)
```
singularity exec cytof_workflow_v3.sif Rscript --vanilla scripts/script2_clust.R meta/pheno_noqc_20_10k.yaml
singularity exec cytof_workflow_v3.sif Rscript --vanilla scripts/script2_1_clust_plots.R meta/pheno_noqc_20_10k.yaml 
```

## Dimensionality Reduction (Figures 9 to 14)
```
singularity exec cytof_workflow_v3.sif Rscript --vanilla scripts/script3_dr.R meta/pheno_noqc_20_10k.yaml
singularity exec cytof_workflow_v3.sif Rscript --vanilla scripts/script3_1_dr_plots.R meta/pheno_noqc_20_10k.yaml
```

## Cluster Merging 1 (Figures 15 to 18)
```
singularity exec cytof_workflow_v3.sif Rscript --vanilla scripts/script4_merging1.R meta/pheno_noqc_20_10k.yaml
```

## Differential Abundance (DA) Analysis using GLMM and edgeR (Figures 20 to 22)
```
singularity exec cytof_workflow_v3.sif Rscript --vanilla scripts/script5_da.R meta/pheno_noqc_20_10k.yaml
```

## Differential State (DS) Analysis using LMM (Figures 23 to 26)
```
singularity exec cytof_workflow_v3.sif Rscript --vanilla scripts/script6_ds.R meta/pheno_noqc_20_10k.yaml
```