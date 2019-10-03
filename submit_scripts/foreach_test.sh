#!/bin/bash -e

# Job Name
#$ -N foreach_test

# Execute the script from the Current Working Directory
#$ -cwd

# Merge the output of the script, and any error messages generated to one file
#$ -j y

# Send the output of the script to a directory called 'UGE-output' uder current
# working directory (cwd)
if [ ! -d "UGE-output" ]; then #Create output directory in case it does NOT
          exist
mkdir UGE-output
fi
#$ -o UGE-output/

# Tell the job your hardware requirements
#$ -l mem_free=100G,h_vmem=100G
#$ -l avx2
#$ -pe threaded 2

# Send mail when the job is submitted, and when the job
# completes
#$ -m be

#  Specify an email address to use
#$ -M rohit.farmer@nih.gov

module load R/3.6.0-goolf-1.7.20

#Rscript --vanilla scripts/script1_arcsinh.R meta/all_phospho_20.yaml
#Rscript --vanilla scripts/script1_1_diag_plots.R meta/all_109_20_1k.yaml
#Rscript --vanilla scripts/script2_clust.R meta/all_phospho_20.yaml
Rscript --vanilla scripts/script2_1_clust_plots.R meta/all_phospho_20.yaml
#Rscript --vanilla scripts/script3_dr.R meta/all_phospho_20.yaml
#Rscript --vanilla scripts/script3_1_dr_plots.R meta/all_109_20_1k.yaml


