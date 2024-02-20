#!/bin/bash
## USAGE: sbatch solo_pipeline.sh <workingdirectory> <input_file_type> <input_file_name> <identifier>

#Input file type is either: 10x or RData
#This script will:
##1. Take an R object and make it an h5AD by calling RDatatoh5ad.R 
### OR ###
##1. Take the 10X object and make it an H5ad by calling 10x_to_h5ad.R
### THEN ###
##2. Calls solo 

#SBATCH -c 1                               # 1 core
#SBATCH -t 0-04:00                         # Runtime of 4 hours, in D-HH:MM format
#SBATCH -p short                           # Run in short partition
#SBATCH -o solo-analysis-pipeline_%j.out           # File to which STDOUT + STDERR will be written, including job ID in filename
#SBATCH --mail-type=END,FAIL                    # Notifies when job ends or if job fails
#SBATCH --mail-user=lauren_wasson@hms.harvard.edu # Email to which notifications will be sent
#SBATCH --mem-per-cpu=8G
#SBATCH -e %j.solo-analysis-pipeline.err


module load R/4.0.1
module load conda2/4.2.13

wd=$1
analysis_type=$2
input_name=$3
identifier=$4

echo "working directory = $wd"
echo "input_type = $input_type"
echo "input_name = $input_name"
echo "identifier = $identifier"

if [[ $analysis_type = "10X" ]]; then
        echo "Running Pipeline from a 10x file"
        cd $wd
	echo "Submitting the following R script:"
	echo "Rscript /n/groups/seidman/lauren/scRNA-seq/scripts/10x_to_h5ad.R $wd $input_name $identifier"
        Rscript /n/groups/seidman/lauren/scRNA-seq/scripts/10x_to_h5ad.R $wd $input_name $identifier 
elif [[ $analysis_type = "RData" ]]; then
        echo "Running Pipeline from an RData file"
        cd $wd
	echo "Submitting the following R script:"
	echo "Rscript /n/groups/seidman/lauren/scRNA-seq/scripts/RDatatoh5ad.R $wd $input_name"
        Rscript /n/groups/seidman/lauren/scRNA-seq/scripts/RDatatoh5ad.R $wd $input_name
else
        echo "Please provide the correct input file type (10X or RData"
fi


h5ad=${identifier}.h5ad
#Run solo
solo /n/data1/hms/genetics/seidman/reichart/h5ad_solo/model_json/model_json ${h5ad} -o /n/groups/seidman/lauren/scRNA-seq/scripts/L592-2_pipelinetest

#Make a CSV of the solo data
#python3 /n/groups/seidman/lauren/scRNA-seq/scripts/s_integration_after_solo.py ${inputdir} ${identifier}

