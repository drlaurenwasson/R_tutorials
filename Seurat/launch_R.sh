#!/bin/bash
#SBATCH -c 2                    # Number of cores requested
#SBATCH -t 6:00:00                    # Runtime in minutes
                                # Or use HH:MM:SS or D-HH:MM:SS, instead of just number of minutes
#SBATCH -p short                # Partition (queue) to submit to
#SBATCH --mem-per-cpu=64G        # 8 GB memory needed (memory PER CORE)
#SBATCH --open-mode=append      # append adds to outfile, truncate deletes first
### In filenames, %j=jobid, %a=index in job array
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-type=END         # Mail when the job ends  
#SBATCH --mail-user=lauren_wasson@hms.harvard.edu   # Email to which notifications will be sent
#write command-line commands below this line
module load R/4.0.1

#Rscript 2021_10_06_Object.R
#Rscript 2021_10_07_Subset.R
#Rscript 2020_07_13_CHD7_FeaturePlots.R
#Rscript 2021_10_21_PediatricDCM_Recluster.R
#Rscript 2021_10_22_HCA_subset_Recluster.R
#Rscript 2021_10_25_noRNAlater_subset_Recluster.R
Rscript Integrate_solo_pbmc.R /n/groups/seidman/lauren/scRNA-seq/CHD/v7 /n/groups/seidman/lauren/scRNA-seq/CHD/v7/Noonan_keyfile.txt 
