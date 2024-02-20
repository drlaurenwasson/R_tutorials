#!/bin/bash
#SBATCH -c 1                               # 1 core
#SBATCH -t 0-8:00                         # Runtime of 8 hours, in D-HH:MM format
#SBATCH -p short                           # Run in short partition
#SBATCH -o solo_onejob_%j.out           # File to which STDOUT + STDERR will be written, including job ID in filename
#SBATCH --mail-type=FAIL                    # Notifies when job ends or if job fails
#SBATCH --mail-user=lauren_wasson@hms.harvard.edu # Email to which notifications will be sent
#SBATCH --mem-per-cpu=32G
#SBATCH -e solo_onejob%j.err
#SBATCH -p gpu
#SBATCH --gres=gpu:1

module load conda2/4.2.13

eval "$(conda shell.bash hook)"
conda activate solo

while IFS=, read -r model_json h5ad output
do
solo ${model_json} ${h5ad} -o ${output}
done < keyfile.txt
