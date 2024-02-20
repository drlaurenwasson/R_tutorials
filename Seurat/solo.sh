#!/bin/bash
#SBATCH -c 1                               # 1 core
#SBATCH -t 0-1:00                         # Runtime of 8 hours, in D-HH:MM format
#SBATCH -p short                           # Run in short partition
#SBATCH -o hostname_sinfo_%j.out           # File to which STDOUT + STDERR will be written, including job ID in filename
#SBATCH --mail-type=FAIL                    # Notifies when job ends or if job fails
#SBATCH --mail-user=lauren_wasson@hms.harvard.edu # Email to which notifications will be sent
#SBATCH --mem-per-cpu=32G
#SBATCH -e solo_%j.err
#SBATCH -p gpu
#SBATCH --gres=gpu:1

input1=$1
input2=$2
input3=$3

module load conda2/4.2.13

eval "$(conda shell.bash hook)"
conda activate solo

solo $input1 $input2 -o $input3
