#!/bin/sh

#SBATCH --job-name serialPDE
#SBATCH --error err_%j.txt
#SBATCH --output out_%j.txt
#SBATCH --mail-user alessandroarduino95@gmail.com
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --partition allgroups
#SBATCH --ntasks 1
#SBATCH --mem 30M
#SBATCH --time 00:10:00

cd $SLURM_SUBMIT_DIR

srun singularity exec tflow_opencv.sif python script.py