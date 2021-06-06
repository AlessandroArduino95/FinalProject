#!/bin/sh

#SBATCH --job-name serialPDE
#SBATCH --error err_%x.txt
#SBATCH --output out_%x.txt
#SBATCH --mail-user alessandroarduino95@gmail.com
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --partition allgroups
#SBATCH --ntasks 1
#SBATCH --mem 1G
#SBATCH --time 00:20:00

cd $SLURM_SUBMIT_DIR

srun singularity exec fpcont.sif python3 /FP/FinalProject/parabol_pde_solver.py