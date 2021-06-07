#!/bin/sh

#SBATCH --job-name PDE_MPI_3B
#SBATCH --error err_%x.txt
#SBATCH --output out_%x.txt
#SBATCH --mail-user alessandroarduino95@gmail.com
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --partition allgroups
#SBATCH --ntasks 3
#SBATCH --mem 2G
#SBATCH --time 01:00:00

spack load openmpi@3.1.4

cd $SLURM_SUBMIT_DIR

singularity exec fpcont.sif \
mpirun -np 3 python3 -m mpi4py /FP/FinalProject/MPI_parabol_pde_solverBigger.py