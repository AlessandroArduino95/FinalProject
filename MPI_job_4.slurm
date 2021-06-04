#!/bin/sh

#SBATCH --job-name PDE_MPI_4
#SBATCH --error err_%j.txt
#SBATCH --output out_%j.txt
#SBATCH --mail-user alessandroarduino95@gmail.com
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --partition allgroups
#SBATCH --ntasks 4
#SBATCH --mem 1G
#SBATCH --time 00:10:00

spack load openmpi@3.1.4

cd $SLURM_SUBMIT_DIR

singularity exec tflow_opencv.sif \
mpirun -np 4 python3 -m mpi4py /FinalProject/MPI_parabol_pde_solver.py