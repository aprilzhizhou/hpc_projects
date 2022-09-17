#!/bin/bash
#SBATCH -A e30514
#SBATCH -p short
#SBATCH --job-name="sor"
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH -t 00:05:00
#SBATCH --mem=1
#SBATCH --mail-user=zhizhou2022@u.northwestern.edu
#SBATCH --mail-type=BEGIN,END,FAIL
 
cd $SLURM_SUBMIT_DIR 
module load mpi 
mpirun -np 2 ./diffusion 128 1.8 1e-9 100000
