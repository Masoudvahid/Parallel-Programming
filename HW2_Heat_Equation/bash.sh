#!/bin/bash
#
#SBATCH --nodes=1 # Number of nodes
#SBATCH --ntasks-per-node=10 # Number of processes per node
#SBATCH --partition=RT_study
#SBATCH --job-name=Example
#SBATCH --comment="Run mpi from config"
#SBATCH --output=out.txt
#SBATCH --error=error.txt
mpiexec ./a.out