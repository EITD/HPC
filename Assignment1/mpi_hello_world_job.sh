#!/bin/bash -l
# The -l above is required to get the full environment with modules

# The name of the script is mpi_hello_world_job
#SBATCH -J mpi_hello_world_job
# Only 1 hour wall-clock time will be given to this job
#SBATCH -t 1:00:00
#SBATCH -A edu24.DD2356
# Number of nodes
#SBATCH -p main
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
#SBATCH -e error_file.e

# Run the executable file 
# and write the output into my_output_file
srun -n 128 ./a.out > mpi_hello_world_output