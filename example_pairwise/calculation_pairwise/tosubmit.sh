#!/bin/bash
#SBATCH --account=p32318  ## Required: your allocation/account name, i.e. eXXXX, pXXXX or bXXXX
#SBATCH --partition=short ## Required: (buyin, short, normal, long, gengpu, genhimem, etc)
#SBATCH --time=04:00:00 ## Required: How long will the job need to run (remember different partitions have restrictions on this parameter)
#SBATCH --nodes=1 ## how many computers/nodes do you need (no default)
#SBATCH --ntasks-per-node=8 ## how many cpus or processors do you need on per computer/node (default value 1)
#SBATCH --mem=64G ## how much RAM do you need per computer/node (this affects your FairShare score so be careful to not ask for more than you need))
#SBATCH --job-name=_NAME ## When you run squeue -u  this is how you can identify the job
#SBATCH --output=output.log ## standard out and standard error goes to this file

# A regular comment in Bash
module purge all
module load mpi/openmpi-2.1.1-gcc-5.1.0

mpirun -np 8 /projects/p31819/develop/crystalCF/crystalCF
