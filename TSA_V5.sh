#!/bin/bash

#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node 
#SBATCH --error=job.%J.err 
#SBATCH --output=job.%J.out 
#SBATCH --mail-user=yfs2333333@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
module load r/4.0.4-py3-4khjixy
Rscript TSA_alpha_V5.R --r $1 > NMA_$SLURM_JOBID.Rout


