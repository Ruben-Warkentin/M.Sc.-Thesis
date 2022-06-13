#!/bin/bash
#
#SBATCH --cpus-per-task=16
#SBATCH --mem 8192
#SBATCH -o slurm.%N.%j.out    # STDOUT
#SBATCH -t 0-01:00            # time (D-HH:MM)
#SBATCH --gres=gpu:1
#SBATCH --account=def-dkwan

module load StdEnv/2018.3
module load cuda/10.0.130
module load namd-multicore/2.13
namd2 +p$SLURM_CPUS_PER_TASK  +idlepoll SMD_run.tcl
