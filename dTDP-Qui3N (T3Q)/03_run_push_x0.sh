#!/bin/bash
#
#SBATCH --cpus-per-task=16
#SBATCH --mem 8192
#SBATCH -o slurm.%N.%j.out    # STDOUT
#SBATCH -t 2-00:00            # time (D-HH:MM)
#SBATCH --gres=gpu:1
#SBATCH --account=def-dkwan

module load StdEnv/2020
module load cuda/11.0
module load namd-multicore/2.14
namd2 +p$SLURM_CPUS_PER_TASK  +idlepoll 03_NAMD_push_x0.tcl
