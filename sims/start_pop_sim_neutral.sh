#!/usr/bin/env bash

#SBATCH -D /home/mstetter/popgen_maize_sim/
#SBATCH -o /home/mstetter/popgen_maize_sim/logs/simulation_out-%j.txt
#SBATCH -e /home/mstetter/popgen_maize_sim/logs/simulation_stderr-%j.txt
#SBATCH -t 8-00:00
#SBATCH -J sim_neutralMaizeEuro
#SBATCH --array=0-99
#SBATCH --nodes=1-1
#SBATCH --ntasks 1
#SBATCH --exclude bigmem10


mkdir -p results/tmp_neutral2
python sim_maize_europeans_neutral.py $SLURM_ARRAY_TASK_ID


