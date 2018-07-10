#!/usr/bin/env bash

#SBATCH -D /home/mstetter/BGS_sims/
#SBATCH -o /home/mstetter/BGS_sims/logs/simulation_out-%j.txt
#SBATCH -e /home/mstetter/BGS_sims/logs/simulation_stderr-%j.txt
#SBATCH -t 14-00:00
#SBATCH -J bgsSim
#SBATCH --array=0-999
#SBATCH --nodes=1-1
#SBATCH --ntasks 1
#SBATCH --exclude bigmem1,bigmem10

echo replicate $SLURM_ARRAY_TASK_ID

#mkdir -p results/tennessen/burnins
#srun python sims/simulate_humans.py demographies/tennessen.csv results/tennessen/ $SLURM_ARRAY_TASK_ID

#mkdir -p results/torres/burnins
#srun python sims/simulate_humans.py demographies/torres.csv results/torres/ $SLURM_ARRAY_TASK_ID

mkdir -p results/maize/burnins
srun python sims/simulate_maize.py demographies/maize.csv results/maize/ $SLURM_ARRAY_TASK_ID

#mkdir -p results/generic/burnins
#srun python sims/simulate_generic.py demographies/generic_models.csv results/generic/ $SLURM_ARRAY_TASK_ID
