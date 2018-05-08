#!/usr/bin/env bash

#SBATCH -D /home/mstetter/BGS_sims/
#SBATCH -o /home/mstetter/BGS_sims/logs/simulation_out-%j.txt
#SBATCH -e /home/mstetter/BGS_sims/logs/simulation_stderr-%j.txt
#SBATCH -t 14-00:00
#SBATCH -J bgsSim
#SBATCH --array=0-1999
#SBATCH --nodes=1-1
#SBATCH --ntasks 1
#SBATCH --exclude bigmem10



#mkdir -p results/tennessen/
#srun python sims/simulate_real_life.py demographies/tennessen.csv results/tennessen/ human $SLURM_ARRAY_TASK_ID

#mkdir -p results/torres/
#srun python sims/simulate_real_life.py demographies/torres.csv results/torres/ human $SLURM_ARRAY_TASK_ID

#mkdir -p results/maize/
#srun python sims/simulate_real_life.py demographies/maize.csv results/maize/ maize $SLURM_ARRAY_TASK_ID

mkdir -p results/generic/
srun python sims/simulate_generic.py demographies/generic_models.csv results/generic/ $SLURM_ARRAY_TASK_ID
