#!/usr/bin/env bash

#SBATCH -D /home/mstetter/BGS_sims/
#SBATCH -o /home/mstetter/BGS_sims/logs/simulation_out-%j.txt
#SBATCH -e /home/mstetter/BGS_sims/logs/simulation_stderr-%j.txt
#SBATCH -t 8-00:00
#SBATCH -J bgsSim
#SBATCH --array=0-499
#SBATCH --nodes=1-1
#SBATCH --ntasks 2
#SBATCH --exclude bigmem10



mkdir -p results/tennessen/
python sims/simulate_humans.py demographies/tennessen.csv results/tennessen/ $SLURM_ARRAY_TASK_ID

#mkdir -p results/torres/
#python sims/simulate_humans.py demographies/torres.csv results/torres/ $SLURM_ARRAY_TASK_ID