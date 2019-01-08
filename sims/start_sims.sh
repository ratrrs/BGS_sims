#!/usr/bin/env bash

#SBATCH -D /group/jrigrp9/stetter_projects/BGS_sims/
#SBATCH -o /group/jrigrp9/stetter_projects/BGS_sims/logs/simulation_out-%j.txt
#SBATCH -e /group/jrigrp9/stetter_projects/BGS_sims/logs/simulation_stderr-%j.txt
#SBATCH -t 14-00:00
#SBATCH -J bgsSim
#SBATCH --array=0-9
#SBATCH --nodes=1-1
#SBATCH --ntasks 1
#SBATCH --exclude bigmem1,bigmem10

echo replicate $SLURM_ARRAY_TASK_ID

#mkdir -p results/tennessen/burnins
#srun python sims/simulate_humans.py demographies/tennessen.csv results/tennessen/ $SLURM_ARRAY_TASK_ID

#mkdir -p results/torres/burnins
#srun python sims/simulate_humans.py demographies/torres.csv results/torres/ $SLURM_ARRAY_TASK_ID

#mkdir -p results/maize/burnins
#srun python sims/simulate_maize.py demographies/maize.csv results/maize/ $SLURM_ARRAY_TASK_ID

mkdir -p results/fixed_n/n200/
srun python sims/fixed_n_sim.py 200 results/fixed_n/n200/ $SLURM_ARRAY_TASK_ID

#mkdir -p results/fixed_n/n400/
#srun python sims/fixed_n_sim.py 400 results/fixed_n/n400/ $SLURM_ARRAY_TASK_ID

#mkdir -p results/fixed_n/n800/
#srun python sims/fixed_n_sim.py 800 results/fixed_n/n800/ $SLURM_ARRAY_TASK_ID

