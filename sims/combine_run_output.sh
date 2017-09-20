#!/usr/bin/env bash

#SBATCH -D /home/mstetter/popgen_maize_sim/
#SBATCH -o /home/mstetter/popgen_maize_sim/logs/simulation_out-%j.txt
#SBATCH -e /home/mstetter/popgen_maize_sim/logs/simulation_stderr-%j.txt
#SBATCH -t 8-00:00
#SBATCH -J combineRes
#SBATCH --nodes=1-1
#SBATCH --ntasks 1
#SBATCH --exclude bigmem10

s=neutral2
cat results/tmp_"$s"/sim_neutralPi_European*.csv > results/sim_neutralPi_European_"$s".csv
cat results/tmp_"$s"/sim_neutralPi_maize*.csv > results/sim_neutralPi_maize_"$s".csv

cat results/tmp_"$s"/sim_singleton_European*.csv > results/sim_singleton_European_"$s".csv
cat results/tmp_"$s"/sim_singleton_maize*.csv > results/sim_singleton_maize_"$s".csv

#rm -rf tmp_"$s" 
