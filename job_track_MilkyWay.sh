#!/bin/bash

#SBATCH --array=0-44
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/scratch/phil1/logs/SymphonyMilkyWay/log_%A_%a.out
#SBATCH --error=/scratch/phil1/logs/SymphonyMilkyWay/log_%A_%a.err
#SBATCH -p kipac

config=configs/MilkyWay/config.txt
suffix=fid
snap_range=0:235

#python3 find_infall_cores.py ${config} ${SLURM_ARRAY_TASK_ID} &&
#python3 print_core_catalogue.py ${config} ${SLURM_ARRAY_TASK_ID} --reset --suffix=${suffix} --snap_range=${snap_range} &&
python3 print_core_catalogue.py ${config} ${SLURM_ARRAY_TASK_ID} --suffix=${suffix} --snap_range=${snap_range} &&
   python3 convert_core_catalogue.py ${config} ${SLURM_ARRAY_TASK_ID} --suffix=${suffix} &&
   echo "done"
