#!/bin/bash

#SBATCH --array=0-44
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --output=/sdf/home/p/phil1/code/src/github.com/phil-mansfield/symphony_pipeline/logs/MilkyWay/track/log_%A_%a.out
#SBATCH --error=/sdf/home/p/phil1/code/src/github.com/phil-mansfield/symphony_pipeline/logs/MilkyWay/track/log_%A_%a.err
#SBATCH -p kipac

config=configs/MilkyWay/config.txt
suffix=k64_n30
snap_range=235:235

#python3 find_infall_cores.py ${config} ${SLURM_ARRAY_TASK_ID} &&
python3 print_core_catalogue.py ${config} ${SLURM_ARRAY_TASK_ID} --reset --suffix=${suffix} --snap_range=${snap_range} &&
   python3 convert_core_catalogue.py ${config} ${SLURM_ARRAY_TASK_ID} --suffix=${suffix} &&
   echo "done"
