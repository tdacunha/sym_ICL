#!/bin/sh
#SBATCH --partition=kipac
#SBATCH --array=0-4
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH -J MilkyWayJR
#SBATCH --output=logs/log.tag.MilkyWayHR.%j.oe
#SBATCH --mem-per-cpu=32G


config=configs/MilkyWayHR/config.txt
go run tag_particles.go ${config} ${SLURM_ARRAY_TASK_ID} &&
   go run xv.go ${config} ${SLURM_ARRAY_TASK_ID} &&
   echo "done"
