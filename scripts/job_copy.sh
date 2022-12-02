#!/bin/sh
#SBATCH --partition=kipac
#SBATCH -t 8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -J copy_files
#SBATCH --output=log.oe
#SBATCH --mem-per-cpu=2G

python3 copy_sparse_snapshots.py
