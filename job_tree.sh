#!/bin/sh
#SBATCH --partition=kipac
#SBATCH -t 8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -J MilkyWay_make_trees
#SBATCH --output=logs/log.MilkyWay.%j.oe
#SBATCH --mem-per-cpu=32G

config=configs/MilkyWay/config.txt
index=-1
go run write_binary_tree.go ${config} ${index} && go run write_tree_header.go ${config} ${index} && go run write_subhalos_file.go ${config} ${index}
