import sys

suite = sys.argv[1]
fname = "job_tree_%s.sh" % suite

fp = open(fname, "w+")
print("""#!/bin/sh
#SBATCH --partition=kipac
#SBATCH -t 8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -J %s_make_trees
#SBATCH --output=logs/log.%s_tree.%%j.oe
#SBATCH --mem-per-cpu=32G

config=configs/%s/config.txt
index=-1
go run write_binary_tree.go ${config} ${index} && go run write_tree_header.go ${config} ${index} && go run write_subhalos_file.go ${config} ${index}
""" % (suite, suite, suite), file=fp)
