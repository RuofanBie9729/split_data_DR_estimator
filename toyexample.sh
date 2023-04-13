#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH -n 1	# Number of cores start small then expand as needed.
#SBATCH -N 1	# Ensure that all cores are on one machine 
#SBATCH -t 7-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -o toyexample.out		# File to which STDOUT will be written 
#SBATCH -e toyexample.err			# File to which STDERR will be written
#SBATCH --mail-type=ALL	# Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ruofan_bie@stat.brown.edu	# Email to which notifications will be sent

cd /home/rbie/split_data_DR
 
R CMD BATCH toyexample.R	# This is your R job/code

