#!/bin/bash

#SBATCH -p cn
#SBATCH -J workflow-1
#SBATCH -o output-%J.o
#SBATCH -n 50
#SBATCH --mem=102400

cd ~/proj/gcap/test-workflow/result

Rscript ~/proj/gcap/test-workflow/one-pair.R
