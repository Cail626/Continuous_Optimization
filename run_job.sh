#!/bin/bash
 #SBATCH --job-name=test
 #SBATCH --ntasks=1
#SBATCH --time=01:30:00
#SBATCH --mem-per-cpu=256
source ./pythonLSTM/bin/activate
module load py3-numpy/1.19.0
python $1 $2 $3 $4
deactivate
