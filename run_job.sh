#!/bin/bash
#SBATCH --qos=long
source ./pythonLSTM/bin/activate
module load py3-numpy/1.19.0
python $1 $2 $3 $4
deactivate
