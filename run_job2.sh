#!/bin/bash
source ./glpkpython/bin/activate
module load py3-numpy/1.19.9.5
module load pyomo/6.3.0
python $1 $2 $3 $4
deactivate
