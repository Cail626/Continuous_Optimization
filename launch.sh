#!/bin/bash

FILE_LIST=("model_1.py" "model_2.py" "lagrangian_by_node.py" "lagrangian_not_node.py")
INSTANCE_LIST=("custom.tsp")
K_LIST=("3" "4" "5")
P_LIST=("0.5" "1" "1.5" "2")
for FILE in ${FILE_LIST}
	for INSTANCE in ${INSTANCE_LIST[*]}
	do
		for K in ${K_LIST[*]}
		do
			for P in ${P_LIST[*]}
			do
				COMMAND="sbatch run_job.sh ${FILE} ${INSTANCE} ${K_LIST} ${P}"
				while ! ${COMMAND}
				do
					sleep 2
				done
			done
		done
	done
done