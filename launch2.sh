#!/bin/bash

FILE_LIST=("model_1.py" "model_2.py" "lagrangian_by_node.py" "lagrangian_not_node.py")
INSTANCE_LIST=("a280.tsp" "eil51.tsp" "pr439.tsp" "bier127.tsp" "eil101.tsp" "eil76.tsp")
K_LIST=("3" "4" "5")
P_LIST=("1" "1.5" "2" "2.5")
for K in ${K_LIST[*]}
	do
	for FILE in ${FILE_LIST[*]}
	do
		for INSTANCE in ${INSTANCE_LIST[*]}
		do
			for P in ${P_LIST[*]}
			do
				run_job2.sh ${FILE} ${INSTANCE} ${K} ${P} &
			done
		done
		wait
	done
done