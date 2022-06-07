import math
import pyomo.environ as pyo
import numpy as np
import time
import sys
import os
global n, K, P

def get_results(filename):
  with open(filename,'r') as f:
        line = f.read()
        value = line.split(" ")[0]
        time = line.split(" ")[1]
        return value, time

model_list = ["model_1.py", "model_2.py", "lagrangian_by_node.py", "lagrangian_not_node.py"]    
instance_list = ["eil51.tsp", "eil76.tsp", "bier127.tsp", "a280.tsp", "pr439.tsp", "eil101.tsp"]
K_list = ["3", "4", "5"]
P_list =  ["1", "1.5," "2", "2.5", "3"]
value = [ [ [ [ [] for m in range(len(model_list)) ] for k in range(len(P_list))] for j in range(len(K_list)) ] for i in range(len(instance_list))]
time = [ [ [ [ [] for m in range(len(model_list)) ] for k in range(len(P_list))] for j in range(len(K_list)) ] for i in range(len(instance_list))]

for instance in range(len(instance_list)):
    for K in range(len(K_list)):
        for P in range(len(P_list)):
            for model in range(len(model_list)):
                filename = "result"+os.sep+model_list[model]+"_"+instance_list[instance].split('.')[0]+"_"+str(K_list[K])+"_"+str(P_list[P])+".txt"
                value[instance[K[P[model]]]], time[instance[K[P[model]]]] = get_results(filename)
                
# TODO GRAPHS
color = [(254/256,97/256,0/256),(255/256,176/256,0),(100/256,143/256,255/256),(120/256,94/256,240/256)] #["red","darkorange","navy","royalblue"]
lines = ["-","--","-.",":"] #  ["solid", "dash", "dot", "dashdoted"]
figvalue=plt.figure()
for i in range(len(model_list)):
  plt.plot([int(re.sub("[^0-9]", "",instance_list[i]))],[value[j[0[0[i]]]] for j in range(len(instance_list))], lines[i], label=model_list[i].split(".")[0],  color=color[i])
plt.xlabel("Instance size")
plt.ylabel("Value")
plt.legend(loc='best')
figvalue.tight_layout()
plt.savefig(image+os.sep()+"Value_instance_K-3_P-1.png",dpi=500)
plt.close(figvalue)#plt.show() #

figvalue=plt.figure()
for i in range(len(model_list)):
  plt.plot([int(re.sub("[^0-9]", "",instance_list[i]))],[time[j[0[0[i]]]] for j in range(len(instance_list))], lines[i], label=model_list[i].split(".")[0],  color=color[i])
plt.xlabel("Instance size")
plt.ylabel("Time")
plt.legend(loc='best')
figvalue.tight_layout()
plt.savefig(image+os.sep()+"Time_instance_K-3_P-1.png",dpi=500)
plt.close(figvalue)#plt.show() #
