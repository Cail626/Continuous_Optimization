import math
import pyomo.environ as pyo
import numpy as np
import time
import sys
import os
global n, K, P
import matplotlib
import matplotlib.pyplot as plt
import re

matplotlib.rcParams['lines.linewidth'] = 4
matplotlib.rcParams['font.size'] = 14   # fontsize of the axes title
matplotlib.rcParams['axes.labelsize'] = 20   # fontsize of the axes title
matplotlib.rcParams['axes.titlesize'] = 20   # fontsize of the axes title


def get_results(filename):
  with open(filename,'r') as f:
        line = f.read()
        value = float(line.split(" ")[0])
        time = float(line.split(" ")[1])
        return value, time

def get_results_node(filename):
  with open(filename,'r') as f:
        line = f.read()
        value_lb = float(line.split(" ")[0])
        value_ub = float(line.split(" ")[1])
        time = float(line.split(" ")[2])
        return value_lb,value_ub, time

def make_graph_instance_size(p,k):
  fig_value_lb=plt.figure()
  for m in range(len(model_list)):
    plt.plot([int(re.sub("[^0-9]", "",instance_list[j])) for j in range(len(instance_list)) ],[value_lb[j][k][p][m] for j in range(len(instance_list))], lines[m], label=model_list[m], color=color[m], alpha = 0.7, marker ='o')
  plt.xlabel("Instance size")
  plt.yscale("log")
  plt.ylabel("Value lower bound")
  plt.legend(loc='best')
  fig_value_lb.tight_layout()
  plt.savefig("image"+os.sep+"Value_lb_instance_K-"+str(K_list[k])+"_P-"+str(P_list[p])+".png",dpi=500)
  plt.close(fig_value_lb)#plt.show() #

  fig_value_ub=plt.figure()
  for m in range(len(model_list)-2):
    plt.plot([int(re.sub("[^0-9]", "",instance_list[j])) for j in range(len(instance_list)) ],[value_ub[j][k][p][m] for j in range(len(instance_list))], lines[m+2], label=model_list[m+2], color=color[m+2], alpha = 0.7, marker ='o')
  plt.xlabel("Instance size")
  plt.ylabel("Value upper bound")
  plt.legend(loc='best')
  fig_value_lb.tight_layout()
  plt.savefig("image"+os.sep+"Value_ub_instance_K-"+str(K_list[k])+"_P-"+str(P_list[p])+".png",dpi=500)
  plt.close(fig_value_ub)#plt.show() #

  figtime=plt.figure()
  for m in range(len(model_list)):
    plt.plot([int(re.sub("[^0-9]", "",instance_list[j])) for j in range(len(instance_list)) ],[time[j][k][p][m] for j in range(len(instance_list))], lines[m], label=model_list[m].split(".")[0], color=color[m], alpha = 0.7, marker ='o')
  plt.xlabel("Instance size")
  plt.ylabel("Time")
  plt.legend(loc='best')
  figtime.tight_layout()
  plt.savefig("image"+os.sep+"Time_instance_K-"+str(K_list[k])+"_P-"+str(P_list[p])+".png",dpi=500)
  plt.close(figtime)#plt.show() #

def make_graph_K(p,instance):
  fig_value_lb=plt.figure()
  for m in range(len(model_list)):
    plt.plot([int(K_list[k]) for k in range(len(K_list)) ],[value_lb[instance][k][p][m] for k in range(len(K_list))], lines[m], label=model_list[m], color=color[m], alpha = 0.7, marker ='o')
  plt.xlabel("Instance size")
  plt.yscale("log")
  plt.ylabel("Value lower bound")
  plt.legend(loc='best')
  fig_value_lb.tight_layout()
  plt.savefig("image"+os.sep+"Value_lb_K_instance-"+str(instance_list[instance])+"_P-"+str(P_list[p])+".png",dpi=500)
  plt.close(fig_value_lb)#plt.show() #

  fig_value_ub=plt.figure()
  for m in range(len(model_list)-2):
    plt.plot([int(K_list[k]) for k in range(len(K_list)) ],[value_ub[instance][k][p][m] for k in range(len(K_list))], lines[m+2], label=model_list[m+2], color=color[m+2], alpha = 0.7, marker ='o')
  plt.xlabel("Instance size")
  plt.ylabel("Value upper bound")
  plt.legend(loc='best')
  fig_value_lb.tight_layout()
  plt.savefig("image"+os.sep+"Value_ub_K_instance-"+str(instance_list[instance])+"_P-"+str(P_list[p])+".png",dpi=500)
  plt.close(fig_value_ub)#plt.show() #

  figtime=plt.figure()
  for m in range(len(model_list)):
    plt.plot([int(K_list[k]) for k in range(len(K_list)) ],[time[instance][k][p][m] for k in range(len(K_list))], lines[m], label=model_list[m].split(".")[0], color=color[m], alpha = 0.7, marker ='o')
  plt.xlabel("Instance size")
  plt.ylabel("Time")
  plt.legend(loc='best')
  figtime.tight_layout()
  plt.savefig("image"+os.sep+"Time_K_instance-"+str(instance_list[instance])+"_P-"+str(P_list[p])+".png",dpi=500)
  plt.close(figtime)#plt.show() #


def make_graph_P(k,instance):
  fig_value_lb=plt.figure()
  for m in range(len(model_list)):
    plt.plot([float(P_list[p]) for p in range(len(P_list)) ],[value_lb[instance][k][p][m] for p in range(len(P_list))], lines[m], label=model_list[m], color=color[m], alpha = 0.7, marker ='o')
  plt.xlabel("Instance size")
  plt.yscale("log")
  plt.ylabel("Value lower bound")
  plt.legend(loc='best')
  fig_value_lb.tight_layout()
  plt.savefig("image"+os.sep+"Value_lb_P_instance-"+str(instance_list[instance])+"_K-"+str(K_list[k])+".png",dpi=500)
  plt.close(fig_value_lb)#plt.show() #

  fig_value_ub=plt.figure()
  for m in range(len(model_list)-2):
    plt.plot([float(P_list[p]) for p in range(len(P_list)) ],[value_ub[instance][k][p][m] for p in range(len(P_list))], lines[m+2], label=model_list[m+2], color=color[m+2], alpha = 0.7, marker ='o')
  plt.xlabel("Instance size")
  plt.ylabel("Value upper bound")
  plt.legend(loc='best')
  fig_value_lb.tight_layout()
  plt.savefig("image"+os.sep+"Value_ub_P_instance-"+str(instance_list[instance])+"_K-"+str(K_list[k])+".png",dpi=500)
  plt.close(fig_value_ub)#plt.show() #

  figtime=plt.figure()
  for m in range(len(model_list)):
    plt.plot([float(P_list[p]) for p in range(len(P_list)) ],[time[instance][k][p][m] for p in range(len(P_list))], lines[m], label=model_list[m].split(".")[0], color=color[m], alpha = 0.7, marker ='o')
  plt.xlabel("Instance size")
  plt.ylabel("Time")
  plt.legend(loc='best')
  figtime.tight_layout()
  plt.savefig("image"+os.sep+"Time_P_instance-"+str(instance_list[instance])+"_K-"+str(K_list[k])+".png",dpi=500)
  plt.close(figtime)#plt.show() #



model_list =  ["model_1", "model_2", "result_not_node","result_by_node"] # ["model_1", "model_2", "result_not_node"]  # 
instance_list = ["eil51.tsp", "eil76.tsp",  "eil101.tsp", "bier127.tsp", "a280.tsp"] # , "pr439.tsp"]
K_list = ["3", "4", "5"] #["3"] # 
P_list =  ["1.0", "1.5", "2.0", "2.5", "3.0"]
value_lb = [ [ [ [ [] for m in range(len(model_list)) ] for k in range(len(P_list))] for j in range(len(K_list)) ] for i in range(len(instance_list))]
value_ub = [ [ [ [ [] for m in range(len(model_list)-2) ] for k in range(len(P_list))] for j in range(len(K_list)) ] for i in range(len(instance_list))]
time = [ [ [ [ [] for m in range(len(model_list)) ] for k in range(len(P_list))] for j in range(len(K_list)) ] for i in range(len(instance_list))]

for instance in range(len(instance_list)):
    for K in range(len(K_list)):
        for P in range(len(P_list)):
            for model in range(len(model_list)):
                filename = "result3"+os.sep+model_list[model]+"_"+instance_list[instance].split('.')[0]+"_"+str(K_list[K])+"_"+str(P_list[P])+".txt"
                if "node" not in model_list[model]:
                  try :
                    value_lb[instance][K][P][model], time[instance][K][P][model] = get_results(filename)
                  except FileNotFoundError:
                    value_lb[instance][K][P][model], time[instance][K][P][model] = 0, 0
                else:
                  try :
                    value_lb[instance][K][P][model], value_ub[instance][K][P][model-2], time[instance][K][P][model] = get_results_node(filename)
                  except FileNotFoundError:
                    value_lb[instance][K][P][model], value_ub[instance][K][P][model-2], time[instance][K][P][model] = 0, 0,0


# MAKE GRAPHS
color = [(254/256,97/256,0/256),(255/256,176/256,0),(100/256,143/256,255/256),(120/256,94/256,240/256)] #["red","darkorange","navy","royalblue"]
lines = ["-","--","-.",":"] #  ["solid", "dash", "dot", "dashdoted"]
if not os.path.exists("image"):
  os.mkdir("image")
make_graph_instance_size(0,0)
make_graph_instance_size(4,0)
make_graph_K(0,0)
make_graph_P(0,0)
