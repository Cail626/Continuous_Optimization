import subprocess
import math

for instance in ["custom.tsp"]:
    for K in ["3", "4", "5"]:
        for P in [0.5,1,1.5,2]:
            bashCmd = ["python", "lagrangian_by_node.py",instance,K,str(P)]
            process = subprocess.Popen(bashCmd)
            process.wait()
            