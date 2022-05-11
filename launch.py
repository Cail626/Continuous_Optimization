import subprocess
import math

for instance in ["aaaaaa", "bbbbb", "ccccccc"]:
    for K in ["3", "4", "5"]:
        nk = math.floor(n/K)
        for P in [0.5,1,1.5,2]:
            bashCmd = ["python", "lagrangian_by_node.py",instance,K,P]
            process = subprocess.Popen(bashCmd)
