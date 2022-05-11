import math

import pyomo.environ as pyo
import numpy as np
import time

global n, K, P, lambda1, lambda2, C9, C13

def read_instance(file_name):
    # Read the instance and compute the instance in dictionary
    #   Pyomo needs dictionary to compute
    cost_matrix = -1
    try:
        file = open("Instances/{}".format(file_name), 'r')

        name = file.readline().split(":")
        comment = file.readline().split(":")
        tsp = file.readline().split(":")
        dimension = int(file.readline().split(":")[-1])
        edge_weight_type = file.readline().split(":")
        info = file.readline().split(":")

        node_positions = np.empty(dimension, dtype='2i')
        for line in range(dimension):
            buffer = ' '.join(file.readline().split())  # remove multiple white space (the .tsp file layout are not consistant)
            buffer = buffer.split()
            node_positions[line] = (int(buffer[1]), int(buffer[2]))

        file.close()
        cost_matrix = np.empty([dimension, dimension], dtype=float)

        x,y = 0,1
        for node1 in range(dimension):
            for node2 in range(dimension):
                cost_matrix[node1][node2] = np.sqrt(np.power(node_positions[node2][x] - node_positions[node1][x], 2) + \
                            np.power(node_positions[node2][y] - node_positions[node1][y], 2) )
    except:
        print("Error reading file.")
    return cost_matrix

def calculate_cost(model):
    global lambda1, lambda2, C9, C13, n
    cost = 0
    for i in range(n):
        C9[i] = (sum(model.Z[i,k] for k in range(i+1))-1)
        for j in range(i):
            for k in range(min(i,j)+1):
                cost += model.C[i,j]*model.Y[i,j,k]
    C13 = (sum(model.Z[k,k] for k in range(n)) - K) 
    return cost - lambda2*C13 - sum([a*b for a,b in zip(lambda1,C9)])

def Constraint_10(model,i,j,k):
    """ Link i and z variables"""
    if k <= min(i, j) and j<i:
        return model.Y[i,j,k] <= model.Z[i,k] 
    return pyo.Constraint.Skip

def Constraint_11(model,i,j,k):
    """ Link y and z variables"""
    if k <= min(i, j) and j<i:
        return model.Y[i,j,k] <= model.Z[j,k]    
    return pyo.Constraint.Skip

  
def Constraint_12(model,k):
    """ Make sure a node is representative if and only if zkk is equal to one, and that each node represents at most P − 1 other nodes, hence leading to subsets withat most P nodes."""
    if k == n-1:
        return pyo.Constraint.Skip

    return sum(model.Z[i,k] for i in range(k+1,n)) <= (P-1) * model.Z[k,k]

def dic_initialize_links():
    y = {}
    for i in range(n):
        for j in range(n):
            for k in range(n):
                y[i,j,k] = 0
    return y

def dic_initialize_subsets():
    z = {} # need dictionary for pyomo
    # z = np.empty(shape=(n,n))
    for i in range(n):
        for k in range(n):
            if k == i%K:
                z[i,k] = 1
            else:
                z[i,k] = 0
    return z

def update_lambdas(lower_bound, upper_bound):
    global lambda1, lambda2,  C9, C13

    t1 = (lower_bound-upper_bound)/sum([C9[i]**2 for i in range(n)])
    t2 = (lower_bound-upper_bound)/(C13**2)
    for i in range(n):
        lambda1[i] += t1*C9[i]
    lambda2 += t2*C13

def solve_lagrangian(instance_name):
    global n, K, P, lambda1, lambda2
    C = read_instance(instance_name)

    n = len(C)
    K = 2
    P = math.ceil(n / K)
    lambda1 = [0 for i in range(n)]
    lambda2 = 0
    lower_bound = 0
    upper_bound = 10**100
    t = None
    
    C_dict = {}
    for i in range(0,n):
        for j in range(0,n):
            C_dict[(i,j)] = C[i][j]


    model = pyo.ConcreteModel()
    model.i = pyo.RangeSet(0,n-1)
    model.j = pyo.RangeSet(0,n-1)
    model.k = pyo.RangeSet(0,n-1)

    model.C = pyo.Param(model.i,model.j,initialize=C_dict)

    Z_dic = dic_initialize_subsets()
    model.Z = pyo.Var(model.i,model.k, domain=pyo.Binary, initialize=Z_dic)
    Y_dic = dic_initialize_links()
    model.Y = pyo.Var(model.i,model.j,model.k,domain=pyo.Binary, initialize=Y_dic)

    model.goal = pyo.Objective(expr = calculate_cost, sense = pyo.maximize)

    model.Constraint_9 = pyo.Constraint(model.i,rule=Constraint_9)
    model.Constraint_10 = pyo.Constraint(model.i,model.j,model.k,rule=Constraint_10)
    model.Constraint_11 = pyo.Constraint(model.i,model.j,model.k,rule=Constraint_11)    
    model.Constraint_12 = pyo.Constraint(model.k,rule=Constraint_12)
    model.Constraint_13 = pyo.Constraint(rule=Constraint_13)

    opt = pyo.SolverFactory('glpk')
    start = time.time()         #the variable that holds the starting time the code for the clock come from https://stackoverflow.com/questions/13893287/python-time-limit
    elapsed = 0                 #the variable that holds the number of seconds elapsed.
    while elapsed < 300 and lower_bound/upper_bound > 0.999:
        update_lambdas()
        upper_bound = solve_relaxation(model) # we maximise thus the upper bound is given by the relaxation
        elapsed = time.time() - start #update the time elapsed
  
def solve_relaxation(model):

    opt.solve(model, tee=True)
    return pyo.value(model.goal)





if __name__ == "__main__":

    # file_name = "a280.tsp"
    # file_name = "eil51.tsp"
    file_name = "custom.tsp"
    solve_lagrangian(file_name)
