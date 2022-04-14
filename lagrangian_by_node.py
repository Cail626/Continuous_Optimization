import math

import pyomo.environ as pyo
import numpy as np

global n, K, P

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
    cost = 0
    for j in range(k,n-1):
        for i in range(j+1, n):
            cost += model.C[i,j]*model.Y[i,j]
            
    C9 = [model.Z[i] for i in range(0,n)]
    C13 = model.Z[k]
    
    return cost - lambda2*C13 - sum([a*b for a,b in zip(lambda1,C9)])

def Constraint_9(model,i):
    """ all nodes must belong to one and only one subset"""
    # return (0,sum(model.Z[i,k] for k in range(i)),1)
    return sum(model.Z[i,k] for k in range(i+1)) == 1

def Constraint_10(model,i,j):
    """ Link i and z variables"""
    if k <= min(i, j) and j<i:
        return model.Y[i,j] <= model.Z[i] 
    return pyo.Constraint.Skip

def Constraint_11(model,i,j):
    """ Link y and z variables"""
    if k <= min(i, j) and j<i:
        return model.Y[i,j] <= model.Z[j]    
    return pyo.Constraint.Skip
  
def Constraint_12(model):
    """ Make sure a node is representative if and only if zkk is equal to one, and that each node represents at most P − 1 other nodes, hence leading to subsets withat most P nodes."""
    if k == n-1:
        return pyo.Constraint.Skip

    return sum(model.Z[i] for i in range(k+1,n)) <= (P-1) * model.Z[k]
    
def Constraint_13(model):
    """ Each subset must be used (K in the partition)."""
    return sum(model.Z[k,k] for k in range(n)) == K

def dic_initialize_links():
    y = {}
    for i in range(n):
        for j in range(n):
            y[i,j] = 0
    return y

def dic_initialize_subsets(k):
    z = {} # need dictionary for pyomo
    # z = np.empty(shape=(n,n))
    for i in range(n):
        if k == i%K:
            z[i] = 1
        else:
            z[i] = 0
    return z

def solve_node(C_dict):
    
    model = pyo.ConcreteModel()
    model.i = pyo.RangeSet(0,n-1)
    model.j = pyo.RangeSet(0,n-1)

    model.C = pyo.Param(model.i,model.j,initialize=C_dict)

    Z_dic = dic_initialize_subsets(k)
    model.Z = pyo.Var(model.i, domain=pyo.Binary, initialize=Z_dic)
    Y_dic = dic_initialize_links()
    model.Y = pyo.Var(model.i,model.j,domain=pyo.Binary, initialize=Y_dic)

    model.goal = pyo.Objective(expr = calculate_cost(model), sense = pyo.maximize)

    model.Constraint_10 = pyo.Constraint(model.i,model.j,rule=Constraint_10)
    model.Constraint_11 = pyo.Constraint(model.i,model.j,rule=Constraint_11)    
    model.Constraint_12 = pyo.Constraint(rule=Constraint_12)
    
    opt = pyo.SolverFactory('glpk')
    opt.options['tmlim'] = 60

    opt.solve(model, tee=True)
    return pyo.value(model.goal)

    

def solve_lagrangian(instance_name):
    global n, K, P,k
    global lambda1, lambda2, C9, C13

    C = read_instance(instance_name)

    n = len(C)
    K = 3
    P = math.ceil(n / K)
    
    lambda1 = [1 for i in range(n)]
    lambda2 = 1

    C_dict = {}
    for i in range(n):
        for j in range(n):
            C_dict[(i,j)] = C[i][j]

    cost = 0
    for k in range(n):
        cost += solve_node(C_dict)

    cost += sum(lambda1) + lambda2*K

    print("cost final: ", cost)


if __name__ == "__main__":

    #file_name = "a280.tsp"
    file_name = "eil51.tsp"
    file_name = "custom.tsp"
    solve_lagrangian(file_name)
