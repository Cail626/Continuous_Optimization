import math
import pyomo.environ as pyo
import numpy as np
import sys

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
            buffer = ' '.join(file.readline().split()) # remove multiple white space (the .tsp file layout are not consistant)
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
    for i in range(n):
        for j in range(n):
            for k in range(K):
                cost += model.C[i,j]*model.Y[i,j,k]
    return cost

def Constraint_2(model,i):
    """ ensure each node is in exactly one subset"""
    # return (0,sum(model.Z[i,k] for k in range(i)),1)
    return sum(model.Z[i,k] for k in range(K)) == 1
 
def Constraint_3(model,i,j,k):
    """ Link y and z variables"""
    return model.Y[i,j,k] <= model.Z[i,k] 

def Constraint_4(model,i,j,k):
    """ Link y and z variables"""
    return model.Y[i,j,k] <= model.Z[j,k]    
  
def Constraint_5(model,k):
    """ Make sure a node is representative if and only if zkk is equal to one, and that each node represents at most P − 1 other nodes, hence leading to subsets withat most P nodes."""
    return (1,sum(model.Z[i,k] for i in range(n)),P)

def dic_initialize_subsets():
    z = {} # need dictionary for pyomo
    for i in range(n):
        for k in range(K):
            if (i+k+1)%K == 0:
                z[i,k] = 1
            else:
                z[i,k] = 0
    return z

def dic_initialize_links():
    y = {}
    for i in range(n):
        for j in range(n):
            for k in range(K):
                y[i,j,k] = 0
    return y

    
def test_Constraint_2(Z_dic):
    for i in range(n):
        if not sum(Z_dic[i,k] for k in range(K)) == 1:
            print("Initial values violate constaint 2 at node "+str(i))

def test_Constraint_5(Z_dic):
    for k in range(K):
        if not(sum(Z_dic[i,k] for i in range(n)) >=1 and sum(Z_dic[i,k] for i in range(n)) <= P):
            print("Initial values violate constaint 5 ins subset "+str(k))

def solve_lagrangian(instance_name):
    global n, K, P

    C = read_instance(instance_name)
    
    #print("C", C)
    
    n = len(C)
    K = 5
    P = math.ceil(n / K)
   
    C_dict = {}
    for i in range(0,n):
        for j in range(0,n):
            C_dict[(i,j)] = C[i][j]


    model = pyo.ConcreteModel()
    model.i = pyo.RangeSet(0,n-1)
    model.j = pyo.RangeSet(0,n-1)
    model.k = pyo.RangeSet(0,K-1)

    model.C = pyo.Param(model.i,model.j,initialize=C_dict)
    Z_dic = dic_initialize_subsets()
    test_Constraint_2(Z_dic)
    test_Constraint_5(Z_dic)
    model.Z = pyo.Var(model.i,model.k, domain=pyo.Binary,initialize=Z_dic)
    Y_dic = dic_initialize_links()
    model.Y = pyo.Var(model.i,model.j,model.k,domain=pyo.Binary,initialize=Y_dic)

    cost = sum(model.C[i,j]*model.Y[i,j,k] for i in range(n) for j in range(n) for k in range(K))

    model.goal = pyo.Objective(expr = cost, sense = pyo.maximize)

    model.Constraint_2 = pyo.Constraint(model.i,rule=Constraint_2)
    model.Constraint_3 = pyo.Constraint(model.i,model.j,model.k,rule=Constraint_3)
    model.Constraint_4 = pyo.Constraint(model.i,model.j,model.k,rule=Constraint_4)    
    model.Constraint_5 = pyo.Constraint(model.k,rule=Constraint_5)
    opt = pyo.SolverFactory('glpk')
    opt.options['tmlim'] = 60
    model.display('test.txt')

    opt.solve(model, tee=True)
    print(pyo.value(model.goal))




if __name__ == "__main__":
    # matrix = read_instance("a280.tsp")
    # matrix = read_instance("eil51.tsp")
    # print(matrix)
    file_name = "eil51.tsp"
    solve_lagrangian(file_name)
