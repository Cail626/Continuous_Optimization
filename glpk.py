import pyomo.environ as pyo
import numpy as np

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
            buffer = file.readline().split()
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

def calculate_cost(C,Y,n):
    cost = 0
    for i in range(n):
        for j in range(n):
            for k in range(0,min(i,j)):
                cost += C[i,j]*Y[i,j,k]
    return cost

def Constaint9(model,i,Z):
    """ all nodes must belong to one and only one subset"""
    return sum(model.Z[i,:i]) == 1


def Constraint_10(model,i,j,k):
    """ Link y and z variables"""
    if(k <= min(i,j)):
        return model.Y[i,j,k] <= model.Z[i,k] 
    return Constraint.Skip

def Constraint_11(model,i,j,k):
    """ Link y and z variables"""
    if(k <= min(i,j)):
        return model.Y[i,j,k] <= model.Z[j,k]    
    return Constraint.Skip
  
def Constraint_12(model,k,n,P):
    """ Make sure a node is representative if and only if zkk is equal to one, and that each node represents at most P − 1 other nodes, hence leading to subsets withat most P nodes."""
    return sum(model.Z[i,k] for i in range(k,n)) <= (P-1)*model.Z[k,k]
    
def constraint13(model,K,n):
    """ Each subset must be used (K in the partition)."""
    return sum(model.Z[k,k] for k in range(n)) == K

def solve_lagrangian():
    
    instance_name = 'eil51' # 'a280.tsp'
    C = read_instance(instance_name)
    n = len(C)
    K = 10
    print(n)
   
    C_dict = {}
    for i in range(0,n):
        for j in range(0,n):
            C_dict[(i,j)] = C[i][j]
    
    model = pyo.ConcreteModel()
    model.i = pyo.RangeSet(0,n-1)
    model.j = pyo.RangeSet(0,n-1)
    model.k = pyo.RangeSet(0,n-1)

    model.C = pyo.Param(model.i,model.j,initialize=C_dict)

    model.Z = pyo.Var(model.i,model.k, domain=pyo.Binary)
    model.Y = pyo.Var(model.i,model.j,model.k,domain=pyo.Binary)

    model.goal = pyo.Objective(expr = calculate_cost(model.C,model.Y,n), sense = pyo.minimize)

    model.Constraint_9 = pyo.Constraint(model.Z,K,n,rule=Constraint_9)
    model.Constraint_10 = pyo.Constraint(model.i,model.j,model.k,rule=Constraint_10)
    model.Constraint_11 = pyo.Constraint(model.i,model.j,model.k,rule=Constraint_11)    
    model.Constraint_12 = pyo.Constraint(model.k,k,n,P,rule=Constraint_12)    
    model.Constraint_13 = pyo. constraint13(model.Z,K,n) 
    opt = pyo.SolverFactory('glpk')
    opt.solve(model) 
    print(pyo.value(model.obj))


solve_lagrangian()
