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

def fix_constraints(Z: np.ndarray) -> np.ndarray:
    """
    Fix the parameter Z for the constraints 9 and 13
    :param Z: Parameter Z of the problem. See paper.
    """

    ### Fix Constraint 9 :

    # Will check if there is many zero on the lower triangular
    # matrix and will keep only the rightmost one.
    # (explanation is following after the example)
    ##  Example :
    # [1, ., ., ., .]
    # [0, 1, ., ., .]
    # [0, 1, 1, ., .]
    # [0, 0, 1, 1, .]
    # [1, 0, 0, 0, 1]
    ##  will become :
    # [1, ., ., ., .]
    # [0, 1, ., ., .]
    # [0, 0, 1, ., .]
    # [0, 0, 0, 1, .]
    # [0, 0, 0, 0, 1]
    # If the line of the lower triangle matrix is full of zero,
    # it puts a zero on diagonal
    ## Example :
    # [0, ., ., ., .]
    # [0, 0, ., ., .]
    # [0, 0, 0, ., .]
    # [0, 0, 0, 0, .]
    # [0, 0, 0, 0, 0]
    ## will become :
    # [1, ., ., ., .]
    # [0, 1, ., ., .]
    # [0, 0, 1, ., .]
    # [0, 0, 0, 1, .]
    # [0, 0, 0, 0, 1]

    for i in range(n):
        null = True
        for k in range(i, -1, -1):
            if Z[i, k] == 1:
                null = False
                Z[i, :k] = 0
        if null:
            Z[i, i] = 1
    ### End constraint 9



    ### Fix Constraint 13 :
    # Will compare the sum of the diagonal (nb_subsets) to the
    # number of subset we need (K)
    #  If there is too many ones on the diagonal, we move the
    # ones on the diagonal on their respective line. We have to
    # not open a new subset, and we choose heuristically
    # the most filled subset (but not full).
    #  If there is too few ones on the diagonal, we add ones
    # on the diagonal.
    #  For the both case we change the matrix from the lower part
    # of the lower triangle matrix
    nb_subsets = sum(Z[k, k] for k in range(n))
    col = np.array([(P - sum(Z[i:, i])) * Z[i, i] for i in range(n)])

    if nb_subsets > K:
        for i in range(n - 1, -1, -1):
            if nb_subsets <= K: break
            if Z[i, i] == 1:
                Z[i, i] = 0
                nb_subsets -= 1

                disp_subsets = np.where(col > 0)[0]
                col_min = disp_subsets[np.argmin(col[disp_subsets])]
                Z[i, col_min] = 1
                col[col_min] -= 1

    elif nb_subsets < K:
        for i in range(n - 1, -1, -1):
            if nb_subsets >= K: break
            if Z[i, i] == 0:
                Z[i, i] = 1
                nb_subsets += 1

                Z[i, :i] = 0

    # End Constraint 13
    return Z

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
