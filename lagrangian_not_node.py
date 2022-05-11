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
    C9 = np.empty(shape=n, dtype=int)
    for i in range(n):
        C9[i] = (sum(pyo.value(model.Z[i,k]) for k in range(i+1))-1)
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
    # if k == n-1:
    #     return pyo.Constraint.Skip

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

    for i in range(n):
        for k in range(n):
            if k == i%K:
                z[i,k] = 1
            else:
                z[i,k] = 0
    return z

def find_feasible_solution(model):
    Z = np.array(pyo.value(model.Z[:,:])).reshape((n,n))
    Y = fix_constraints(Z)

    cost = 0
    for i in range(n):
        for j in range(i):
            for k in range(j + 1):
                cost += model.C[i, j] * Y[i, j, k]
    return cost

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
                nb_subsets -= 1

                col[i] = 0
                # all subset moved out
                for it in range(i, n):
                    if Z[it, i] == 1:
                        Z[it, i] = 0
                        disp_subsets = np.where(col > 0)[0]
                        col_min = disp_subsets[np.argmin(col[disp_subsets])]

                        Z[it, col_min] = 1
                        col[col_min] -= 1


    elif nb_subsets < K:
        for i in range(n - 1, -1, -1):
            if nb_subsets >= K: break
            if Z[i, i] == 0:
                Z[i, i] = 1
                nb_subsets += 1

                Z[i, :i] = 0

    # End Constraint 13

    return output_Y

def output_Y(Z: np.ndarray) -> np.ndarray:
    # Apply Constraints 10 and 11 to Y
    Y = np.zeros(shape=(n,n,n), dtype=int)

    for i in range(n):
        for j in range(i):
            for k in range(j+1):
                Y[i, j, k] = min(Z[i,k], Z[j,k])
    return Y

def update_lambdas(lower_bound, upper_bound):
    global lambda1, lambda2,  C9, C13

    t1 = (lower_bound-upper_bound)/sum([C9[i]**2 for i in range(n)])
    t2 = (lower_bound-upper_bound)/(C13**2)
    for i in range(n):
        lambda1[i] += t1*C9[i]
    lambda2 += t2*C13

def solve_lagrangian(instance_name, debug=False, verbose=True):
    global n, K, P
    global lambda1, lambda2, C9, C13
    
    ### INITIALIZE PARAMETERS ###
    C = read_instance(instance_name)

    n = len(C)  # Number of nodes
    K = 3  # Number of subsets
    P = math.ceil(n / K)  # Maximum nodes by subset

    # Cost matrix
    #TODO: Directly have dict from read_instance function
    C_dict = {} # Cost matrix, pyomo needs dictionary
    for i in range(0,n):
        for j in range(0,n):
            C_dict[(i,j)] = C[i][j]

    # Lambda relaxation
    lambda1 = [1 for _ in range(n)]
    lambda2 = 1

    ### PYOMO MODEL ###

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

    model.Constraint_10 = pyo.Constraint(model.i,model.j,model.k,rule=Constraint_10)
    model.Constraint_11 = pyo.Constraint(model.i,model.j,model.k,rule=Constraint_11)    
    model.Constraint_12 = pyo.Constraint(model.k,rule=Constraint_12)

    ### SOLVING ###
    opt = pyo.SolverFactory('glpk')  # GLPK OPTION
    #opt.options['tmlim'] = iteration_time  # LIMITATION COMPUTATION TIME

    if debug:
        model.display('before_solving_lagrangian_not_node.txt')

    start = time.time()  # the variable that holds the starting time the code for the clock come from https://stackoverflow.com/questions/13893287/python-time-limit
    elapsed = 0  # the variable that holds the number of seconds elapsed.

    while elapsed < 300 and lower_bound / upper_bound > 0.999:
        upper_bound = solve_relaxation(model)  # we maximise thus the upper bound is given by the relaxation
        elapsed = time.time() - start  # update the time elapsed
        lower_bound = find_feasible_solution(model)
        update_lambdas(lower_bound, upper_bound)

    if verbose:
        print(pyo.value(model.goal))

    if debug:
        model.display('after_solving_lagrangian_not_node.txt')



if __name__ == "__main__":

    #file_name = "a280.tsp"
    #file_name = "eil51.tsp"
    file_name = "custom.tsp"
    solve_lagrangian(file_name)
