# -*- coding: utf-8 -*-
import math

import pyomo.environ as pyo
import numpy as np
import time
import sys
import os

global n, K, P, k_global, C
global lambda1, lambda2, theta

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
    global k_global

    cost = 0
    for j in range(k_global,n-1):
        for i in range(j+1, n):
            cost += model.C[i,j] * model.Y[i,j]
            
    C9 = [model.Z[i] for i in range(n)]
    C13 = model.Z[k_global]
    
    return cost - lambda2 * C13 - sum([a*b for a,b in zip(lambda1,C9)])

def Constraint_9(model,i):
    """ all nodes must belong to one and only one subset"""
    return sum(model.Z[i,k] for k in range(i+1)) == 1

def Constraint_10(model,i,j):
    """ Link y and z variables"""
    if k_global <= min(i, j) and j<i:
        return model.Y[i,j] <= model.Z[i] 
    return pyo.Constraint.Skip

def Constraint_11(model,i,j):
    """ Link y and z variables"""
    if k_global <= min(i, j) and j<i:
        return model.Y[i,j] <= model.Z[j]    
    return pyo.Constraint.Skip

def Constraint_12(model):
    """ Make sure a node is representative if and only if zkk is equal to one, and that each node represents at most P − 1 other nodes, hence leading to subsets withat most P nodes."""

    return sum(model.Z[i] for i in range(k_global+1,n)) <= (P-1) * model.Z[k_global]
    
def Constraint_13(model):
    """ Each subset must be used (K in the partition)."""
    return sum(model.Z[k,k] for k in range(n)) == K

def dic_initialize_links():
    y = {}
    for i in range(n):
        for j in range(n):
            y[i,j] = 0
    return y

def dic_initialize_subsets():
    z = {} # need dictionary for pyomo
    # z = np.empty(shape=(n,n))
    for i in range(n):
        if k_global == i%K:
            z[i] = 1
        else:
            z[i] = 0
    return z

def fix_constraints(Z):
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

    return output_Y(Z)

def output_Y(Z):
    # Apply Constraints 10 and 11 to Y
    Y = np.zeros(shape=(n,n,n), dtype=int)

    for i in range(n):
        for j in range(i):
            for k in range(j+1):
                Y[i, j, k] = min(Z[i,k], Z[j,k])
    return Y

def compute_C9(Z):
    res = np.zeros(shape=n, dtype=int)

    for i in range(n):
        for k in range(i+1):
            res[i] += Z[i,k]

    return res-1

def compute_C13(Z):
    return sum(Z[k,k] for k in range(n)) - K

def update_lambdas(lower_bound, upper_bound, Z):
    global lambda1, lambda2, theta
    C9 = compute_C9(Z)
    C13 = compute_C13(Z)

    if sum([C9[i]**2 for i in range(n)]) == 0:
        t1 = 0
    else:
        t1 = theta*(upper_bound-lower_bound)/sum([C9[i]**2 for i in range(n)])

    if C13 == 0:
        t2 = 0
    else:
        t2 = theta*(upper_bound-lower_bound)/(C13**2)

    for i in range(n):
        lambda1[i] += t1*C9[i]
    lambda2 += t2*C13


def solve_node(C_dict, debug=False):
    ### PYOMO MODEL ###

    model = pyo.ConcreteModel()
    model.i = pyo.RangeSet(0, n - 1)
    model.j = pyo.RangeSet(0, n - 1)

    model.C = pyo.Param(model.i, model.j, initialize=C_dict)

    Z_dic = dic_initialize_subsets()
    model.Z = pyo.Var(model.i, domain=pyo.Binary, initialize=Z_dic)
    Y_dic = dic_initialize_links()
    model.Y = pyo.Var(model.i, model.j, domain=pyo.Binary, initialize=Y_dic)

    model.goal = pyo.Objective(expr=calculate_cost(model), sense=pyo.maximize)

    model.Constraint_10 = pyo.Constraint(model.i, model.j, rule=Constraint_10)
    model.Constraint_11 = pyo.Constraint(model.i, model.j, rule=Constraint_11)
    model.Constraint_12 = pyo.Constraint(rule=Constraint_12)

    if debug:
        model.display('before_solving_lagrangian_not_node.txt')

    ### SOLVING ###
    opt = pyo.SolverFactory('glpk')  # GLPK OPTION
    # opt.options['tmlim'] = 60  # LIMITATION COMPUTATION TIME
    results = opt.solve(model, tee=False)
    # print(results.solver.goal)

    if results.solver.termination_condition == pyo.TerminationCondition.optimal:
        print(model.solutions.load_from(results))

    if debug:
        model.display('after_solving_lagrangian_not_node.txt')

    return pyo.value(model.goal), pyo.value(model.Z[:]), np.array(pyo.value(model.Y[:,:])).reshape(n,n).tolist()

def find_feasible_solution(Z):
    global C

    Y = fix_constraints(Z)

    cost = 0
    for i in range(n):
        for j in range(i):
            for k in range(j + 1):
                cost += C[i, j] * Y[i, j, k]
    return cost

def solve_lagrangian(p,instance_name, debug=False):
    global n, K, P, k_global, C
    global lambda1, lambda2, theta

    ### INITIALIZE TIME MEASURE ###
    start = time.time()  # the variable that holds the starting time the code for the clock come from https://stackoverflow.com/questions/13893287/python-time-limit
    elapsed = 0  # the variable that holds the number of seconds elapsed.

    ### INITIALIZE PARAMETERS ###
    C = read_instance(instance_name)

    n = len(C) # Number of nodes
    P = int(math.ceil(n / K) * p) # Maximum nodes by subset

    # Cost matrix
    C_dict = {}
    for i in range(n):
        for j in range(n):
            C_dict[(i,j)] = C[i][j]

    # Lambdas for relaxation
    lambda1 = [1 for _ in range(n)]
    lambda1_init = lambda1.copy()
    lambda2 = 1
    lambda2_init = lambda2
    # Theta
    theta = 0.5

    # Sure lower and upper bounds
    lower_bound, upper_bound = 0, np.sum(C)
    best_lower_bound, best_upper_bound = 0, np.sum(C)
    #init_lower_bound, init_upper_bound = 0, np.sum(C)

    n_it = 0  # nb of iterations # n_init is used to optimize the theta
    min_it = 5  # nb of iterations after the algo is checking the divergence

    Z = np.ndarray(shape=(n, n))
    Z_init = Z.copy()
    Y = np.ndarray(shape=(n, n, n))
    Y_init = Y.copy()

    ### COMPUTE COST BY NODE ###
    while elapsed < 3600 and lower_bound / upper_bound < 0.999:

        lambda1_buffer, lambda2_buffer = lambda1.copy(), lambda2

        ## UPPER BOUND
        upper_bound = 0
        for k_global in range(n):
            upper_bound_buffer, Z[:,k_global], Y[:,:,k_global] = solve_node(C_dict, debug=debug)
            upper_bound += upper_bound_buffer
        upper_bound += sum(lambda1) + lambda2 * K

        ## LOWER BOUND
        Z_low = Z.copy()
        lower_bound = find_feasible_solution(Z_low) # Z reshaped to be a numpy array
        update_lambdas(lower_bound, upper_bound, Z)

        if n_it == 1:
            init_upper_bound = upper_bound

        if debug:
            print("elapsed time = %f"%elapsed)
            print("lambda1 = " + str(lambda1))
            print("lambda2 = " + str(lambda2))
            print("lower_bound = " + str(lower_bound))
            print("upper_bound = " + str(upper_bound))

            if lambda1_buffer == lambda1 and lambda2_buffer == lambda2:
                print("stalling...", end="")

                if theta >= 0.001:
                    theta *= 0.7
                    print("change theta value to %f"%theta)
                else:
                    print("end with theta value of %f"%theta)
                    break
        else:
            if lambda1_buffer == lambda1 and lambda2_buffer == lambda2:
                if theta >= 0.001:
                    theta *= 0.7
                else:
                    break
        
        if lower_bound > best_lower_bound:
            best_lower_bound = lower_bound
        if upper_bound < best_upper_bound:
            best_upper_bound = upper_bound

        if debug:
            if n_it > min_it and upper_bound > init_upper_bound:
                print("-----------------")
                print("Divergence detected")
                Z = Z_init.copy()
                Y = Y_init.copy()
                lambda1 = lambda1_init.copy()
                lambda2 = lambda2_init
                n_it = 0
                theta *= 0.5
                print("New theta: ", theta)
                print("-----------------")
        else:
            if n_it > min_it and upper_bound > init_upper_bound:
                Z = Z_init.copy()
                Y = Y_init.copy()
                lambda1 = lambda1_init.copy()
                lambda2 = lambda2_init
                n_it = 0
                theta *= 0.5

        elapsed = time.time() - start  # update the time elapsed
        n_it += 1

    ### OUTPUT RESULT FILE ###
    print("Saving file")

    # Create output folder if it does not exist
    if not os.path.exists("result"):
        os.mkdir("result")

    # Saving
    with open("result"+os.sep+"result_by_node_"+file_name.split('.')[0]+"_"+str(K)+"_"+str(P)+".txt",'w') as f:
        elapsed = time.time() - start
        f.write(str(best_lower_bound)+" "+ str(best_upper_bound)+" "+str(elapsed))


if __name__ == "__main__":
    global K
    
    # file_name = "a280.tsp"
    # file_name = "eil51.tsp"
    # file_name = "custom.tsp"
    
    file_name = sys.argv[1]
    K = int(sys.argv[2])
    p = float(sys.argv[3])
    
    solve_lagrangian(p,file_name, debug=False)

