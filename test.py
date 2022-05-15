import numpy as np
import math
global n,K,P

def test_all_constraints(Z):
    for k in range(K):
        if not(sum(Z[i,k] for i in range(k+1,n)) <= (P-1) * Z[k,k]):
            print("Error constraint 12 and column k={k}: Z[:,{k}]={Z}".format(k=k,Z=Z[:,k]))
            # print(" sum(Z[i,k] for i in range(k+1,n)) !<= (P-1) * Z[k,k]\n {x1} !<= {x2}".format(x1=sum(Z[i,k] for i in range(k+1,n)),x2=(P-1) * Z[k,k]))
    for i in range(n):
        if not (sum(Z[i, k] for k in range(i + 1)) == 1):
            print("Error constraint 9")
            # print("> on node {i}:\n sum(Z[i,k] for k in range(i+1)) = {x} ≠ 1".format(i=i, x=sum(
            #     Z[i, k] for k in range(i + 1))))
    if not (sum(Z[k, k] for k in range(n)) == K):
        print("Error constraint 13")
        # print("> sum(Z[k,k] for k in range(n)) ≠ K \n {x1} ≠ {x2}".format(x1=sum(Z[k, k] for k in range(n)), x2=K))
    print("Test constraint end")

def test_constraint_12(Z):
    for k in range(K):
        if not(sum(Z[i,k] for i in range(k+1,n)) <= (P-1) * Z[k,k]):
            print("Error constraint 12 and column k={k}: Z[:,{k}]={Z}".format(k=k,Z=Z[:,k]))
            # print(" sum(Z[i,k] for i in range(k+1,n)) !<= (P-1) * Z[k,k]\n {x1} !<= {x2}".format(x1=sum(Z[i,k] for i in range(k+1,n)),x2=(P-1) * Z[k,k]))
    print("Test constraint 12 end")


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

    print("Z = \n"+str(Z))
    return output_Y(Z)

def output_Y(Z: np.ndarray) -> np.ndarray:
    # Apply Constraints 10 and 11 to Y
    Y = np.zeros(shape=(n,n,n), dtype=int)

    for i in range(n):
        for j in range(i):
            for k in range(j+1):
                Y[i, j, k] = min(Z[i,k], Z[j,k])
    return Y


if __name__ == "__main__":
    global n, K, P

    n = 5
    K = 2
    P = math.ceil(n / K)

    Z = [
        [0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1],
    ]

    Z = np.array(Z)

    test_constraint_12(Z)
    Y = fix_constraints(Z)
    # for k in range(n):
    #     print("\nY[:,:,"+str(k)+")] = \n"+str(Y[:,:,k]))

    test_all_constraints(Z)




