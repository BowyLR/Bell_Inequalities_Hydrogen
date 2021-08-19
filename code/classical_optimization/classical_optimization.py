# Importing python packages
# import itertools
# import numpy as np
from numba import jit

# Importing functions
# from classical_optimization.Bell_inequality import calc_Bell_inequality


@jit(nopython=True) 
def classical_optimization(coeffs, indices, possible_configurations, N, m):

    # Initializing the matrix M, 0'th row does not correspond to a measurement
    M_c = [ [ 1 for _ in range(m+1) ] for _ in range(N) ]

    # Initializing Bell inequality
    I = 1e6

    # Looping over all possible configurations
    for conf in possible_configurations:

        # Updating correlation matrix
        for i in range(N):
            for j in range(1, m+1):
                M_c[i][j] = conf[j-1+i*(N-1)]

        # Calculating the inequality
        # Initializing a list containing the products of the correlation matrix
        M_list = []

        # Looping over the indices
        for idxs in indices:

            M_prod = 1

            # Calculating the product of all terms
            for j in range(N):
                M_prod *= M_c[j][idxs[j]]
            
            # Adding the terms to M_list
            M_list.append( M_prod )

        # Calculating the new inequality
        I_new = 0
        for i in range( len(coeffs) ):
            I_new += M_list[i] * coeffs[i]

        # Checking if we have a new minimum
        if I_new < I:
            I = I_new

    # Returning the lowest Bell inequality
    return I




# def classical_optimization(coeffs, indices, N, m):

#     # Initializing the matrix M, 0'th row does not correspond to a measurement
#     M_c = np.ones((m+1, N))  

#     # Obtaining all possible configurations of the correlation matrix
#     possible_configurations = list(itertools.product([1, -1], repeat=m*N))

#     # Initializing Bell inequality
#     I = 1e6

#     # Looping over all possible configurations
#     for conf in possible_configurations:
#         # Updating correlation matrix
#         M_c[1:, :] = np.reshape(conf, (m, N))

#         # Calculating the Bell inequality
#         I_new = calc_Bell_inequality(M_c, coeffs, indices)

#         # Checking if we have a new minimum
#         if I_new < I:
#             I = I_new

#     # Returning the lowest Bell inequality
#     return I



