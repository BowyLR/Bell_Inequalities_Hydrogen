# Importing python packages
# import numpy as np
from numba import jit

# Importing functions
# from classical_optimization.Bell_inequality import calc_Bell_inequality


@jit(nopython=True) 
def classical_optimization(coeffs, indices, possible_configurations, N, m):

    # Initializing the matrix M, 0'th row does not correspond to a measurement
    M_c = [ [ 1 for _ in range(N) ] for _ in range(m+1) ]  # The plus 1 corresponds to the identity matrix

    # Initializing Bell inequality
    I = 1e6

    # Looping over all possible configurations
    for conf in possible_configurations:

        # Updating correlation matrix
        for i in range(1, m+1):
            for j in range(N):
                M_c[i][j] = conf[j+(i-1)*N]

        # Initializing a list containing the products of the correlation matrix
        M_list = []

        # Looping over the indices
        for idxs in indices:

            M_prod = 1

            # Calculating the product of all terms
            for k in range(len(idxs)):
                M_prod *= M_c[idxs[k]][k]
            
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



# def classical_optimization(coeffs, indices, possible_configurations, N, m):

#     # Initializing the matrix M, 0'th row does not correspond to a measurement
#     M_c = np.ones((m+1, N))  

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



