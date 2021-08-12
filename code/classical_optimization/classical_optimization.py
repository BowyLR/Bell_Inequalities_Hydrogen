# Importing python packages
import itertools
import numpy as np

# Importing functions
from classical_optimization.Bell_inequality import calc_Bell_inequality


def classical_optimization(coeffs, indices, N, m):

    # Initializing the matrix M, 0'th row does not correspond to a measurement
    M_c = np.ones((m+1, N))  

    # Obtaining all possible configurations of the correlation matrix
    possible_configurations = list(itertools.product([1, -1], repeat=m*N))

    # Initializing Bell inequality
    I = 1e6

    # Looping over all possible configurations
    for conf in possible_configurations:
        # Updating correlation matrix
        M_c[1:, :] = np.reshape(conf, (m, N))

        # Calculating the Bell inequality
        I_new = calc_Bell_inequality(M_c, coeffs, indices)

        # Checking if we have a new minimum
        if I_new < I:
            I = I_new

    # Returning the lowest Bell inequality
    return I



