# Importing python funcitons
import numpy as np
import itertools
from sympy import symbols

# Importing Pauli matrices
from pauli_matrices.pauli_matrices import I, X, Y, Z


def get_correlation_matrix(N, m, model='hydrogen'):
    # Defining measurement angles
    indices = list(itertools.product([i for i in range(m)], repeat=N))

    # Obtaining variables in front of the coefficients
    variables = {}
    for idx in indices:

        # Obtaining the variable indices
        lower_string = ''

        for j in idx:
            lower_string += str(j)

        # Finilizing the variable
        var_name = 'x_{' + lower_string + '}'
        variables[var_name] = symbols(var_name)

        var_name = 'y_{' + lower_string + '}'
        variables[var_name] = symbols(var_name)

    # Initializing the correlation matrix
    M = [ [ [] for _ in range(m+1) ] for _ in range(N) ]

    # Computing its elements
    for i in range(N):
        for j in range(m+1):
            if j == 0:
                M[i][j] = I
            else:
                M[i][j] = variables['x_{'+str(i)+str(j-1)+'}']*X + variables['y_{'+str(i)+str(j-1)+'}']*Z

    # Converting the correlation matrix to a numpy array and returing it
    return np.array(M), variables