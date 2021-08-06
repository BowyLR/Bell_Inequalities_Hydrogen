# Importing python funcitons
import numpy as np
from sympy import symbols

# Importing Pauli matrices
from pauli_matrices.pauli_matrices import I, X, Y, Z


def get_correlation_matrix(measurement_angles, N, m, model='hydrogen'):

    # Initializing the correlation matrix
    M = [ [ [] for _ in range(m+1) ] for _ in range(N) ]

    # Computing its elements
    for i in range(N):
        for j in range(m+1):
            if j == 0:
                M[i][j] = I
            else:
                M[i][j] = np.cos(measurement_angles[i, j-1])*X + np.sin(measurement_angles[i, j-1])*Z

    # Converting the correlation matrix to a numpy array and returing it
    return np.array(M)
    

def get_correlation_matrix_symbolic(measurement_angles, N, m, model='hydrogen'):
    # Initializing the variables
    variables = {}

    # Obtaining variables in front of the coefficients
    for i in range(N):
        for j in range(m): 

            # Obtaining the variable indices
            lower_string = str(i) + str(j)

            # Finilizing the cosine variable
            var_name = 'x_{' + lower_string + '}'
            variables[var_name] = symbols(var_name)

            # Finilizing the sine variable
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