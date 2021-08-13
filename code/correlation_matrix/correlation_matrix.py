# Importing python funcitons
import numpy as np

# Importing Pauli matrices
from pauli_matrices.pauli_matrices import I, X, Y, Z

# Defining a vector containing the pauli matrices and their respective string
pauli_matrices = [I, X, Y, Z]
pauli_string = ['I', 'X', 'Y', 'Z']


def get_measurement_angles(theta, N, m):

    # Initializing the M matrix
    measurement_angles = np.zeros((N, m))

    # Adding the angle between the different parties
    measurement_angles[1:,:] = np.transpose([theta]*m)

    # Adding angles of the measurements of a respective party
    # for j in range(m):
    #     measurement_angles[:,j] += np.pi/2/(m-1)*(j)
    
    for j in range(m):
        measurement_angles[:,j] += (j+1)*np.pi/m

    return measurement_angles


def get_correlation_matrix(theta, N, m, basis='XY', extra_Z_gate=True):

    # Defining the measurement basis
    measurement_basis = []
    for k in basis:
        boolean = [pauli_string[l] == k for l in range(len(pauli_string))]
        measurement_basis.append(pauli_matrices[np.where(boolean)[0][0]])

    # Computing the measurement angles
    measurement_angles = get_measurement_angles(theta, N, m)
    
    # Initializing the correlation matrix
    M = [ [ [] for _ in range(m+1) ] for _ in range(N) ]

    # Computing its elements
    for i in range(N):

        for j in range(m+1):
            
            if j == 0:
                M[i][j] = I

            else:
                # Calculating the correlation matrix
                M[i][j] = np.cos(measurement_angles[i, j-1])*measurement_basis[0] + np.sin(measurement_angles[i, j-1])*measurement_basis[1]

    # Adding the extra m+1 measurement Z-gate if wanted
    if extra_Z_gate:
        for j in range(N):
            M[j][:] += [Z]

    # Converting the correlation matrix to a numpy array and returing it
    return np.array(M)
    
