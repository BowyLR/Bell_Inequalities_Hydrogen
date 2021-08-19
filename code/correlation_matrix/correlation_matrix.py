# Importing python funcitons
import numpy as np

# Importing Pauli matrices
from pauli_matrices.pauli_matrices import I, X, Y, Z

# Defining a vector containing the pauli matrices and their respective string
pauli_matrices = [I, X, Y, Z]
pauli_string = ['I', 'X', 'Y', 'Z']


# def get_measurement_angles(theta, N, m):

#     # Initializing the M matrix
#     measurement_angles = np.zeros((N, m))

#     # Adding the angle between the different parties
#     measurement_angles[1:,:] = np.transpose([theta]*m)

#     # Adding angles of the measurements of a respective party
#     # for j in range(m):
#     #     measurement_angles[:,j] += np.pi/2/(m-1)*(j)
    
#     for j in range(m):
#         measurement_angles[:,j] += (j+1)*np.pi/m

#     return measurement_angles


def get_correlation_matrix(theta, N, m, basis='XY', extra_Z_gate=True):

    # Returning error in the case where the length of theta does not fit the number of quantum observables
    if type(theta) == float or type(theta) == int:
        if N*m-1 != 1:
            raise ValueError('length of the array \'theta\' is not equal to N*m-1.')

    elif len(theta) != N*m-1:
        raise ValueError('length of the array \'theta\' is not equal to N*m-1.')

    # Defining the measurement basis
    measurement_basis = []
    for k in basis:
        boolean = [pauli_string[l] == k for l in range(len(pauli_string))]
        measurement_basis.append(pauli_matrices[np.where(boolean)[0][0]])

    # Initializing empty array for the angles of the quantum observables
    measurement_angles = np.zeros(m*N)

    # Adding angles to the array and reshaping it. The 0'th element is always fixed
    measurement_angles[1:] = theta
    measurement_angles = np.array( np.reshape( measurement_angles, (N, m) ) )
    
    # Initializing the correlation matrix
    M = [ [ [] for _ in range(m+1) ] for _ in range(N) ]

    # # Calculating the elements of the correlation matrix with quantum observables
    for i in range(N):

        for j in range(m+1):
            
            if j == 0:
                M[i][j] = I

            else:
                M[i][j] = np.cos(measurement_angles[i, j-1])*measurement_basis[0] + np.sin(measurement_angles[i, j-1])*measurement_basis[1]

    # Adding the extra m+1 measurement Z-gate if needed
    if extra_Z_gate:
        for j in range(N):
            M[j][:] += [Z]

    # Converting the correlation matrix to a numpy array and returing it
    return np.array(M)
    
