# Importing python funcitons
import numpy as np

# Importing function which creates angles
from correlation_matrix.measurement_angles import get_measurement_angles

# Importing Pauli matrices
from pauli_matrices.pauli_matrices import I, X, Y, Z

# Defining a vector containing the pauli matrices and their respective string
pauli_matrices = [I, X, Y, Z]
pauli_string = ['I', 'X', 'Y', 'Z']


def get_correlation_matrix(theta, N, m, basis='XY', extra_gate=True):
    
    # Returning valueerror if needed
    if (type(theta) == float or type(theta) ==  int):
        if N != 2:
            raise ValueError('Number of of angles in theta should be equal to N-1')
    
    else: 
        if len(theta) != N-1:
            raise ValueError('Number of of angles in theta should be equal to N-1')

    # Obtaining measurement angles
    angles = get_measurement_angles(theta, N, m)

    # Printing possible errors
    if len(angles) != np.sum(m):
        raise ValueError('length of the array \'theta\' is not equal to N*m-1.')

    if len(m) != N:
        raise ValueError('Number of elements in measurement array m not the same as the number of parties N.')

    # Defining the measurement basis
    measurement_basis = []
    for k in basis:
        boolean = [pauli_string[l] == k for l in range(len(pauli_string))]
        measurement_basis.append(
            pauli_matrices[np.where(boolean)[0][0]]
        )

    # Initializing empty array for the angles of the quantum observables
    measurement_angles = [ [ 0 for _ in range(np.max(m[i])) ] for i in range(N) ]

    # Intializing a counter 
    counter = 0

    # Looping over the different indixes of the matrix
    for i in range(N):

        for j in range(m[i]):

            # Adding angles to the matrix
            measurement_angles[i][j] = angles[counter]

            # Updating counter
            counter += 1

    # Initializing the correlation matrix
    M = [ [ [] for _ in range(m[i]+1) ] for i in range(N) ]

    # # Calculating the elements of the correlation matrix with quantum observables
    for i in range(N):

        for j in range(m[i]+1):
            
            if j == 0:
                M[i][j] = I

            else:
                M[i][j] = np.cos(measurement_angles[i][j-1])*measurement_basis[0] + np.sin(measurement_angles[i][j-1])*measurement_basis[1]
    
    # Adding the extra m+1 measurement Z-gate if needed
    if extra_gate:
        for i in range(1, len(pauli_string)):
            if pauli_string[i] not in basis:
                for j in range(N):
                    M[j][:] += [ pauli_matrices[i] ]

    # Converting the correlation matrix to a numpy array and returing it
    return M




