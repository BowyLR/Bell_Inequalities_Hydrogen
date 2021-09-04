# Importing python packages
import itertools


# NOT USED YET
def get_indices(N, m, num_qubits):

    # Initializing indices
    indices = []

    # Looping over the possible configurations
    for j in range(N - num_qubits + 1):

        # Calculating the possible confiurations
        temp_indices = list(itertools.product([ i for i in range(m+1) ], repeat=num_qubits))

        # Adding the terms to a new array
        for i in range( (m+1)**num_qubits ):

            if j == 0 or [*temp_indices[i]] != [0]*num_qubits:
                # Adding terms
                indices.append( 
                    [*[0]*j, *temp_indices[i], *[0]*(N-num_qubits-j)]
                )

    # Obtaining indices of the identity terms
    identity_indices = [0, 1, *[ j*(m+1)**num_qubits for j in range(1, N - num_qubits + 3)]]
    additions = [0, 0, *[-1*j  for j in range(N-num_qubits+2)]]
    identity_indices = [identity_indices[j] + additions[j] for j in range(N - num_qubits + 3)]

    # Seperating the indices to be usable for the classical optimzation
    opt_indices = []
    for j in range(len(identity_indices)-1):
        opt_indices.append(
            indices[identity_indices[j] : identity_indices[j+1]]
        )

    return indices, opt_indices

