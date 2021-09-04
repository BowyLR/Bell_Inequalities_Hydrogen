# Importing python packages
import numpy as np

# Importing functions
from pauli_matrices.pauli_chain import get_Pauli_chain

# Importing Pauli matrices
from pauli_matrices.pauli_matrices import I, X, Y, Z


# Initializing used lists containing the Pauli matrices
pauli_matrices = [I, X, Y, Z]
pauli_strings = ['I', 'X', 'Y', 'Z']


def get_matrix_form(N, qubit_hamiltonian, num_qubits=[]):
    
    # Initializing the Hamiltonian
    H = np.zeros( (2**N, 2**N) )

    # Checking if the number of qubits is not given, if so, make it so we capture the full Hamiltonian
    if num_qubits == []:
        num_qubits = N

    # Looping over the keys of the Hamiltonian
    for key in qubit_hamiltonian.terms.keys():

        # Initializing the indices
        idx = [0] * N

        # Calculating the index distance between the first and last qubit of the key
        if len(key) == 0:
            dist = 0   # Evades indexing error
        else:
            dist = key[-1][0] - key[0][0]

        # Creating a criteria that make sit so we only have num_qubits neirest neighbouring qubits
        if dist < num_qubits:

            # Updating the possible indices
            for i in range(len(key)):

                for j in range(len(pauli_strings)):

                    if key[i][1] == pauli_strings[j]:
                        idx[ key[i][0] ] = j
            
            # Defining the Pauli chain
            p_chain = get_Pauli_chain( [ pauli_matrices[k] for k in idx ] )

            # Adding terms to the Hamiltonian
            H += np.real( qubit_hamiltonian.terms[key] * np.array(p_chain) )

    # Returning the Hamiltonian
    return H