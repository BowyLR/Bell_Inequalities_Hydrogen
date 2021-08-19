# Importing python packages
import numpy as np

# Importing functions
from pauli_matrices.pauli_chain import get_Pauli_chain

# Importing Pauli matrices
from pauli_matrices.pauli_matrices import I, X, Y, Z


# Initializing used lists containing the Pauli matrices
pauli_matrices = [I, X, Y, Z]
pauli_strings = ['I', 'X', 'Y', 'Z']


def get_matrix_form(N, qubit_hamiltonian):
    # Initializing the hamiltonian
    H = np.zeros((2**N, 2**N))

    # Calculating the Hamiltonian
    for key in qubit_hamiltonian.terms.keys():
        # Initializing the indices
        idx = [0]*N

        # Updating the possible indices
        for i in range(len(key)):
            for j in range(len(pauli_strings)):
                if key[i][1] == pauli_strings[j]:
                    idx[key[i][0]] = j

        # Defining the Pauli chain
        matrices = [pauli_matrices[k] for k in idx]
        p_chain = get_Pauli_chain(matrices)

        # Adding terms to the Hamiltonian
        H += np.real( qubit_hamiltonian.terms[key] * np.array(p_chain) )

    return H