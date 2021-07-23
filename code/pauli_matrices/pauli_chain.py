# Importing python packages
import numpy as np

def get_Pauli_chain(list_w_matrices):

    # Initializing the Pauli chain
    Pauli_chain = list_w_matrices[0]

    # Adding the additional terms
    for matrix in list_w_matrices[1:]:
        Pauli_chain = np.kron(Pauli_chain, matrix)

    # returning the final Pauli chain
    return Pauli_chain