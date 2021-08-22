# Importing python packages
import sys
import numpy as np
from openfermion.transforms import jordan_wigner, bravyi_kitaev

# Importing pauli matrices
from pauli_matrices.pauli_matrices import I, X, Y, Z

# Importing functions
from hamiltonian.fermionic_hamiltonian import get_fermionic_hamiltonian
from hamiltonian.matrix_hamiltonian import get_matrix_form

class get_hamiltonian():

    def __init__(self, qubit_transform = 'JW'):

        # Implementing fail safe for unwanted models
        possible_qubit_transforms = ['JW', 'BK']
        if np.sum([(qubit_transform == j) for j in possible_qubit_transforms]) != 1:

            sys.exit('retrieved '+qubit_transform+' but this name is incompatible. \nPossible qubit transformations are \'JW\' (Jordan-Wigner) and \'BK\' (Bravyi-Kitaev)')

        # Defining the number of qubits
        self.N = 4

        # Calculating the fermionic Hamiltonian
        self.fermionic_hamiltonian = get_fermionic_hamiltonian(N=self.N)

        # Obtaining the qubit Hamiltonian
        if qubit_transform == possible_qubit_transforms[0]:
            self.qubit_hamiltonian = jordan_wigner(self.fermionic_hamiltonian)

        elif qubit_transform == possible_qubit_transforms[1]:
            self.qubit_hamiltonian = bravyi_kitaev(self.fermionic_hamiltonian)

        # Obtaining the matrix form of the Hamiltonian
        self.matrix_form = get_matrix_form(self.N, self.qubit_hamiltonian)



        


