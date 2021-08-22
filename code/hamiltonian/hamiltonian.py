# Importing python packages
import numpy as np
from openfermion.transforms import jordan_wigner, bravyi_kitaev

# Importing pauli matrices
from pauli_matrices.pauli_matrices import I, X, Y, Z

# Importing functions
from hamiltonian.fermionic_hamiltonian import get_fermionic_hamiltonian
from hamiltonian.matrix_hamiltonian import get_matrix_form


class get_hamiltonian():

    def __init__(self, geometry, basis, multiplicity, charge, qubit_transform='JW'):

        # Implementing fail safe for unwanted models
        possible_qubit_transforms = ['JW', 'BK']
        if np.sum([(qubit_transform == j) for j in possible_qubit_transforms]) != 1:

            raise ValueError('retrieved '+ qubit_transform +
                            ' as an input, but this name is incompatible. ' + 
                            '\nPossible qubit transformations are \'JW\' (Jordan-Wigner) and \'BK\' (Bravyi-Kitaev)')

        # Calculating the fermionic Hamiltonian
        self.fermionic_hamiltonian = get_fermionic_hamiltonian(geometry, basis, multiplicity, charge)

        # Obtaining the qubit Hamiltonian
        if qubit_transform == possible_qubit_transforms[0]:
            self.qubit_hamiltonian = jordan_wigner(self.fermionic_hamiltonian)

        elif qubit_transform == possible_qubit_transforms[1]:
            self.qubit_hamiltonian = bravyi_kitaev(self.fermionic_hamiltonian)

        # Obtaining the number of qubits in the Hamiltonian
        self.N = self.qubit_hamiltonian.many_body_order()

        # Obtaining the matrix form of the Hamiltonian
        self.matrix_form = get_matrix_form(self.N, self.qubit_hamiltonian)
