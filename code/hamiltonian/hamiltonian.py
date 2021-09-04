# Importing python packages
import numpy as np
from openfermion.transforms import jordan_wigner, bravyi_kitaev

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

        # Initializing the qmaximum and minimum indices of the qubits
        min_arg = 1e6
        max_arg = 0

        # Looping over the keys
        for key in self.qubit_hamiltonian.terms.keys():
            
            for j in key:

                if j[0] > max_arg:
                    max_arg = j[0]
                if j[0] < min_arg:
                    min_arg = j[0]

        # Obtaining the number of qubits in the Hamiltonian
        self.N = max_arg - min_arg + 1   # The plus 1 corresponds to accounting for the lowest index as well

        # Obtaining the matrix form of the Hamiltonian
        self.matrix_form = get_matrix_form(self.N, self.qubit_hamiltonian)
