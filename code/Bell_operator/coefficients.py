# Importing python packages
import numpy as np
import itertools
# from sympy import solve, Eq

# Importing Pauli matrices
from pauli_matrices.pauli_matrices import I, X, Y, Z

# Importing functions
from pauli_matrices.pauli_chain import get_Pauli_chain


# Defining list with all pauli matrices
Pauli_matrices = [I, X, Y, Z]


def get_coefficients(Bell_terms, H, N):

    # Obtaining the indices for all possible the projector matrices
    indices = list(itertools.product([ i for i in range(len(Pauli_matrices)) ], repeat=N))

    # Initializing an empty array and vector that are used to obtain the system of equations
    syst_equations = np.zeros((len(indices), len(Bell_terms)))
    syst_vector = np.zeros(len(indices))

    for i in range(len(indices)):
        
        # Obtaining the projector matrix
        projector_choices = [ Pauli_matrices[indices[i][p]] for p in range(N) ]
        projector = get_Pauli_chain(projector_choices)

        # Calculating coefficients of the Hamiltonian
        syst_vector[i] = np.real(np.trace( np.matmul(H, projector) ))

        # Setting up the system of equations
        for j in range(len(Bell_terms)):
            syst_equations[i,j] = np.real(np.trace( np.matmul(Bell_terms[j], projector) ))

    # Calculating the coefficients with the least squared difference method and returning it
    return np.linalg.lstsq(syst_equations, syst_vector, rcond=None)[0]
 
