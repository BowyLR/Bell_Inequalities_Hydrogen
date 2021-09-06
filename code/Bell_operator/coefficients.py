# Importing python packages
import numpy as np
import itertools

# Importing Pauli matrices
from pauli_matrices.pauli_matrices import I, X, Y, Z

# Importing functions
from pauli_matrices.pauli_chain import get_Pauli_chain


# Defining list with all pauli matrices
Pauli_matrices = [I, X, Y, Z]


def get_coefficients(Bell_terms, H, N):

    r""" 
        Obtains the coefficients gamma which are used to construct the Bell operator and calculate the classical bound.
        For this we first set up the system of equations by calculating setting up a projector.
        This projector is a kronecker product of all possible conbination of the pauli matrices and the identity for N parties.
        We then calculate the system of equations via
            Tr(B*P) = Tr(H*P),
        where B are the elements in the Bell_terms list and P is the projector.
        We loop over the elements in Bell_terms and over all possible projectors
        This results in a 4^N x len(Bell_terms) matrix which we solve to find gamma with a least squared approximation.

        
        Parameters
        ----------
        Bell_terms: list
            list of all possible kronecker products of the quantum observables in M.
        H: numpy array
            2^N x 2^N numpy array
            Hamiltonian describing the system.
        N: integer
            number of parties
        
            
        Returns:
        --------
        gamma: list
            value of all the coefficients gamma
    """

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
    gamma = np.linalg.lstsq(syst_equations, syst_vector, rcond=None)[0]
    return gamma
 
