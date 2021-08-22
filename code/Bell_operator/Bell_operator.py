# Importing python packages
from sympy import symbols
import itertools
import numba as nb

# Importing functions
from pauli_matrices.pauli_chain import get_Pauli_chain


def get_Bell_terms(M, N, m):
    # Initializing a list which contains all posible terms from the Bell operator
    Bell_terms = []

    # Creating indices for all measurements
    indices = nb.typed.List( 
        list(itertools.product([ i for i in range(m+1) ], repeat=N)) 
    )

    # Calculating the Bell operator
    for measurement_indices in indices:

        # Constructing the measurement choices  
        measurement_choices = [ M[p, measurement_indices[p]] for p in range(N) ]

        # Obtaining Pauli chain
        Bell_terms.append( get_Pauli_chain(measurement_choices) )

    # Returning the Bell terms
    return Bell_terms, indices
    