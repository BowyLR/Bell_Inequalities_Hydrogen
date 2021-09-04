# Importing python packages
import itertools
import numpy as np
import numba as nb

# Importing functions
from pauli_matrices.pauli_chain import get_Pauli_chain
from Bell_operator.coefficients import get_coefficients
from correlation_matrix.correlation_matrix import get_correlation_matrix


def get_Bell_terms(M, N, m):

    # Initializing a list which contains all posible terms from the Bell operator
    Bell_terms = []

    # Creating indices for all measurements
    indices = nb.typed.List( 
        list(itertools.product(
            *[[i for i in range(m[i] + 1)] for i in range(N)]) 
        )
    )

    # Calculating the Bell operator
    for measurement_indices in indices:

        # Constructing the measurement choices  
        measurement_choices = [ M[p][measurement_indices[p]] for p in range(N) ]

        # Obtaining Pauli chain
        Bell_terms.append( get_Pauli_chain(measurement_choices) )

    # Returning the Bell terms
    return Bell_terms, indices
    



def calc_Bell_operator(theta, H, rho, N, m, basis, extra_gate):

    # Obtaining the general correlation matrix
    M = get_correlation_matrix(theta, N, m, basis=basis, extra_gate=True)

    # Adding one to number of measurements to account for the extra z_gate
    if extra_gate:
        m = nb.typed.List(
            [j + 1 for j in m]
        )

    # Calcuating a list containing all the possible matrices from kronecker products of the observables in the correlation matrix
    Bell_terms, indices = get_Bell_terms(M, N, m)

    # Calculating the coefficients
    coeffs = get_coefficients(Bell_terms, H, N)

    # Initializing the Bell operator
    B = np.zeros((2**N, 2**N), dtype='complex128')

    # Adding terms to the Bell operator
    for j in range(len(coeffs)):
        B += coeffs[j]*Bell_terms[j]

    # Calculting the quantum value
    beta_Q = np.trace(np.matmul(B, rho))
    
    # Discarding the imaginary part (its always zero) and returning it
    return np.real(beta_Q)
