# Importing python packages
import itertools
import numpy as np
import numba as nb

# Importing functions
from pauli_matrices.pauli_chain import get_Pauli_chain
from Bell_operator.coefficients import get_coefficients
from quantum_observables.quantum_observables import get_quantum_observables


def get_Bell_terms(M, N, m):

    r""" 
        Obtains all possible kronecker products of the quantum observables of each party.
        These quantum observables contain an identity matrix for each party to account for lower body correlators.
        It calculates all these kronecker products by first calculating the indices to 
        acces a matrix containing all quantum observables and calculating the Kronecker product between the indiced matrices.
        We calculate these indices with pythons itertools package.
        It stores each matrix product in a list called Bell_terms.

        
        Parameters
        ----------
        M: list
            matrix containing the quantum observables of each party.
            zeroth element is the identity matrix
        N: integer
            number of parties
        m: list
            number of measurements per party.
            the length of m should be the same as the number of parties for the function to work
        
            
        Returns:
        --------
        Bell_terms: list
            list of all possible kronecker products of the quantum observables in M.
        indices: list  
            indices used to index the matrix M.
            each element of indices is of length N and contains the indices for each party to acces the i'th quantum observable in M.
            calculated with itertools.product, see its documentation for more information
    """

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

    r""" 
        Calculates the ground state of the Bell operator with the classical bound included.
        It does so by first constructing the quantum observables with the get_quantum_observables function.
        Using the quantum obervables, we calculate all possible products between the obervables of 
        different parties using the get_Bell_terms function.
        We then solve the system of equations to obtain the coefficients gamma, 
        which we then use to construct the Bell operator.
        This function essentially returns the ground state of the Hamiltonian H as of now.
        Still WIP to be able to change the Bell operator.

        
        Parameters
        ----------
        theta: list
            angle between the first measurement of each party relative to the first measurement of the zeroth party.
            should be of the length of N-1.
        H: numpy array
            2^N x 2^N numpy array
            Hamiltonian describing the system.
        rho: numpy array
            2^N x 2^N numpy array
            density operator of the ground state of the system
        N: integer
            number of parties
        m: list
            number of measurements per party.
            the length of m should be the same as the number of parties for the function to work
        basis: string
            basis in which we perform our measurements
            should be of length 2 and one can only choose from the X, Y or Z basis.
            e.g. 'XY' is a valid option for basis, but 'Z' is not
        extra_gate: boolean
            determines if we should add an extra gate to the measurements which is not present in the basis.
            e.g. is basis is 'XY' and extra_gate is True, we add an extra measurement in the Z direction.
        
            
        Returns:
        --------
        beta_Q: float
            ground state of the Bell operator without the classical bound included
    """

    # Obtaining the general correlation matrix
    M = get_quantum_observables(theta, N, m, basis=basis, extra_gate=extra_gate)

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
