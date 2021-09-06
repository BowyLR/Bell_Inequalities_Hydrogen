# Importing python functions
import itertools
import numba as nb
import numpy as np

# Importing functions
from Bell_operator.Bell_operator import get_Bell_terms
from Bell_operator.coefficients import get_coefficients
from classical_optimization.classical_optimization import classical_optimization
from quantum_observables.quantum_observables import get_quantum_observables


def calc_classical_bound(theta, H, N, m, basis, extra_gate):

    r""" 
        Calculates the classical bound
        It does so by first constructing the quantum observables with the get_quantum_observables function.
        Using the quantum obervables, we calculate all possible products between the obervables of 
        different parties using the get_Bell_terms function.
        We then solve the system of equations to obtain the coefficients gamma.
        Using these coefficients gamma, we calculate the classical bound by comparing all possible configurations
        of the local deterministic stategies.

        
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
        beta_C: float
            classical bound
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

    # Obtaining all possible configurations of the correlation matrix
    possible_configurations = nb.typed.List(itertools.product([1, -1], repeat=np.sum(m)))

    # Calculating and returning the classical bound
    return classical_optimization(coeffs, nb.typed.List(indices), possible_configurations, N, m)