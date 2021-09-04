# Importing python functions
import itertools
import numba as nb
import numpy as np

# Importing functions
from Bell_operator.Bell_operator import get_Bell_terms
from Bell_operator.coefficients import get_coefficients
from classical_optimization.classical_optimization import classical_optimization
from correlation_matrix.correlation_matrix import get_correlation_matrix


def calc_classical_bound(theta, H, N, m, basis, extra_gate):

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

    # Obtaining all possible configurations of the correlation matrix
    possible_configurations = nb.typed.List(itertools.product([1, -1], repeat=np.sum(m)))

    # Calculating and returning the classical bound
    return classical_optimization(coeffs, nb.typed.List(indices), possible_configurations, N, m)