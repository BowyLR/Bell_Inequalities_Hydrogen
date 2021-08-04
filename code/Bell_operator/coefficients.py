# Importing python packages
import numpy as np
import itertools
from sympy import solve, Eq

# Importing Pauli matrices
from pauli_matrices.pauli_matrices import I, X, Y, Z

# Importing functions
from pauli_matrices.pauli_chain import get_Pauli_chain

# Obtaining the system of equations
pauli_matrices = [I, X, Y, Z]


def get_system_of_equations(B, H, N):

    # Initializing the equations list
    eqs = []

    # Calculating new set of indices
    indices = list(itertools.product([i for i in range(len(pauli_matrices))], repeat=N))

    # Looping over the equations
    for projector_indices in indices:

        # Constructing the measurement choices
        measurement_choices = [pauli_matrices[p] for p in projector_indices]
        Projector = get_Pauli_chain(measurement_choices)

        # calculating systems of equations
        eqs.append( 
            Eq( 
                np.trace( np.matmul(B, Projector) ), 
                np.trace( np.matmul(H, Projector) )
            ) 
        )

    return eqs


def get_coefficients(B, H, var_dict, N):
    
    # Obtaining the system of equations
    eqs = get_system_of_equations(B, H, N)
        
    # Solving the system of equations
    ans = solve(eqs, [var_dict[i] for i in var_dict.keys()])

    # Updating the dictionary
    for key in ans.keys():
        var_dict[str(key)] = ans[key]

    return var_dict