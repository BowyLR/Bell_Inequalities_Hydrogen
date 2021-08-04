# Importing python packages
from sympy import symbols
import itertools

# Importing functions
from pauli_matrices.pauli_chain import get_Pauli_chain


def get_Bell_term(measurement_choices, measurement_indices):
            
    # Obtaining Pauli chain
    Pauli_chain = get_Pauli_chain(measurement_choices)

    # Obtaining the variable indices
    lower_string = ''

    for j in measurement_indices:
        lower_string += str(j)

    # Finilizing the variable
    var_name = 'c_{' + lower_string + '}'
    var = symbols(var_name)

    # Returning the Bell operator term
    return var * Pauli_chain, var_name, var


def get_Bell_operator(M, m, N):
    # Initializing the Bell operator, the Bell inequality and a string containing the variables
    B = 0
    var_dict = {}

    # Creating indices for all measurements
    indices = list(itertools.product([ i for i in range(m+1) ], repeat=N))

    # Calculating the Bell operator
    for measurement_indices in indices:

        # Constructing the measurement choices
        measurement_choices = [ M[p, measurement_indices[p]] for p in range(N) ]

        # Calculating the Bell terms and the corresponding variable
        B_temp, var_name, var_temp = get_Bell_term(measurement_choices, measurement_indices)
        
        # Adding and storing the values
        B += B_temp
        var_dict[var_name] = var_temp

    # Returning the Bell operator in its matrix form and a dictionary containing the coefficients
    return B, var_dict