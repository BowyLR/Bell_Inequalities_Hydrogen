# Importing python packages
import numpy as np

def calc_inequality_terms(M, var_dict, indices):
    # Extracting the coefficients and the matrix products
    coeffs = []
    M_elements = []

    for idx in indices:
        # Initializing coefficients and matrix products
        name = ''
        M_prod = 1

        # Looping over the m measurements
        for j in range(len(idx)):
            name += str(idx[j])
            M_prod *= M[idx[j], j]

        # Storing values
        coeffs.append(var_dict['c_{'+name+'}'])
        M_elements.append(M_prod)
    
    # Calculating the Bell inequality
    return np.sum([coeffs[j]*M_elements[j] for j in range(len(indices))])
