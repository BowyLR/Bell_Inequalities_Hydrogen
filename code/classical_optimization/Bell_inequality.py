# Importing python packages
import numpy as np

def calc_Bell_inequality(M_c, coeffs, indices):
    # Initializing a list containing the products of the correlation matrix
    M_list = []

    # Looping over the indices
    for idxs in indices:
        M_prod = 1

        # Calculating the product of all terms
        for j in range(len(idxs)):
            M_prod *= M_c[idxs[j], j]
        
        # Adding the terms to M_list
        M_list.append( M_prod )

    # Calculating the new inequality and returning it
    return np.sum( [M_list[i] * coeffs[i] for i in range(1, len(coeffs))] )
