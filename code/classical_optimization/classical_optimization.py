# Importing python packages
# import numpy as np
import numba as nb

# Importing functions
# from classical_optimization.Bell_inequality import calc_Bell_inequality


@nb.jit( nopython=True ) 
def classical_optimization(coeffs, indices, possible_configurations, N, m):

    # Initializing the matrix M, 0'th row does not correspond to a measurement
    M_c = [ [ nb.types.int64(1) for _ in range(N) ] for _ in range(m+1) ] # Initial column corresponds to the identity matrices

    # Initializing Bell inequality
    I = nb.types.float64( 1e6 )

    # Looping over all possible configurations
    for conf in possible_configurations:
        
        # Updating correlation matrix
        for i in range(1, m+1):

            for j in range(N):

                M_c[i][j] = nb.types.int64( conf[j+(i-1)*N] )

        # Initializing a list containing the products of the correlation matrix
        M_list = [] 

        # Looping over the indices
        for idxs in indices:

            M_prod = nb.types.int64( 1 )

            # Calculating the product of all terms
            for k in range( len(idxs) ):

                M_prod *= M_c[idxs[k]][k]
            
            # Adding the terms to M_list
            M_list.append( M_prod )

        # Calculating the new inequality
        I_new = nb.types.float64( 0 )
        for i in range( len(coeffs) ):
            
            I_new +=  M_list[i] * nb.types.float64( coeffs[i] ) 

        # Checking if we have a new minimum
        if I_new < I:
            I = I_new

    # Returning the lowest Bell inequality
    return I



# def classical_optimization(coeffs, indices, possible_configurations, N, m):

#     # Initializing the matrix M, 0'th row does not correspond to a measurement
#     M_c = np.ones((m+1, N))  

#     # Initializing Bell inequality
#     I = 1e6

#     # Looping over all possible configurations
#     for conf in possible_configurations:
#         # Updating correlation matrix
#         M_c[1:, :] = np.reshape(conf, (m, N))

#         # Calculating the Bell inequality
#         I_new = calc_Bell_inequality(M_c, coeffs, indices)

#         # Checking if we have a new minimum
#         if I_new < I:
#             I = I_new

#     # Returning the lowest Bell inequality
#     return I



