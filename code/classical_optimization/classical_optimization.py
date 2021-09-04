# Importing python packages
import numba as nb


@nb.jit( nopython=True ) 
def classical_optimization(coeffs, indices, possible_configurations, N, m):

    # Initializing the matrix M, 0'th row does not correspond to a measurement
    M_c = [ [ nb.types.int64(1) for _ in range(m[i]+1) ] for i in range(N) ] # Initial column corresponds to the identity matrices

    # Initializing Bell inequality
    I = nb.types.float64( 1e6 )

    # Looping over all possible configurations
    for conf in possible_configurations:
        
        # Updating correlation matrix
        counter = 0
        for i in range(N):

            for j in range(1, m[i]+1):

                # Adding possible configuration
                M_c[i][j] = nb.types.int64( conf[counter] )

                # Updating the counter
                counter += 1

        # Initializing a list containing the products of the correlation matrix
        M_list = [] 

        # Looping over the indices
        for idxs in indices:

            M_prod = nb.types.int64( 1 )

            # Calculating the product of all terms
            for k in range(N):

                M_prod *= M_c[k][idxs[k]]
            
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



