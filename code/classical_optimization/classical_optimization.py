# Importing python packages
import numba as nb


@nb.jit( nopython=True ) 
def classical_optimization(gamma, indices, possible_configurations, N, m):

    r""" 
        Calculates the minimal classical bound by looping over all local deterministic strategies.

        
        Parameters
        ----------
        gamma: list
            value of all the coefficients gamma
        indices: list  
            indices used to index the matrix M.
            each element of indices is of length N and contains the indices for each party to acces the i'th quantum observable in M.
            calculated with itertools.product, see its documentation for more information
        possible_configurations: list
            all possible local deterministic strategies 
        N: integer
            number of parties
        m: list
            number of measurements per party.
            the length of m should be the same as the number of parties for the function to work
        
          
        Returns:
        --------
        I: float
            minimal I 
    """

    # Initializing the matrix M, 0'th row does not correspond to a measurement
    M_c = [ [ nb.types.int64(1) for _ in range(m[i]+1) ] for i in range(N) ] # Initial column corresponds to the identity matrices

    # Initializing Bell inequality
    I = nb.types.float64( 1e6 )

    # Looping over all possible configurations
    for conf in possible_configurations:
        
        # Updating correlation matrix
        counter = 0
        for i in range(N):

            for j in range(1, m[i]+1):   # Starting at 1 excludes the identity matrix

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

        # initializing the new inequality
        I_new = nb.types.float64(0)

        # Calculating the new inequality
        for i in range( len(gamma) ):
            
            I_new +=  M_list[i] * nb.types.float64( gamma[i] ) 

        # Checking if we have a new minimum
        if I_new < I:
            I = I_new

    # Returning the lowest Bell inequality
    return I



