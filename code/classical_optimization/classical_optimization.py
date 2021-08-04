# Importing python packages
import itertools
import numpy as np

# Importing functions
from classical_optimization.Bell_inequality import calc_inequality_terms

def classical_optimization(var_dict, N, m):
    # Initializing the matrix M, 0'th row does not correspond to a measurement
    M_c = np.ones((m+1, N))  

    # Initializing temporary indices
    indices = list(itertools.product([i for i in range(m+1)], repeat=N))

    # Obtaining all possible configurations of the correlation matrix
    possible_configurations = list(itertools.product([1, -1], repeat=m*N))

    # Initializing Bell inequality
    I = 1e6

    # Looping over all possible configurations
    for conf in possible_configurations:
        # Updating correlation matrix
        M_c[1:, :] = np.reshape(conf, (m, N))

        # Calculating the Bell inequality
        I_new = calc_inequality_terms(M_c, var_dict, indices)

        # Checking if we have a new minimum
        if I_new < I:
            I = I_new

    return I




# def classical_optimization(var_dict, N, m, recursive=True):

#     # Initializing the matrix M, 0'th row does not correspond to a measurement
#     M_c = np.ones((m+1, N))   # DOes this need to be initialized in a random manner?

#     # Initializing temporary indices
#     temp_indices = list(itertools.product([i for i in range(m+1)], repeat=N))


#     # Obtaining the classical bound
#     if recursive:

#         # Initializing E_i
#         E_i = 0

#         # Calculating the possible outcomes for the i'th row of the matrix M
#         possible_configurations = list(itertools.product([1, -1], repeat=m))

#         # Calculating the i'th indices
#         for idx in tqdm(range(N)):

#             # This needs to be done for each iteration of idx, the above does not have to be done each iteration
#             indices = []

#             # First condition ensures that we do not have any double counting
#             # Second condition ensures that we do not take a constant non measurable value into account
#             # Last condition makes it so we do not capture any measurements from other observables
#             for j in temp_indices:

#                 if np.sum(j[:idx]) == 0 and np.sum(j) != 0 and j[idx] != 0: 
#                     indices.append([*[0] * (idx), *j[idx:]])

#             # Updating correlation matrix and setting a very high value for h_i such that we always start with a lowest configuration after one iteration
#             M_c_temp = M_c.copy()
#             h_i = 1e6

#             # Looping over the different configurations
#             for conf in possible_configurations:

#                 # Updating the correlation matrix
#                 M_c_temp[1:, idx] = conf

#                 # Calculating h_i
#                 h_i_temp = calc_inequality(M_c_temp, var_dict, indices)

#                 # Checking if the new h_i is smaller than the old one, if yes, then we adapt the new value
#                 if h_i_temp < h_i:
#                     M_c = M_c_temp.copy()
#                     h_i = h_i_temp

#             # Adding the last configuration of h_i to E_i
#             E_i += h_i

#         # Obtaining the classical bound
#         beta_C = -E_i

#     else:
#         # Obtaining all possible indices and configurations
#         indices = temp_indices
#         possible_configurations = list(itertools.product([1, -1], repeat=m*N))

#         # Initializing Bell inequality
#         I = 1e6

#         for conf in tqdm(possible_configurations):
#             # Updating correlation matrix
#             M_c[1:, :] = np.reshape(conf, (m, N))

#             # Calculating the Bell inequality
#             I_new = calc_inequality(M_c, var_dict, indices)

#             # Checking if we have a new minimum
#             if I_new < I:
#                 I = I_new

#         beta_C = -I

#     # Returning the classical bound
#     return beta_C



    
    