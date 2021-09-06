# importing python functions
import numpy as np


def get_measurement_angles(theta, N, m):

    r""" 
        Calculates the angles of the first measurements of each party relative to the first measurement of the first party.
        The angles per party are as follows:
            theta = j*pi/m+theta,
        where 0<=j<m is the jth measurement of a given party and theta is an angle relative to the first party.

        
        Parameters
        ----------
        theta: list
            angle between the first measurement of each party relative to the first measurement of the zeroth party.
            should be of the length of N-1.
        N: integer
            number of parties
        m: list
            number of measurements per party.
            the length of m should be the same as the number of parties for the function to work
        
          
        Returns:
        --------
        angles: list
            angles of all quantum observables relative to 0
    """

    # Adding 0 to the initial angle
    theta = [0, *theta]

    # Initializing empty angles array
    angles = []

    # Creating angles
    for i in range(N):

        for j in range(m[i]):

            # Adding one additional angle to the list
            angles.append(
                j*np.pi/m[i] + theta[i]
            )
    
    # returning the angles
    return angles
