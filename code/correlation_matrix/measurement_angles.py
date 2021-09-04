# importing python functions
import numpy as np

def get_measurement_angles(theta, N, m):

    # Adding 0 to the angles theta
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
