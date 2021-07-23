# Importing python packages
import numpy as np
from openfermion.ops import FermionOperator

def get_fermionic_hamiltonian(N=4):
    # Let us define a molecular Hamiltonian that we can use to feed it to openfermion
    h11 = h22 = -1.252477
    h33 = h44 = -0.475934
    v1221 = v2112 = 0.674493
    v3443 = v4334 = 0.697397
    v1331 = v2442 = v2332 = v2442 = v3113 = v1441 = v3223 = v4224 = 0.663472
    v1313 = v2424 = v3241 = v3421 = v1423 = v1234 = 0.181287

    # Initializing single particle and Coulomb operators
    h_single_and_coulomb = np.array([
        [h11, v1221, v1331-v1313, v1441       ],
        [0,   h22,   v2332,       v2442-v2424 ],
        [0,   0,     h33,         v3443       ],
        [0,   0,     0,           h44         ]
    ]) 

    # Initializing the Hamiltonian
    hamiltonian_fermionic_op = 0

    # Adding the Coulom interactions to the Hamiltonian
    for i in range(N):
        # Adding single Hamiltonian terms
        hamiltonian_fermionic_op += FermionOperator(str(i)+'^ '+str(i), h_single_and_coulomb[i, i])

        # Adding coulomb terms
        for j in range(i+1, N):
            hamiltonian_fermionic_op += FermionOperator(str(j) + '^ ' + str(j) + ' ' + str(i) + '^ ' + str(i), h_single_and_coulomb[i, j])

    # Adding the double excitation operators
    hamiltonian_fermionic_op += (FermionOperator('0^ 1^ 3 2', v1234) + FermionOperator('2^ 3^ 1 0', v1234) +
                FermionOperator('0^ 3^ 1 2', v1234) + FermionOperator('2^ 1^ 3 0', v1234)
    )

    return hamiltonian_fermionic_op