# Importing python packages
from openfermion import MolecularData, get_fermion_operator, FermionOperator
from openfermionpyscf import run_pyscf

# Initializing the openfermionpyscf software
run_scf = 1
run_mp2 = 0
run_cisd = 0
run_ccsd = 0
run_fci = 0
delete_input = True
delete_output = True


def get_fermionic_hamiltonian(geometry, basis, multiplicity, charge):

    r""" 
        Calculates the Hamiltonian in terms of fermionic CAPs for a given input geometry,
        basis from which we approximate the atomic orbitals, multiplicity and charge.
        It uses the openfermioncpyscf package for this.
        For more information, see their documentation on how to use this package.
        We remove the electrostatic repulsion between the nuclei from the Hamiltonian,
        such that the Hamiltonian is in line with the literature.

        
        Parameters
        ----------
        geometry: list
            contains the atoms included in the molecule alongside with their respective cartesian coordinates
        basis: string
            basis which we use to approximate the atomic orbitals
            see https://github.com/pyscf/pyscf/tree/master/pyscf/gto/basis for all possible basis.
        multiplicity: integer
            spin multiplicity
        charge: integer
            charge of the molecule in e
        
          
        Returns:
        --------
        fermionic_hamiltonian: object
            Hamiltonian in terms of the fermionic CAPs without the electrostatic repulsion of the nuclei.
    """

    # Defining the molecule
    mol = MolecularData(geometry, basis, multiplicity,
                         charge, str(geometry[-1][-1][-1])
        )

    # Running pyscf
    molecule = run_pyscf(mol,
        run_scf=run_scf,
        run_mp2=run_mp2,
        run_cisd=run_cisd,
        run_ccsd=run_ccsd,
        run_fci=run_fci
    )

    # Obtaining the fermionic Hamiltonian
    molecular_hamiltonian = molecule.get_molecular_hamiltonian()
    fermionic_hamiltonian = get_fermion_operator(molecular_hamiltonian)

    # Removing the nuclear repuslion 
    fermionic_hamiltonian += -FermionOperator('', fermionic_hamiltonian.terms[()])

    return fermionic_hamiltonian