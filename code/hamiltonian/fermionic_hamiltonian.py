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

    # Defining the molecule
    description = str(charge)
    mol = MolecularData(geometry, basis, multiplicity,
                         charge, description)

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