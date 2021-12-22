# Bell_Inequalities_Hydrogen
We investigate if the family of Bell inequalities constructed from spin Hamiltonians, acquired through mapping the fermionic Hamiltonian describing the Hydrogen molecule the Jordan-Wigner and the Bravyi-Kitaev transformation, are equivalent or inequivalent to eachother. 
Moreover, we also investigate if the ground states of the corresponding spin Hamiltonians show nonlocal correlations.
To get a better understanding of the code and how the simulations are done, it is recommended to read the [report](./report_final_version.pdf)
The repository is devided into two folders.
## Code
The [code](./code) folder contains the code used to do the simulations:
* [toy_hamiltonian.ipynb](./code/toy_hamiltonian.ipynb): construction of Bell inequalities for a simple toy Hamiltonian taht results in the CHSH inequality.
* [molecular_hydrogen.ipynb](./code/molecular_hydrogen.ipynb): construction of Bell inequalities from the spin Hamiltonians of the Hydrogen molecule.
* [plots_report.ipynb](./code/plots_report.ipynb): contains all the code that generates the plots from the report
The rest of this folder contains other folders that contain written functions used throughout the simulations.
The respective functions all contain a detailed description of what they do, what their input and what their outputs are.


## Data
Contains the data saved from the simulations and the figures.
