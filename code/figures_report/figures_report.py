# Importing python functions
import numpy as np
import numba as nb

# Importing matplotlib packages
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Importing figure properties
from figures_report.figure_properties import *

# Importing functions
from hamiltonian.hamiltonian import get_hamiltonian
from Bell_operator.Bell_operator import calc_Bell_operator

# Importing Pauli matrices
from pauli_matrices.pauli_matrices import X, Z


def figure_1():
    # Defining the number of qubits
    N = 2

    # Defining the Hamiltonian
    H = np.sqrt(2)*(np.kron(Z, Z) + np.kron(X, X)) 

    # Calculating the eigenvalues and eigenstates
    eig_vals, eig_vecs = np.linalg.eigh(H)

    # extracting the ground state energy and the respective eigenstate
    psi_G = np.zeros((2**N, 1))
    psi_G[:,0] = eig_vecs[:,0]
    rho = np.matmul(psi_G, np.transpose(psi_G))

    # Defining number of measurements and the angles of the N*m-1 measurements
    m = nb.typed.List(
        [2, 2]
    )

    # Defining measurements
    theta = [np.pi/4]

    # Initializing basis and extra Z gate
    basis = 'XZ'
    extra_gate = False

    # Calculating eigenvalue of the Bell operator
    beta_Q= calc_Bell_operator(theta, H, rho, N, m, basis, extra_gate)

    # loading files
    beta_C = np.load(save_file_dir+'toy_model/beta_C.npy')
    stored_angles = np.load(save_file_dir+'toy_model/theta.npy')

    # Obtaining number of steos
    n = np.linspace(1, len(beta_C), len(beta_C))

    # Obtaining lines at which the new basinhopping steps start.
    starts = []
    for j in range(len(n)-1):
        if beta_C[j]-beta_C[j+1] > .1:
            starts.append(n[j]+1)

    # Initlazing the text locations and linewidth 
    x = -10
    lw = .75

    # initializing figures
    fig = plt.figure(figsize=(1/0.75 * lwidth, 1/0.85 * lwidth/2))
    gs = gridspec.GridSpec(
        nrows=1,
        ncols=2,
        figure=fig,
        width_ratios=[1, 1],
        height_ratios=[1],
        wspace=0.2,
    )

    # Adding first figure
    ax = fig.add_subplot(gs[0, 0])
    ax.plot( n, np.array(beta_C) , label=r'-$\beta_C$', linewidth=lw)
    ax.plot( n, np.ones(len(beta_C))*beta_Q, 
    label=r'$\mathrm{Tr}\left( \mathcal{H} \rho \right)$', linewidth=lw)
    for j in starts:
        ax.axvline(j, color='k', linestyle='--', alpha=.2)
    ax.legend()
    ax.set_xlabel(r'Iteration number')
    ax.set_ylim((-3.5, -1.5))
    ax.set_xlim((0, n[-1]))
    ax.set_yticks([-3.5, -3, -2.5, -2, -1.5])
    ax.set_yticklabels([-3.5, -3, -2.5, -2, -1.5])
    ax.text(x, -1.375, '$\mathbf{a}$',  size = MEDIUM_SIZE)

    # Adding second figure
    ax = fig.add_subplot(gs[0, 1])
    angles = np.array(stored_angles)
    ax.plot(n, angles, label=r'$\theta$', linewidth=lw)
    for j in starts:
        ax.axvline(j, color='k', linestyle='--', alpha=.2)
    ax.set_yticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
    ax.set_yticklabels([0, r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$'])
    ax.set_xlabel(r'Iteration number')
    ax.set_ylim((0, np.pi))
    ax.set_xlim((0, n[-1]))
    ax.set_ylabel(r'$\theta$')
    ax.text(x, np.pi+.175, '$\mathbf{b}$',  size = MEDIUM_SIZE)

    # Saving the figure
    fig.savefig(save_fig_dir+'toy_model/results_toy_model.png', dpi = 300, bbox_inches='tight')   




def figure_2():

    # Set parameters to make a simple molecule.
    geometry = [('H', (0., 0., 0.)), ('H', (0., 0., .7414))] # in Angstroms
    basis = 'sto-3G'
    multiplicity = 1
    charge = 0

    # Defining qubit transform
    qubit_transform = 'JW'

    # Obtainig Hamiltonian
    hamiltonian = get_hamiltonian(geometry, basis, multiplicity, charge, qubit_transform = qubit_transform)

    # Extracting parameters
    H = hamiltonian.matrix_form
    N = hamiltonian.N

    # Calculating the eigenvalues and eigenstates
    eig_vals, eig_vecs = np.linalg.eigh(H)

    # Extracting the ground state energy and the respective eigenstate
    psi_G = np.zeros((2**N, 1))
    psi_G[:,0] = eig_vecs[:,0]
    rho = np.matmul(psi_G, np.transpose(psi_G))


    # Defining small error
    eps = 1e-3

    # Defining number of measurements and the angle between the two parties
    m = nb.typed.List(
        [2,2,2,2]
    )
    theta = np.random.rand(N-1)*(np.pi-2*eps)+eps

    # # Initializing basis and extra Z gate
    basis = 'XY'
    extra_gate = True

    # Calculating eigenvalue of the Bell operator
    beta_Q = calc_Bell_operator(theta, H, rho, N, m, basis, extra_gate)

    # Initlazing the text locations and linewidth
    x_2 = -12
    x_1 = -17
    lw = .75
    y_1 = -1.785
    # scale = 1.3
    scale = 1

    # Loading Jordan-Wigner data
    beta_C_JW = np.load(save_file_dir+'hydrogen_mol/beta_C_JW.npy')
    angles_JW = np.load(save_file_dir+'hydrogen_mol/theta_JW.npy')
    n_JW = np.linspace(1, len(beta_C_JW), len(beta_C_JW))


    # Loading Bravyi-Kitaev data
    beta_C_BK = np.load(save_file_dir+'hydrogen_mol/beta_C_BK.npy')
    angles_BK = np.load(save_file_dir+'hydrogen_mol/theta_BK.npy')
    n_BK = np.linspace(1, len(beta_C_BK), len(beta_C_BK))

    # initializing figures
    fig = plt.figure(figsize=(1/0.75 * lwidth/scale, 1/0.85 * lwidth/scale))
    gs = gridspec.GridSpec(
        nrows=2,
        ncols=2,
        figure=fig,
        width_ratios=[1, 1],
        height_ratios=[1,1],
        wspace = 0.12,
        hspace = 0.12
    )

    # Plotting minimal classical bound of the JW transform
    ax = fig.add_subplot(gs[0, 0])
    ax.plot( n_JW, np.array(beta_C_JW) , label=r'$\beta_C$', linewidth=lw)
    ax.plot( n_JW, np.ones(len(beta_C_JW))*beta_Q,
        label=r'$\mathrm{Tr}\left( \mathcal{H} \rho \right)$', linewidth=lw
    )
    ax.set_xticklabels([])
    ax.set_ylim((-2, -1.8))
    ax.set_xlim((0, n_JW[-1]))
    lines = [70, 103, 163, 210]
    for j in lines:
        ax.axvline(j, color='k', linestyle='--', alpha=.2)
    ax.set_yticks([-2, -1.95, -1.9, -1.85, -1.8])
    ax.set_yticklabels([-2, -1.95, -1.9, -1.85, -1.8])
    ax.text(x_1, y_1, '$\mathbf{a}$',  size = MEDIUM_SIZE)

    # Plotting optimal angles of the JW transform
    ax = fig.add_subplot(gs[1, 0])
    for j in range(3):
        ax.plot(n_JW, np.array(angles_JW)[:,j], label=r'$\theta^{(' +str(j)+')}$', linewidth=lw)
    lines = [70, 103, 163, 210]
    for j in lines:
        ax.axvline(j, color='k', linestyle='--', alpha=.2)
    ax.set_yticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
    ax.set_yticklabels([0, r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$'])
    ax.set_xlabel(r'Iteration number')
    ax.set_ylim((0, np.pi))
    ax.set_xlim((0, n_JW[-1]))
    ax.set_ylabel(r'$\theta^{(i)}$')
    ax.text(x_1, np.pi+.1, '$\mathbf{c}$',  size = MEDIUM_SIZE)

    # Plotting the classical bound of the BK transform
    ax = fig.add_subplot(gs[0, 1])
    ax.plot( n_BK, np.array(beta_C_BK) , label=r'$-\beta_C$', linewidth=lw)
    ax.plot( n_BK, np.ones(len(beta_C_BK))*beta_Q,
        label=r'$\mathrm{Tr}\left( \mathcal{H} \rho \right)$', linewidth=lw
    )
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_ylim((-2, -1.8))
    ax.set_xlim((0, n_BK[-1]))
    ax.legend()
    lines = [43, 78, 121, 165]
    for j in lines:
        ax.axvline(j, color='k', linestyle='--', alpha=.2)
    ax.set_yticks([-2, -1.95, -1.9, -1.85, -1.8])
    ax.set_yticklabels([])
    ax.text(x_2, y_1, '$\mathbf{b}$',  size = MEDIUM_SIZE)

    # Removing 1 from number of steps (forgot to save initial value)
    n_BK = np.linspace(1, len(beta_C_BK)-1, len(beta_C_BK)-1)

    # Plotting angles of the BK transform
    ax = fig.add_subplot(gs[1, 1])
    for j in range(3):
        ax.plot(n_BK, np.array(angles_BK)[:,j], label=r'$\theta^{(' +str(j+1)+')}$', linewidth=lw)
    lines = [43, 78, 121, 165]
    for j in lines:
        ax.axvline(j, color='k', linestyle='--', alpha=.2)
    ax.set_yticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
    ax.set_yticklabels([0, r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$'])
    ax.set_xlabel(r'Iteration number')
    ax.set_ylim((0, np.pi))
    ax.set_xlim((0, n_BK[-1]))
    ax.legend(loc='lower right')
    ax.set_yticklabels([])
    ax.text(x_2, np.pi+.1, '$\mathbf{d}$',  size = MEDIUM_SIZE)

    # Saving the figure
    fig.savefig(save_fig_dir+'hydrogen_mol/results_hydrogen.png', dpi = 300, bbox_inches='tight')   





def figure_3():
    # Set parameters to make a simple molecule.
    geometry = [('H', (0., 0., 0.)), ('H', (0., 0., .7414))] # in Angstroms
    basis = 'sto-3G'
    multiplicity = 1
    charge = 0

    # Defining qubit transform
    qubit_transform = 'JW'

    # Obtainig Hamiltonian
    hamiltonian = get_hamiltonian(geometry, basis, multiplicity, charge, qubit_transform = qubit_transform)

    # Extracting parameters
    H = hamiltonian.matrix_form
    N = hamiltonian.N

    # Calculating the eigenvalues and eigenstates
    eig_vals, eig_vecs = np.linalg.eigh(H)

    # Extracting the ground state energy and the respective eigenstate
    psi_G = np.zeros((2**N, 1))
    psi_G[:,0] = eig_vecs[:,0]
    rho = np.matmul(psi_G, np.transpose(psi_G))


    # Defining small error
    eps = 1e-3

    # Defining number of measurements and the angle between the two parties
    m = nb.typed.List(
        [2,2,2,2]
    )
    theta = np.random.rand(N-1)*(np.pi-2*eps)+eps

    # # Initializing basis and extra Z gate
    basis = 'XY'
    extra_gate = True

    # Calculating eigenvalue of the Bell operator
    beta_Q = calc_Bell_operator(theta, H, rho, N, m, basis, extra_gate)

    # Initlazing the text locations and linewidth
    x = -10
    y = -1.83
    lw = .75

    # initializing figures
    fig =plt.figure(figsize=(1/0.75 * lwidth, 1/0.85 * lwidth/2))
    gs = gridspec.GridSpec(
        nrows=1,
        ncols=2,
        figure=fig,
        width_ratios=[1, 1],
        height_ratios=[1],
        wspace=0.2,
    )

    # Loading Jordan-Wigner data
    beta_C = np.load(save_file_dir+'hydrogen_mol/beta_C_JW.npy')
    beta_C1 = np.load(save_file_dir+'hydrogen_mol/beta_C_JW_m_9.npy')
    beta_C2 = np.load(save_file_dir+'hydrogen_mol/beta_C_JW_m_10.npy')
    beta_C3 = np.load(save_file_dir+'hydrogen_mol/beta_C_JW_m_11.npy')

    # Plotting the Jordan-Wigner data
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(np.linspace(1, len(beta_C), len(beta_C)) , beta_C , label=r'$n_{m=3} = 0$', linewidth=lw)
    ax.plot(np.linspace(1, len(beta_C1), len(beta_C1)) , beta_C1 , label=r'$n_{m=3} = 1$', linewidth=lw)
    ax.plot(np.linspace(1, len(beta_C2), len(beta_C2)) , beta_C2 , label=r'$n_{m=3} = 2$', linewidth=lw)
    ax.plot(np.linspace(1, len(beta_C3), len(beta_C3)) , beta_C3 , label=r'$n_{m=3} = 3$', linewidth=lw)
    ax.plot(
        np.linspace(1, len(beta_C3), len(beta_C3)), np.ones(len(beta_C3))*beta_Q, 'k',
        label=r'$\mathrm{Tr}\left( \mathcal{H} \rho \right)$', linewidth=lw
    )
    ax.set_xlabel(r'Iteration number')
    ax.set_xlim((1, 200))
    ax.set_ylim((-2.02, -1.84))
    # ax.set_ylabel(r'$-\beta_C$')
    ax.text(x, y, '$\mathbf{a}$',  size = MEDIUM_SIZE)

    # Loading Bravyi-Kitaev data
    beta_C = np.load(save_file_dir+'hydrogen_mol/beta_C_BK.npy')
    beta_C1 = np.load(save_file_dir+'hydrogen_mol/beta_C_BK_m_9.npy')
    beta_C2 = np.load(save_file_dir+'hydrogen_mol/beta_C_BK_m_10.npy')
    beta_C3 = np.load(save_file_dir+'hydrogen_mol/beta_C_BK_m_11.npy')

    # Plotting the Bravyi-Kitaev data
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(np.linspace(1, len(beta_C), len(beta_C)) , beta_C , label=r'$-\beta_C(n_m = 0)$', linewidth=lw)
    ax.plot(np.linspace(1, len(beta_C1), len(beta_C1)) , beta_C1 , label=r'$-\beta_C(n_m = 1)$', linewidth=lw)
    ax.plot(np.linspace(1, len(beta_C2), len(beta_C2)) , beta_C2 , label=r'$-\beta_C(n_m = 2)$', linewidth=lw)
    ax.plot(np.linspace(1, len(beta_C3), len(beta_C3)) , beta_C3 , label=r'$-\beta_C(n_m = 3)$', linewidth=lw)
    ax.plot(
        np.linspace(1, len(beta_C3), len(beta_C3)), np.ones(len(beta_C3))*beta_Q, 'k',
        label=r'$\mathrm{Tr}\left( \mathcal{H} \rho \right)$', linewidth=lw
    )
    ax.legend(loc='best', bbox_to_anchor=(0.5,0.45))
    ax.set_xlabel(r'Iteration number')
    ax.set_xlim((1, 200)) 
    ax.set_ylim((-2.02, -1.84))
    ax.text(x, y, '$\mathbf{b}$',  size = MEDIUM_SIZE)

    # Saving the figure
    fig.savefig(save_fig_dir+'hydrogen_mol/results_hydrogen_varying_m.png', dpi = 300, bbox_inches='tight')   