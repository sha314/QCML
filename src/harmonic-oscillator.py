import numpy as np


# X operator
LL=8  # in the limit [-L,L] with unit step
def QFT(j, k, n_qubits):
    NN = 2**n_qubits
    return np.exp(1j*2*np.pi*j*k/NN)/np.sqrt(NN)

def IQFT(j, k, n_qubits):
    NN = 2**n_qubits
    return np.exp(-1j*2*np.pi*j*k/NN)/np.sqrt(NN)


def get_Fmat(n_qubits):
    NN = 2**n_qubits
    diag = np.arange(-NN/2, NN/2, 1)
    return np.array([[QFT(j, k, n_qubits) for j in diag] for k in diag])


def get_IFmat(n_qubits):
    NN = 2**n_qubits
    diag = np.arange(-NN/2, NN/2, 1)
    return np.array([[IQFT(j, k, n_qubits) for j in diag] for k in diag])


def getHamiltonian(n_qubits):
    Xmat = get_X_operator(n_qubits)

    Fmat = get_Fmat(n_qubits)
    IFmat = get_IFmat(n_qubits)
   
    XmatSq = Xmat**2
    PmatSq = np.dot(np.dot(IFmat, XmatSq), Fmat)
    HamiltonianP = (XmatSq + PmatSq)/2 
    return HamiltonianP


def get_X_operator(n_qubits):
    """
    n_qubits : number of qubits
    """
    NN = 2**n_qubits
    diag = np.arange(-NN/2, NN/2, 1)
    # diag = diag.reshape((-1,1))
    # print(np.identity(NN))
    Xmat = np.identity(NN)* diag
    # print(Xmat)
    return Xmat/2



def get_Unitary_Ux(LL, t):
    X_hat = get_X_operator(LL)
    Ux = np.zeros(X_hat.shape)
    for i in range(X_hat.shape[0]):
        Ux[i,i] = np.exp(-1j * t * X_hat[i,i]**2)
        pass

    # Ux = np.exp(-1j*t*X_hat**2)
    # return Ux/Ux[0,0]
    return Ux

import numpy as np
import scipy as sp
import itertools
import functools as ft
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import re 



def decompose_ham_to_pauli(H):
    """Decomposes a Hermitian matrix into a linear combination of Pauli operators.

    Args:
        H (array[complex]): a Hermitian matrix of dimension 2**n x 2**n.

    Returns:
        tuple[list[float], list[string], list [ndarray]]: a list of coefficients, 
        a list of corresponding string representation of tensor products of Pauli 
        observables that decompose the Hamiltonian, and
        a list of their matrix representation as numpy arrays
    """

    n = int(np.log2(len(H)))
    N = 2 ** n

    # Sanity Checks
    if H.shape != (N, N):
        raise ValueError(
            "The Hamiltonian should have shape (2**n, 2**n), for any qubit number n>=1"
        )

    if not np.allclose(H, H.conj().T):
        raise ValueError("The Hamiltonian is not Hermitian")

    sI = np.eye(2, 2, dtype=complex)
    sX = np.array([[0, 1], [1, 0]], dtype=complex)
    sZ = np.array([[1, 0], [0,-1]], dtype=complex)
    sY = complex(0,-1)*np.matmul(sZ,sX)
    paulis = [sI, sX, sY, sZ]
    paulis_label = ['I', 'X', 'Y', 'Z']
    obs = []
    coeffs = []
    matrix = []
    
    for term in itertools.product(paulis, repeat=n):
        matrices = [pauli for pauli in term]
        coeff = np.trace(ft.reduce(np.kron, matrices) @ H) / N 
        coeff = np.real_if_close(coeff).item()
        
        # Hilbert-Schmidt-Product
        if not np.allclose(coeff, 0): 
            coeffs.append(coeff)
            obs.append(''.join([paulis_label[[i for i, x in enumerate(paulis) 
            # if np.all(x == t)][0]]+str(idx) for idx, t in enumerate(reversed(term))]))  ## I0X1 format 
            if np.all(x == t)][0]] for idx, t in enumerate(reversed(term))]))  ## IX format
            matrix.append(ft.reduce(np.kron, matrices))

    return obs, coeffs , matrix



if __name__ == "__main__":
    print(get_Fmat(2))
    print(getHamiltonian(2))
    print(decompose_ham_to_pauli(getHamiltonian(2)))

    

    pass

