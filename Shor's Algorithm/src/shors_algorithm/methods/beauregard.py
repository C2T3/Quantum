""" 
    This file implements the subroutines used in Beauregard's paper for Shor's algorithm.
    By : Mathis Beaudoin (Summer 2024)
"""


from .utils.draper import phi_add, c_phi_add, cc_phi_add
from .utils.QFT import build_QFT, build_QFT_dag
from qiskit import QuantumCircuit


def cc_phi_add_mod(circuit : QuantumCircuit, a : int, N : int, c_1 : int, c_2 : int, ancillary : int,
                target_qubits : list) -> None:
    """Computes |b + c_1*c_2*a mod N> where "a" is a classical integer, |b> a quantum integer, N a classical integer greater than "a" and |b> and c_1/c_2 are control qubits. An ancillary qubit initialized in state |0> will also be needed and will be restored at the end of the computation.

    Args:
        circuit (QuantumCircuit): The circuit containing the quantum integer, the control qubits and the ancillary qubit.
        a (int): The classical integer to add to the quantum integer.
        N (int): A classical integer > a,b.
        c_1/c_2 (int): Control qubits.
        ancillary (int): An ancillary qubit to help with the computation.
        target_qubits (list): Qubits in the circuit that contain the quantum integer.
    """
    qft_size = len(target_qubits)
    a_mod_N = a % N

    cc_phi_add(circuit, a_mod_N, c_1, c_2, target_qubits)
    phi_add(circuit, N, target_qubits, inverse=True)

    circuit.append(build_QFT_dag(qft_size), target_qubits)
    circuit.cx(target_qubits[-1], ancillary)
    circuit.append(build_QFT(qft_size), target_qubits)

    c_phi_add(circuit, N, ancillary, target_qubits)
    cc_phi_add(circuit, a_mod_N, c_1, c_2, target_qubits, inverse=True)

    circuit.append(build_QFT_dag(qft_size), target_qubits)
    circuit.x(target_qubits[-1])
    circuit.cx(target_qubits[-1], ancillary)
    circuit.x(target_qubits[-1])
    circuit.append(build_QFT(qft_size), target_qubits)

    cc_phi_add(circuit, a_mod_N, c_1, c_2, target_qubits)


def cc_phi_add_mod_dag(circuit : QuantumCircuit, a : int, N : int, c_1 : int, c_2 : int, ancillary : int, target_qubits : list) -> None:
    """Inverse of the cc_phi_add_mod function.

    Args:
        circuit (QuantumCircuit): The circuit containing the quantum integer, the control qubits and the ancillary qubit.
        a (int): The classical integer to add to the quantum integer.
        N (int): A classical integer > a,b.
        c_1/c_2 (int): Control qubits.
        ancillary (int): An ancillary qubit to help with the computation.
        target_qubits (list): Qubits in the circuit that contain the quantum integer.
    """
    qft_size = len(target_qubits)
    a_mod_N = a % N

    cc_phi_add(circuit, a_mod_N, c_1, c_2, target_qubits, inverse=True)

    circuit.append(build_QFT_dag(qft_size), target_qubits)
    circuit.x(target_qubits[-1])
    circuit.cx(target_qubits[-1], ancillary)
    circuit.x(target_qubits[-1])
    circuit.append(build_QFT(qft_size), target_qubits)

    cc_phi_add(circuit, a_mod_N, c_1, c_2, target_qubits)

    c_phi_add(circuit, N, ancillary, target_qubits, inverse=True)

    circuit.append(build_QFT_dag(qft_size), target_qubits)
    circuit.cx(target_qubits[-1], ancillary)
    circuit.append(build_QFT(qft_size), target_qubits)

    phi_add(circuit, N, target_qubits)

    cc_phi_add(circuit, a_mod_N, c_1, c_2, target_qubits, inverse=True)


def c_mult_mod(circuit : QuantumCircuit, a : int, N : int, c : int, ctrl_qubits : list, 
             ancillary : int, target_qubits : list, inverse = False) -> None:
    """Computes |b + c*a*x mod N>. Subroutine for rev_mult_mod. 

    Args:
        circuit (QuantumCircuit): The circuit containing the quantum integer, the control qubits and the ancillary qubit.
        a (int): The classical integer to add to the quantum integer.
        N (int): A classical integer > a,b.
        c (int): A control qubit.
        ctrl_qubits (list): The qubits where x is stored.
        ancillary (int): An ancillary qubit to help with calculations.
        target_qubits (list): Qubits in the circuit that contain the quantum integer.
    """
    if inverse:
        for i, ctrl in enumerate(ctrl_qubits[::-1]):
            cc_phi_add_mod_dag(circuit, (2**(len(ctrl_qubits) - 1 - i))*a, N, c, ctrl, ancillary, target_qubits)
    else:
        for i, ctrl in enumerate(ctrl_qubits):
            cc_phi_add_mod(circuit, (2**i)*a, N, c, ctrl, ancillary, target_qubits)


def rev_mult_mod(circuit : QuantumCircuit, a : int, N : int, c : int, ctrl_qubits : list, 
                 ancillary : int, target_qubits : list) -> None:
    """Reversible modular multiplier used in mod_exp.

    Args:
        circuit (QuantumCircuit): The circuit containing the quantum integer, the control qubits and the ancillary qubit.
        a (int): The classical integer to add to the quantum integer.
        N (int): A classical integer > a,b.
        c (int): A control qubit.
        ctrl_qubits (list): The qubits where x is stored.
        ancillary (int): An ancillary qubit to help with calculations.
        target_qubits (list): Qubits in the circuit that contain the quantum integer.
    """
    c_mult_mod(circuit, a, N, c, ctrl_qubits, ancillary, target_qubits)
    circuit.append(build_QFT_dag(len(target_qubits)), target_qubits)
    circuit.swap(ctrl_qubits, target_qubits[:len(ctrl_qubits)])
    circuit.append(build_QFT(len(target_qubits)), target_qubits)
    c_mult_mod(circuit, pow(a, -1, N), N, c, ctrl_qubits, ancillary, target_qubits, inverse=True)


def mod_exp(circuit : QuantumCircuit, a : int, N : int, ctrl_qubits : list, target_qubits : list, 
            n : int) -> None:
    """Adds the modular exponentiaton gate to a quantum circuit. This the black box in the order 
    finding algorithm.

    Args:
        circuit (QuantumCircuit): The circuit to add the modular exponentiation.
        a (int): A value between 2 and N-1
        N (int): The number to be factored.
        ctrl_qubits (list): The register where the estimation takes place.
        target_qubits (list): The register where to modular exponentiation takes place. We expect the first n qubits to be where the result of each step in the modular exponentiation is stored and the next n + 2 qubits where calculations are done.
        n (int): The number of bits/qubits to write N.
    """
    for i, ctrl in enumerate(ctrl_qubits):
        rev_mult_mod(circuit, a**(2**i), N, ctrl, target_qubits[:n], target_qubits[-1], target_qubits[n:-1])