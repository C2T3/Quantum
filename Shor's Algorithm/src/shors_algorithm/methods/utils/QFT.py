""" 
    This file contains the code to build the QFT and inverse QFT.
    By : Mathis Beaudoin (Summer 2024)
"""


from qiskit import QuantumCircuit
from qiskit.circuit.library.standard_gates import PhaseGate
from math import pi


def build_QFT(n_qubits : int, swap_beginning=False, swap_end=False) -> QuantumCircuit:
    """Builds the n qubit QFT.

    Args:
        n_qubits (int): The number of qubits.
        insert_barriers (bool, optional): To separate each part of the QFT. Defaults to False.
        do_swaps (bool, optional): Add swaps at the end of the circuit. Defaults to True.

    Returns:
        QuantumCircuit: The quantum circuit for the QFT.
    """
    qft_circuit = QuantumCircuit(n_qubits)

    if swap_beginning:
        iterations = int(n_qubits/2)
        for k in range(iterations):
            qft_circuit.swap(k,n_qubits-k-1)

    for i in range(n_qubits, 0, -1):
        qft_circuit.h(i-1)

        for j in range(1,i): 
            phase = pi / 2**(j)
            ctrl_u = PhaseGate(phase).control(1)
            qft_circuit.append(ctrl_u, [i-j-1, i-1])
            
    if swap_end:
        iterations = int(n_qubits/2)
        for k in range(iterations):
            qft_circuit.swap(k,n_qubits-k-1)

    return qft_circuit


def build_QFT_dag(n_qubits : int, swap_beginning=False, swap_end=False) -> QuantumCircuit:
    """Builds the inverse QFT.

    Args:
        n_qubits (int): The number of qubits.
        insert_barriers (bool, optional): To separate each part of the inverse QFT. Defaults to False.
        do_swaps (bool, optional): Add swaps at the end of the circuit. Defaults to True.

    Returns:
        QuantumCircuit: The inverse QFT.
    """
    return build_QFT(n_qubits, swap_beginning, swap_end).inverse()