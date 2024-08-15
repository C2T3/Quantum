""" 
    This file contains the code to build the different types of adders based on Draper's work.
    By : Mathis Beaudoin (Summer 2024)
"""


from qiskit import QuantumCircuit
from qiskit.circuit.library.standard_gates import PhaseGate
from math import pi


def phi_add(circuit : QuantumCircuit, a : int, target_qubits : list, inverse = False) -> None:
    """Adds a classical integer "a" to a quantum integer |b>. Little-endian representation is used 
    (MSB is the last qubit in "target_qubits"). |b> must have enough qubits to contain the sum in 
    order to prevent overflow. phi_add(a, |b>) = |a + b>. 

    Args:
        circuit (QuantumCircuit): The circuit containing the quantum integer.
        a (int): The classical integer to add to the quantum integer.
        target_qubits (list): Qubits in the circuit that contain the quantum integer.
        inverse (bool, optional): If True, will produce the inverse circuit. Defaults to False.
    """
    inv = -1 if inverse else 1
    for i, qubit in enumerate(target_qubits):
        circuit.p(2*pi*a*inv / 2**(i+1), qubit)


def c_phi_add(circuit : QuantumCircuit, a : int, c : int, target_qubits : list, 
              inverse = False) -> None:
    """Adds a classical integer "a" to a quantum integer |b> controlled by the value k of qubit c. 
    Little-endian representation is used (MSB is the last qubit in "target_qubits"). |b> must have 
    enough qubits to contain the sum in order to prevent overflow. c_phi_add(a, |c>, |b>) = |c> |b + 
    k*a>. 

    Args:
        circuit (QuantumCircuit): The circuit containing the quantum integer and the control qubit.
        a (int): The classical integer to add to the quantum integer.
        c (int): The control qubit.
        target_qubits (list): Qubits in the circuit that contain the quantum integer.
        inverse (bool, optional): If True, produces the inverse circuit. Defaults to False.
    """
    inv = -1 if inverse else 1
    for i,qubit in enumerate(target_qubits):
        gate = PhaseGate(2*pi*a*inv / 2**(i+1)).control(1)
        circuit.append(gate, [c, qubit])


def cc_phi_add(circuit : QuantumCircuit, a : int, c_1 : int, c_2 : int, target_qubits : list, 
               inverse = False) -> None:
    """Adds a classical integer "a" to a quantum integer |b> controlled by the value k_1 and k_2 of 
    qubit c_1 and c_2. Little-endian representation is needed (MSB is the last qubit in 
    "target_qubits"). |b> must have enough qubits to contain the sum in order to prevent overflow. 
    cc_phi_add(a, |c_1>, |c_2>, |b>) = |c_1> |c_2> |b + k_1*k_2*a>.

    Args:
        circuit (QuantumCircuit): The circuit containing the quantum integer and the control qubits.
        a (int): The classical integer to add to the quantum integer.
        c_1/c_2 (int): The control qubits.
        target_qubits (list): Qubits in the circuit that contain the quantum integer.
        inverse (bool, optional): If True, produces the inverse circuit. Defaults to False.
    """
    inv = -1 if inverse else 1
    for i,qubit in enumerate(target_qubits):
        gate = PhaseGate(2*pi*a*inv / 2**(i+1)).control(2)
        circuit.append(gate, [c_1, c_2, qubit])


def q_adder(circuit : QuantumCircuit, ctrl_qubits : list, target_qubits : list, 
            inverse = False) -> None:
    """Adds two quantum integers together. Little-endian representation is used (MSB is the last qubit 
    in "ctrl_qubits" and "target_qubits"). |b> must have enough qubits to contain the sum in order to 
    prevent overflow. q_adder(|a>, |b>) = |a> |a + b>.

    Args:
        circuit (QuantumCircuit): The circuit containing the quantum integers.
        ctrl_qubits (list): The qubits in the circuit containing the quantum integer to add.
        target_qubits (list): The qubits in the circuit containing the quantum integer where the 
        addition will occur.
        inverse (bool): If True, produces the inverse circuit. Defaults to False.
    """
    inv = -1 if inverse else 1
    size = len(ctrl_qubits)
    for i in range(size):
        for j in range(size - i):
            gate = PhaseGate(2*pi*inv / 2**(j+1)).control(1)
            circuit.append(gate, [ctrl_qubits[size - 1 - (j+i)], target_qubits[size - 1 - i]])