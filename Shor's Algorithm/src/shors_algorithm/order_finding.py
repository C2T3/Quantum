""" 
    This file implements the order finding algorithm and its subroutines.
    By : Mathis Beaudoin (Summer 2024)
"""


from typing import Callable
from .methods.beauregard import mod_exp
from .methods.utils.QFT import build_QFT_dag
from math import floor, log2, pi
from qiskit import QuantumCircuit, transpile, QuantumRegister, ClassicalRegister
from qiskit_aer import AerSimulator
from qiskit.providers.fake_provider import Fake127QPulseV1
from qiskit_aer.noise import NoiseModel
from continuedfractions.continuedfraction import ContinuedFraction


def build_regular_order_finding(N : int, basis : int, t : int, method : Callable, upper_reg : list, lower_reg : list) -> QuantumCircuit:
    """Builds the regular order finding circuit with multiple qubits to store the value of the period. 

    Args:
        N (int): The integer to factor.
        basis (int): The basis used. 
        t (int): The number of qubits (accuracy) used to store the period.
        method (Callable): The modular exponentiation method used.
        upper_reg (list): The qubits where the value of the period will be stored.
        lower_reg (list): The qubits where the modular exponentiaton is done.

    Returns:
        QuantumCircuit: The quantum circuit corresponding to the order finding for the given parameters.
    """
    n = floor(log2(N)) + 1

    qc = QuantumCircuit(len(upper_reg) + len(lower_reg), t)
    qc.h(upper_reg)
    qc.x(lower_reg[0])

    method(qc, basis, N, upper_reg, lower_reg, n)

    qc.append(build_QFT_dag(t), upper_reg)

    qc.measure(upper_reg, upper_reg)

    return qc


def build_one_qubit_order_finding(N : int, basis : int, t : int, method : Callable, lower_reg : list) -> QuantumCircuit:
    n = floor(log2(N)) + 1
    measurments = []

    upper_reg = [0]
    qr = QuantumRegister(1 + len(lower_reg))
    cr = ClassicalRegister(t + 1)
    qc = QuantumCircuit(qr, cr)
    qc.x(lower_reg[0])

    for i in range(t):
        qc.h(upper_reg)
        method(qc, basis**(2**(t-1-i)), N, upper_reg, lower_reg, n)
        for j in range(i):
            if measurments[-1 - j] == 1:
                qc.p(-pi / 2**(j+1), upper_reg)

        qc.h(upper_reg)

        with qc.if_test((cr[i], 1)) as else_:
            qc.x(upper_reg)
            measurments.append(1)
        with else_:
            qc.id(upper_reg)
            measurments.append(0)

    qc.measure(0, t)

    return qc


def build_order_finding(N : int, basis : int, t : int, method : Callable, upper_reg : list, lower_reg : list, one_qubit = False) -> QuantumCircuit:
    """Builds the order finding circuit for the given parameters.

    Args:
        N (int): The number to factor.
        basis (int): The basis used.
        t (int): If one_qubit is False, corresponds to the number of qubits (accuracy) used to store the period. Otherwise, the number of modular exponentiation gates applied in the one qubit trick.
        method (Callable): The modular exponentiation method used.
        upper_reg (list): If one_qubit is True, the qubits where the value of the period will be stored. Otherwise, the qubit used in the one qubit trick for the measurements.
        lower_reg (list): The qubits where the modular exponentiaton is done.
        one_qubit (bool, optional): If True, will do the one qubit trick. Defaults to False.

    Returns:
        QuantumCircuit: The order finding circuit.
    """
    if one_qubit:
        return build_one_qubit_order_finding(N, basis, t, method, lower_reg)
    else:
        return build_regular_order_finding(N, basis, t, method, upper_reg, lower_reg)


def order_finding_circuit_beauregard(N : int, basis : int, t : int, one_qubit = False) -> QuantumCircuit:
    """Builds the circuit for the order finding algorithm using beauregard's modular exponentiation.

    Args:
        N (int): The integer to factor.
        basis (int): The basis used.
        t (int): If one_qubit is False, corresponds to the number of qubits (accuracy) used to store the period. Otherwise, the number modular exponentiation gates applied in the one qubit trick.
        one_qubit (bool, optional): If True, will produce the one qubit trick circuit. Defaults to False.

    Returns:
        QuantumCircuit: The circuit for the order finding algorithm.
    """
    n = floor(log2(N)) + 1
    upper_reg = None if one_qubit else [i for i in range(t)] 
    lower_reg = [i for i in range(1, 2*n+3)] if one_qubit else [i for i in range(t, t + 2*n + 2)]

    return build_order_finding(N, basis, t, mod_exp, upper_reg, lower_reg, one_qubit)


def order_finding(N : int, t : int, method : Callable, basis : int, upper_reg = None, lower_reg = None, one_qubit = False, noise = False) -> tuple[int, int, int]:
    """Runs the order finding circuit for the given parameters.

    Args:
        N (int): The number to factor.
        t (int):  If one_qubit is False, corresponds to the number of qubits (accuracy) used to store the period. Otherwise, the number of modular exponentiation gates applied in the one qubit trick.
        method (Callable): The modular exponentiation method used.
        basis (int): The basis used.
        upper_reg (list, optional):  If one_qubit is True, the qubits where the value of the period will be stored. Otherwise, the qubit used in the one qubit trick for the measurements. Defaults to None.
        lower_reg (list, optional): The qubits where the modular exponentiaton is done. Defaults to None.
        one_qubit (bool, optional): If True, will do the one qubit trick. Defaults to False.
        noise (bool, optional): If True, will use Fake127QPulseV1 as a backend to add noise to the simulation. Defaults to False.

    Returns:
        tuple[int, int, int]: The order, the number of qubits used and the circuit's depth.
    """
    qc = method(N, basis, t, one_qubit) if lower_reg is None else build_order_finding(N, basis, t, method, upper_reg, lower_reg, one_qubit) 

    noise_model = NoiseModel.from_backend(Fake127QPulseV1()) if noise else None

    # simulator = AerSimulator(noise_model = noise_model, device = "GPU")
    simulator = AerSimulator(noise_model = noise_model)
    transpiled_qc = transpile(qc, simulator)
    measurement = simulator.run(transpiled_qc, shots = 1).result().get_counts()
    
    measurement_value = 0 
    for i, bit in enumerate(list(measurement)[0]):
        if bit == "1":
            measurement_value += 1/(2**(i+1))

    for convergent in ContinuedFraction(measurement_value).convergents.values():
        denom = convergent.denominator
        if (basis**denom) % N == 1:
            return denom, qc.num_qubits, qc.depth()

    return 0, qc.num_qubits, qc.depth()


def hardcoded_order_finding(N : int, basis : int, t : int, hardcoded_method : Callable, one_qubit = False, noise = False) -> tuple[int, int, int]:
    """Runs the order finding circuit with a hardcoded circuit (beauregard, pavlidis, ...).

    Args:
        N (int): The number to factor.
        basis (int): The basis used.
        t (int): If one_qubit is False, corresponds to the number of qubits (accuracy) used to store the period. Otherwise, the number of modular exponentiation gates applied in the one qubit trick.
        hardcoded_method (Callable): A function that builds a specific version of the order finding circuit (beauregard, pavlidis, ...)
        one_qubit (bool, optional): If True, will do the one qubit trick. Defaults to False.
        noise (bool, optional): If True, will use Fake127QPulseV1 as a backend to add noise to the simulation. Defaults to False.

    Returns:
        tuple[int, int, int]: The order, the number of qubits used and the circuit's depth.
    """
    qc = hardcoded_method(N, basis, t, one_qubit)

    noise_model = NoiseModel.from_backend(Fake127QPulseV1()) if noise else None

    # simulator = AerSimulator(noise_model = noise_model, device = "GPU")
    simulator = AerSimulator(noise_model = noise_model)
    transpiled_qc = transpile(qc, simulator)
    measurement = simulator.run(transpiled_qc, shots = 1).result().get_counts()

    measurement_value = 0 
    for i, bit in enumerate(list(measurement)[0]):
        if bit == "1":
            measurement_value += 1/(2**(i+1))

    for convergent in ContinuedFraction(measurement_value).convergents.values():
        denom = convergent.denominator
        if (basis**denom) % N == 1:
            return denom, qc.num_qubits, qc.depth()
        
    return 0, qc.num_qubits, qc.depth()