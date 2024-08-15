""" 
    This file implements different subroutines found in Pavlidis' paper for Shor's algorithm. This method uses a lot of qubits and cannot realistically used. Therefore, we leave its subroutines here without using them elsewhere.
    By : Mathis Beaudoin (Summer 2024)
"""


from .utils.QFT import build_QFT, build_QFT_dag
from .utils.draper import phi_add, c_phi_add, q_adder
from qiskit import QuantumCircuit
from math import log, floor


def q_accumulator(circuit : QuantumCircuit, a : int, ctrl_qubits : list, target_qubits : list) -> None:
    """Adds the multiplication of the quantum integer |x> and the classical integer a to the quantum 
    integer |b>. Little-endian representation is used. q_accumulator(a, |x>, |b>) = |x>|b + a*x>. 

    Args:
        circuit (QuantumCircuit): The quantum circuit containing the quantum integers.
        a (int): A classical integer.
        ctrl_qubits (list): The qubits in the circuit containing the quantum integer to multiply.
        target_qubits (list): The qubits in the circuit containing the quantum integer where the 
        addition will occur.
    """
    for i, qubit in enumerate(ctrl_qubits):
        c_phi_add(circuit, a*(2**i), qubit, target_qubits)


def divider(z : int, d : int, n : int, initialize=True, inverse=False) -> QuantumCircuit:
    """Produces a quantum circuit to compute the quotient (q) and remainder (r) of z/d. z and d have 
    to be able to be represented on 2*n bits for the divider to work. #r is on the first sequence of n 
    qubits, q is on the third sequence of n qubits.

    Args:
        z (int): The numerator to divide.
        d (int): The denominator.
        n (int): The number of bits, chosen so that z, d, q and r can be written using 2*n bits.
        inverse (bool, optional): If True, will produce the inverse circuit. Defaults to False.
        initialize (bool, optional): If True, initializes a section of qubits with the value z. 
        Defaults to True.

    Returns:
        QuantumCircuit: A quantum circuit to calculate q and r.
    """
    # Initialization
    l = 1 + floor(log(d, 2))
    m = floor(((2**n)*(2**l-d)-1)/d)
    d_norm = d << (n-l)

    reg_0, reg_1, reg_2, reg_3, reg_4, reg_5, reg_6 = [[k*n + i for i in range(n)] for k in range(0, 7)]
    aqbit = 7*n

    n_qft = build_QFT(n)
    two_n_qft = build_QFT(2*n)
    n_qft_dag = build_QFT_dag(n)
    two_n_qft_dag = build_QFT_dag(2*n)

    # Initialize circuit with z written over reg_0 and reg_1
    qc = QuantumCircuit(7*n + 1)
    if initialize:
        qc.initialize(f'{z:b}'.zfill(2*n), reg_0 + reg_1)

    # Add n_2 to reg_4
    qc.append(n_qft, reg_4)
    q_adder(qc, reg_6[:n-l] + reg_1[:l], reg_4)
    q_adder(qc, reg_0[l:] + reg_6[n-l:], reg_4)

    # Add n_10 to reg_5
    qc.append(n_qft, reg_5)
    q_adder(qc, reg_6[:n-l] + reg_0[:l], reg_5)
    qc.append(n_qft_dag, reg_5)

    # Add n_1 to aqbit
    qc.cx(reg_5[-1], aqbit)

    # Add AND(-n_1 , d_norm - 2^n) to reg_5 to form n_adj
    qc.append(n_qft, reg_5)
    c_phi_add(qc, d_norm - 2**n, aqbit, reg_5)
    qc.append(n_qft_dag, reg_5)

    # Add n_1 to n_2 in reg_4
    c_phi_add(qc, 1, aqbit, reg_4)
    qc.append(n_qft_dag, reg_4)

    # Add m*(n_2 + n_1) + n_adj to reg_5/6
    qc.append(two_n_qft, reg_5 + reg_6)
    q_accumulator(qc, m, reg_4, reg_5 + reg_6)
    qc.append(two_n_qft_dag, reg_5 + reg_6)

    # Copy HIGH(m*(n_2 + n_1) + n_adj) in reg_2
    qc.cx(reg_6, reg_2)

    # Add n_2 + n_1 in reg_2, than subtract n_1 to have q_1 in reg_2
    qc.append(n_qft, reg_2)
    q_adder(qc, reg_4, reg_2)
    c_phi_add(qc, -1, aqbit, reg_2)
    qc.append(n_qft_dag, reg_2)

    # # Copy q_1 in reg_3
    qc.cx(reg_2, reg_3)
    
    # Add dr to reg_0/1
    qc.append(two_n_qft, reg_0 + reg_1)
    q_accumulator(qc, -d, reg_2, reg_0 + reg_1)
    phi_add(qc, -d, reg_0 + reg_1)
    qc.append(two_n_qft_dag, reg_0 + reg_1)

    # # Add q to reg_2
    qc.append(n_qft, reg_2)
    qc.x(reg_0[-1])
    c_phi_add(qc, 1, reg_0[-1], reg_2)
    qc.x(reg_0[-1])
    qc.append(n_qft_dag, reg_2)

    # Restore z into reg0/1
    qc.append(two_n_qft, reg_0 + reg_1)
    phi_add(qc, d, reg_0 + reg_1)
    q_accumulator(qc, d, reg_3, reg_0+reg_1)
    qc.append(two_n_qft_dag, reg_0 + reg_1)
    
    # Restore reg_3
    qc.append(n_qft, reg_3)
    q_adder(qc, reg_4, reg_3, inverse=True)
    c_phi_add(qc, 1, aqbit, reg_3)
    qc.append(n_qft_dag, reg_3)
    qc.cx(reg_6, reg_3)
    
    # Restore reg_6
    qc.append(two_n_qft, reg_5 + reg_6)
    q_accumulator(qc, -m, reg_4, reg_5 + reg_6)
    qc.append(two_n_qft_dag, reg_5 + reg_6)

    # Restore reg_5 and add n_10 to reg_4
    qc.append(n_qft, reg_4)
    c_phi_add(qc, -1, aqbit, reg_4)
    q_adder(qc, reg_6[:n-l] + reg_1[:l], reg_4, inverse=True)
    q_adder(qc, reg_0[l:] + reg_6[n-l:], reg_4, inverse=True)
    q_adder(qc, reg_6[:n-l] + reg_0[:l], reg_4)
    qc.append(n_qft_dag, reg_4)
    qc.append(n_qft, reg_5)
    c_phi_add(qc, -(d_norm - 2**n), aqbit, reg_5)
    q_adder(qc, reg_4, reg_5, inverse=True)
    qc.append(n_qft_dag, reg_5)

    # Restore aqbit
    qc.cx(reg_4[-1], aqbit)

    # Restore reg_4
    qc.append(n_qft, reg_4)
    q_adder(qc, reg_6[:n-l] + reg_0[:l], reg_4, inverse=True)
    qc.append(n_qft_dag, reg_4)

    # Calculate r
    qc.append(two_n_qft, reg_0 + reg_1)
    q_accumulator(qc, -d, reg_2, reg_0 + reg_1)
    qc.append(two_n_qft_dag, reg_0 + reg_1)

    return qc.inverse() if inverse else qc


def modular_acc(y : int, a : int, N : int, n : int, initialize=True) -> QuantumCircuit:
    """Produces a modular accumulator to use in the modular multiplier. Calculates a*y mod N, leaves y on the first n qubits and a*y mod N on the next segment of n qubits. Little-endian representation is used.

    Args:
        y (int): A quantum integer.
        a (int): A classical integer.
        N (int): A classicacl integer.
        n (int): The number of bits/qubits to represent the integers.
        initialize (bool, optional): If True, will initialize the circuit with the value y. Defaults to True.

    Returns:
        QuantumCircuit: An quantum accumulator to use in the modular multiplier.
    """
    reg_a = [i for i in range(n)]
    reg_b = [i for i in range(n, 3*n)]
    reg_c = [i for i in range(n, 8*n + 1)]
    
    reg_d = [i for i in range(n, 2*n)]
    reg_e = [i for i in range(8*n + 1, 9*n + 1)]

    reg_f = [i for i in range(8*n + 1, 10*n + 1)]
    reg_g = [i for i in range(3*n, 8*n + 1)]
    
    qc = QuantumCircuit(10*n + 1)
    if initialize:
        qc.initialize(f'{y:b}'.zfill(n), reg_a)

    qc.append(build_QFT(2*n), reg_b)
    q_accumulator(qc, a, reg_a, reg_b)
    qc.append(build_QFT_dag(2*n), reg_b)

    qc.append(divider(0, N, n, initialize=False), reg_c)

    qc.cx(reg_d, reg_e)
    
    qc.append(divider(0, N, n, inverse=True, initialize=False), reg_f + reg_g)

    qc.append(build_QFT(2*n), reg_f)
    q_accumulator(qc, -a, reg_a, reg_f)
    qc.append(build_QFT_dag(2*n), reg_f)

    return qc


def modular_mult(y : int, a : int, N : int, n : int, initialize=True) -> QuantumCircuit:
    """Produces a circuit to compute the modular multiplication a*y mod N. The result is on the first n qubits (little-endian representation).

    Args:
        y (int): A quantum integer.
        a (int): A classical integer.
        N (int): A classical integer.
        n (int): The number of bits/qubits to represent y, a, N and a*y mod N.

    Returns:
        QuantumCircuit: A quantum circuit to compute a*y mod N.
    """
    qc = QuantumCircuit(10*n + 1)

    if initialize:
    	qc.initialize(f'{y:b}'.zfill(n), [i for i in range(n)])

    qc.append(modular_acc(y, a, N, n, initialize=False), [i for i in range(10*n + 1)])

    qc.swap([i for i in range(n)], [i for i in range(n, 2*n)])

    a_inverse = pow(a, -1, N)

    qc.append(modular_acc(y, a_inverse, N, n, initialize=False).inverse(), [i for i in range(10*n + 1)])

    return qc
