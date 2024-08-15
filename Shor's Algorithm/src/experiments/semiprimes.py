""" 
    This file implements functions to test Shor's algorithm on semiprimes.
    By : Mathis Beaudoin (Summer 2024)
"""

#For import purposes
import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.abspath(os.path.join(current_dir, '..'))
sys.path.append(src_dir)

#Path for Windows
# experiments_path = src_dir + "\\experiments\\"

#Path for Linux
experiments_path = src_dir + "/experiments/"

from shors_algorithm.shor import classical_shor, shors_algorithm
from shors_algorithm.order_finding import order_finding_circuit_beauregard
import numpy as np
import csv
from sympy import primefactors
import time
from datetime import datetime


def factorize_semiprime_sympy(N : int, path : str, add = False) -> None:
    """Factorizes a semiprime number using sympy.primefactors and adds a new line with the data to an existing csv file.

    Args:
        N (int): The semiprime number to factor.
        path (str): The path to the csv file.
        add (bool, optional): If True, will add a new line with the data to the csv file. Defaults to False.
    """
    start = time.time()
    factors = primefactors(N)
    total_time = time.time() - start

    if add:
        with open(path, mode='a', newline='') as file:
            writer = csv.writer(file)
            # The header in this case is : ["N", "factor_found", "total_time"]
            writer.writerow([N, factors[0], total_time])


def factorize_semiprime_classical_shor(N : int, path : str, add = False) -> None:
    """Factorizes a semiprime number using the classical version of Shor's algorithm and adds a new line with the data to an existing csv file.

    Args:
        N (int): The semiprime number to factor.
        path (str): The path to the csv file.
        add (bool, optional): If True, will add a new line with the data to the csv file. Defaults to False.
    """
    results = classical_shor(N)

    if add:
        with open(path, mode='a', newline='') as file:
            writer = csv.writer(file)
            # The header in this case is : ["N", "factor_found", "n_attempts", "processing_time", "total_time"]
            writer.writerow([N, results["factor_found"], results["n_attempts"], results["processing_time"], results["total_time"]])


def factorize_semiprime_quantum_shor(N : int, t : int, one_qubit : bool, noise : bool, path : str, add = False) -> None:
    """Factorizes a semiprime number using the quantum version of Shor's algorithm (beauregard) and adds a new line with the data to an existing csv file.

    Args:
        N (int): The semiprime number to factor.
        t (int): If one_qubit is False, corresponds to the number of qubits (accuracy) used to store the period. Otherwise, the number of modular exponentiation gates applied in the one qubit trick.
        one_qubit (bool): If True, will do the one qubit trick.
        noise (bool): If True, adds noise to the simulation (Fake127QPulseV1).
        path (str): The path to the csv file.
        add (bool, optional): If True, will add a new line with the data to the csv file. Defaults to False.
    """
    results = shors_algorithm(N, t, order_finding_circuit_beauregard, hardcoded=True, one_qubit=one_qubit, noise=noise)

    if add:
        with open(path, mode='a', newline='') as file:
            writer = csv.writer(file)
            # The header in this case is : ["N", "t", "one_qubit", "factor_found", "n_qubits", "depth", "n_attempts", "simulation_time", "total_time"]
            writer.writerow([N, t, one_qubit, results["factor_found"], results["n_qubits"], results["depth"], results["n_attempts"], results["simulation_time"], results["total_time"], noise])


def test_semiprimes(n_reps : int, path : str, method : str, limit = None, t = None, one_qubit = None, noise = False) -> None:
    """From a pre-computed list of the first 86 157 semiprimes (all semiprimes up to 500 000), factorizes them up to a certain limit and uploads related data to a new csv file. For the classical methods, t, one_qubit and noise don't need to be provided.

    Args:
        n_reps (int): The number of times a semiprime will be factored (in order to do an average later in the analysis).
        path (str): The path where the new csv file will be created (from experiments folder).
        method (str): Either "sympy", "classical_shor" or "quantum_shor".
        limit (int, optional): The number of semiprimes tested. For example, limit = 10 -> the first 10 semiprimes are tested. If limit = None, will factorize all 86 157 semiprimes from the list. Defaults to None.
        t (int, optional): If one_qubit is False, corresponds to the number of qubits (accuracy) used to store the period. Otherwise, the number of modular exponentiation gates applied in the one qubit trick. Defaults to None.
        one_qubit (bool, optional): If True, will do the one qubit trick. Defaults to None.
        noise (bool, optional): If True, adds noise to the simulation (Fake127QPulseV1). Defaults to False.
    """
    #Importing semiprimes
    load = np.load(experiments_path + 'list_semiprimes.npz')
    semiprimes = load['odd_semiprimes'][:limit] if limit is not None else load['odd_semiprimes']

    #Choosing the appropriate header
    if method == "sympy": header = ["N", "factor_found", "total_time"]
    elif method == "classical_shor": header = ["N", "factor_found", "n_attempts", "processing_time", "total_time"]
    elif method == "quantum_shor": header = ["N", "t", "one_qubit", "factor_found", "n_qubits", "depth", "n_attempts", "simulation_time", "total_time", "noise"]

    #Creating a new csv file
    title = method + "_(" + time.strftime("%Y_%m_%d") + "_" + datetime.now().strftime("%H_%M") + ")" + ".csv"
    with open(experiments_path + path + title, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(header)

    #Factorizing semiprimes
    for semiprime in semiprimes:
        for _ in range(n_reps):
            if method == "sympy": factorize_semiprime_sympy(semiprime, experiments_path + path + title, True)
            elif method == "classical_shor": factorize_semiprime_classical_shor(semiprime, experiments_path + path + title, True)
            elif method == "quantum_shor": factorize_semiprime_quantum_shor(int(semiprime), t, one_qubit, noise, experiments_path + path + title, True)



### SCRIPT TO RUN ON SIMULATOR ###
# test_semiprimes(10, "results/", "quantum_shor", t = 10, one_qubit = False, noise = False)