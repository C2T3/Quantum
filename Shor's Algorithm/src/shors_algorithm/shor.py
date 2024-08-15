""" 
    This file implements Shor's algorithm and its subroutines.
    By : Mathis Beaudoin (Summer 2024)
"""


from .order_finding import order_finding, hardcoded_order_finding
from math import gcd, floor, log
from random import randint
from typing import Callable
import time


def pure_power(N : int) -> tuple[int, int]:
    """Checks if N is a pure power (N = a^b for a >= 1 and b >= 2).

    Args:
        N (int): A positive integer.

    Returns:
        tuple[int, int]: a and b if N is a pure power, otherwise N and 1.
    """
    l = floor(log(N, 2)) + 1
    y = log(N,2)

    for b in range(2, l):
        x = y/b
        u1 = floor(pow(2,x))
        u2 = u1 + 1

        if pow(u1, b) == N:
            return u1, b
        
        elif pow(u2, b) == N:
            return u2, b
        
    return N, 1


def shors_algorithm(N : int, t : int, method : Callable, upper_reg = None, lower_reg = None, hardcoded = False, one_qubit = False, noise = False) -> dict[int, int, bool, int, int, int, int, float, float, bool]:
    """Runs Shor's algorithm to find a factor of N.

    Args:
        N (int): The number to factor.
        t (int): If one_qubit is False, corresponds to the number of qubits (accuracy) used to store the period. Otherwise, the number of modular exponentiation gates applied in the one qubit trick.
        method (Callable): If hardcoded is True, a function that builds a specific version of the order finding circuit (beauregard, pavlidis, ...). Otherwise, the modular exponentiation method used.
        upper_reg (list, optional): The qubits where the value of the period will be stored if one_qubit is False. Otherwise, the qubit used in the one qubit trick for the measurements. If hardcoded is True, this field can be left blank. Defaults to None.
        lower_reg (list, optional): The qubits where the modular exponentiaton is done. If hardcoded is True, this field can be left blank. Defaults to None.
        hardcoded (bool, optional): If True, will use a hardcoded version of the order finding circuit (beauregard, pavlidis, ...). Defaults to False.
        one_qubit (bool, optional): If True, will do the one qubit trick. Defaults to False.
        noise (bool, optional): If True, will use Fake127QPulseV1 as a backend to add noise to the simulation. Defaults to False.
        
    Returns:
        dict: - "N" (int): The number that is being factored.
            - "t" (int): If one_qubit is False, corresponds to the number of qubits (accuracy) used to store the period. Otherwise, the number of modular exponentiation gates applied in the one qubit trick.
            - "one_qubit" (bool): If True, the algorithm used the one qubit trick.
            - "factor_found" (int, tuple): Either one integer which is for sure a factor of N or a tuple with a base and exponent (base, exponent) in case N is a pure power.
            - "n_qubits" (int): The number of qubits used when a factor is found.
            - "depth" (int): The circuit's depth when a factor is found.
            - "n_attempts" (int): The number of attempts made to find a factor.
            - "simulation_time (float): The time it took for the simulation to run when a factor is found.
            - "total_time" (float): The time it took for the algorithm to terminate.
            - "noise" (bool): If True, noise was added to the simulation (Fake127QPulseV1).
    """
    results = {"N" : N, "t" : t, "one_qubit" : one_qubit, "factor_found" : 0.0, "n_qubits" : 0.0, "depth" : 0.0, "n_attempts" : 0.0, "simulation_time" : 0.0, "total_time" : 0.0, "noise" : noise}
    n_attempts = 1

    start_total_time = time.time()
    while True:
        #If N is even
        if not (N & 1): 
            results["factor_found"] = 2
            results["n_attempts"] = n_attempts
            results["total_time"] = time.time() - start_total_time
            return results

        #If N is a pure power
        x,y = pure_power(N)
        if x != N: 
            results["factor_found"] = (x,y)
            results["n_attempts"] = n_attempts
            results["total_time"] = time.time() - start_total_time
            return results

        #Choose random basis
        a = randint(1, N-1)
        Gcd = gcd(a, N)
        if Gcd > 1: 
            results["factor_found"] = Gcd
            results["n_attempts"] = n_attempts
            results["total_time"] = time.time() - start_total_time
            return results
        
        #Order finding
        start_simulation_time = time.time()
        order, n_qubits, depth = hardcoded_order_finding(N, a, t, method, one_qubit, noise) if hardcoded else order_finding(N, t, method, a, upper_reg, lower_reg, one_qubit, noise)

        #If order is even
        if not (order & 1) and order != 0:
            z = a ** (order // 2)
            if z % N != N-1:
                #Verify the factors
                factors = []
                factor_1 = gcd(z + 1, N)
                factor_2 = gcd(z - 1, N)

                if factor_1 != 1 and factor_1 != N and (N / factor_1).is_integer(): 
                    factors.append(factor_1)
                if factor_2 != 1 and factor_2 != N and (N / factor_2).is_integer(): 
                    factors.append(factor_2)

                #If factors are found
                if len(factors) != 0:
                    results["factor_found"] = factors[0]    #We only save one factor, even if two can be found in one iteration
                    results["n_qubits"] = n_qubits
                    results["depth"] = depth
                    results["n_attempts"] = n_attempts
                    results["simulation_time"] = time.time() - start_simulation_time
                    results["total_time"] = time.time() - start_total_time
                    return results

        n_attempts += 1


def classical_shor(N : int) -> dict[int, int, int, float, float]:
    """Classical version of Shor's algorithm. The only difference is that the period is found classically.

    Args:
        N (int): The number to factor.

    Returns:
        dict: - "N" (int): The number that is being factored. 
            - "factor_found" (int, tuple): Either one integer which is for sure a factor of N or a tuple with a base and exponent (base, exponent) in case N is a pure power.
            - "n_attempts" (int): The number of attempts made to find a factor.
            - "processing_time (float): The time it took for the order finding portion of the algorithm when a factor is found.
            - "total_time" (float): The time it took for the algorithm to terminate.
    """
    results = {"N" : N, "factor_found" : 0.0, "n_attempts" : 0.0, "processing_time" : 0.0, "total_time" : 0.0}
    n_attempts = 1

    start_total_time = time.time()
    while True:
        #If N is even
        if not (N & 1): 
            results["factor_found"] = 2
            results["n_attempts"] = n_attempts
            results["total_time"] = time.time() - start_total_time
            return results

        #If N is a pure power
        x,y = pure_power(N)
        if x != N: 
            results["factor_found"] = (x,y)
            results["n_attempts"] = n_attempts
            results["total_time"] = time.time() - start_total_time
            return results

        #Choose random basis
        a = randint(1, N-1)
        Gcd = gcd(a, N)
        if Gcd > 1: 
            results["factor_found"] = Gcd
            results["n_attempts"] = n_attempts
            results["total_time"] = time.time() - start_total_time            
            return results
        
        #Order finding
        start_processing_time = time.time()
        order = 0
        for i in range(2, N):
            if a**i % N == 1:
                order = i
                break
        
        #If order is even
        if not (order & 1):
            z = a ** (order // 2)
            if z % N != N-1:
                #Verify the factors
                factors = []
                factor_1 = gcd(z + 1, N)
                factor_2 = gcd(z - 1, N)

                if factor_1 != 1 and factor_1 != N and (N / factor_1).is_integer(): 
                    factors.append(factor_1)
                if factor_2 != 1 and factor_2 != N and (N / factor_2).is_integer(): 
                    factors.append(factor_2)

                #If factors are found
                if len(factors) != 0:
                    results["factor_found"] = factors[0]    #We only save one factor, even if two can be found in one iteration
                    results["n_attempts"] = n_attempts
                    results["processing_time"] = time.time() - start_processing_time
                    results["total_time"] = time.time() - start_total_time
                    return results

        n_attempts += 1