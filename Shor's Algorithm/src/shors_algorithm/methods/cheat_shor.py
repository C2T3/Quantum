""" 
    This file implements a version of Shor's algorithm where we can cheat to factor very large numbers when
    some factors are known in advance. This is done by choosing a special basis where the period is small. 
    By : Mathis Beaudoin (Summer 2024)
"""


from math import gcd
import time
import random


def pretend_to_factor(p : int, q : int) -> dict[int, list, float]:
    """Manipulated version of Shor's algorithm to factor very large numbers quickly. Will factor N = p*q.

    Args:
        p (int): A first factor.
        q (int): A second factor.

    Returns:
        dict : - "N" (int) : The number factored (p*q).
               - "factors" (list) : The factors found (should technically be p and q).
               - "total_time" (float) : The time it took to factor N.
    """
    start = time.time_ns()

    N = p*q
    p_inv = pow(p, -1, q)
    q_inv = pow(q, -1, p)

    signs = [(1,-1), (-1,1)]
    possible_a = []
    for sign in signs:
        possible_a.append(sign[0]*p*p_inv + sign[1]*q*q_inv)
    possible_a[0] += N

    a = random.choice(possible_a)

    factor_1 = gcd(a + 1, N)
    factor_2 = gcd(a - 1, N)

    factors = []
    if (N / factor_1).is_integer(): factors.append(factor_1)
    if (N / factor_2).is_integer(): factors.append(factor_2)

    end = time.time_ns()

    return {"N" : N, "factors" : factors, "total_time" : float(end - start)}