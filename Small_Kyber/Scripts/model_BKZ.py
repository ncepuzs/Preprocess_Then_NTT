from math import *
import functools

log_infinity = 9999

def delta_BKZ(b):
    """ The root hermite factor delta of BKZ-b
    """
    return ((pi*b)**(1./b) * b / (2*pi*exp(1)))**(1./(2.*b-2.))


def svp_plausible(b):
    """ log_2 of best plausible Quantum Cost of SVP in dimension b
    """
    return b *log(sqrt(4./3))/log(2)   # .2075 * b 


def svp_quantum(b):
    """ log_2 of best plausible Quantum Cost of SVP in dimension b
    """
    return b *log(sqrt(13./9))/log(2)   # .265 * b  [Laarhoven Thesis]


def svp_classical(b):
    """ log_2 of best known Quantum Cost of SVP in dimension b
    """
    return b *log(sqrt(3./2))/log(2)    # .292 * b [Becker Ducas Laarhoven Gama]


def nvec_sieve(b):
    """ Number of short vectors outputted by a sieve step of blocksize b
    """
    return b *log(sqrt(4./3))/log(2)    # .2075 * b



def BKZ_first_length(q, nq, n1, b):
    """ Simulate the length of the shortest expected vector in the first b-block
        after a BKZ-b reduction, ignoring the q-vector if some remains.
    """

    d=nq+n1
    delta = delta_BKZ(b)
    l = q**(1.*nq/d)) * delta**d
    return l

def BKZ_last_block_length(q, nq, n1, b):
    """ Simulate the length of the expected Gram-Schmidt vector at position d-b (d = n+m)
        after a BKZ-b reduction.
    """

    (_, _, L) = construct_BKZ_shape(q, nq, n1, b)
    return exp(L[nq + n1 - b])

