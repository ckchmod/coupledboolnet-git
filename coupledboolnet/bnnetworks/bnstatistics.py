from coupledboolnet.bnnetworks.bn import bitstoints, bitstointsrobust, steadystates
import numpy as np
import math
from itertools import combinations


def steadystatesrobust(states):

    if (len(states.shape) == 2):
        genes = states.shape[1]
        states = bitstoints(states)
        bins = np.array([i for i in range(2 ** genes)])
        ssd = np.array([sum(states == i) for i in range(2 ** genes)]) / states.size

        return (ssd, bins)

    elif (len(states.shape) == 3):
        numcells = states.shape[0]
        genes = states.shape[1]
        size = states.shape[2]
        bins = np.array([i for i in range(2 ** genes)])
        ssd = np.zeros((numcells, len(bins)), dtype=float)

        bitints = 2 ** np.arange(genes)[::-1]
        for n in range (states.shape[0]):
            statestemp = bitints.dot(states[n,:,:])
            ssd[n, :] = np.array([sum(statestemp == i) for i in range(2 ** genes)]) / size

        return (ssd, bins)

def KLDold(p,q):
    return np.sum(np.where(p != 0, p * np.log(p / q), 0))

def KLDcompute(P,Q):
    temp = np.multiply(P, np.log(np.divide(P,Q)))

    temp[np.isnan(temp)] = 0
    temp[np.isinf(temp)] = 0
    #temp = np.where(np.isinf(temp), temp, 0)
    #temp = np.where(np.isnan(temp), temp, 0)
    return sum(temp)

def binom(n,k):
    return math.factorial(n) // math.factorial(k) // math.factorial(n - k)

def kldpairwise(ssD):
    """
    KLD Symmetric
    """
    pairwise = np.array(list(combinations([i for i in range(ssD.shape[0])], 2)))
    KLDMatrix = np.zeros(pairwise.shape[0])

    for i in range(len(KLDMatrix)):
        P = ssD[pairwise[i, 0],:]
        Q = ssD[pairwise[i, 1],:]
        KLDMatrix[i] = .5 * (KLDcompute(P, Q) + KLDcompute(Q, P))

    return(KLDMatrix)

def lyapunovexp(k, p):
    lambdalyap = np.log(2*k*p*(1-p))
    return lambdalyap