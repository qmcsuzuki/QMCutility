import numpy as np

"""
make Kronecker points as n*d ndarray
n: number of points
alpha: d-dim generating vector alpha=(a1,a2,...,ad)
"""
def Kronecker_points(n,alpha,randomize=False):
    indices = np.arange(n)[:, np.newaxis]
    if randomize:
        shift = np.random.rand(alpha.shape[0])
        return (indices * alpha + shift) % 1.0
    else:
        return (indices * alpha) % 1.0

def Kronecker_rule(n,alpha,integrand,randomize=False):
    points = Kronecker_points(n,alpha,randomize)
    return np.sum(integrand(points))/n

def Kronecker_type_rule(n,alpha,weights,integrand,randomize):
    points = Kronecker_points(n,alpha,randomize)
    return np.sum(weights*integrand(points))/n

from math import factorial
def SugiharaMorota_rule(n,alpha,rate,integrand,randomize=False):
    C = factorial(2*rate+1)/factorial(rate)**2
    weights = np.array([((i/n)*(1-i/n))**convergence_order for i in range(n)])*C
    return Kronecker_type_rule(n,alpha,weights,integrand,randomize)
