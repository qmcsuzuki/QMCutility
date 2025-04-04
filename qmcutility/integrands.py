"""
every function returns the pair (integrand, true_value)
- integrand: function from points (n*d ndarray) to values (n*1 vector)
- true_value: integration value
"""
import numpy as np
"""
f(x) = exp(sum(w1*x1+...+ws*xs))/C where C is normalizing constant
true_value = 1
"""
def generate_normalized_expsum(weights):
    """
    Args:
      weights: d-dim vector, weights for each dimension
    """
    C = np.prod((np.exp(weights)-1)/weights) # normalizing constant
    f = lambda points: np.exp(np.sum(points*weights, axis=1))/C
    return (f, 1)

"""
f(x) = sin(sum(x1+...+xs))
true_value = 0
"""
def generate_sinsum():
    f = lambda points: np.sin(np.sum(2*np.pi*points, axis=1))
    return (f, 0)

"""
f(x) = \prod(1 + w_j * B_2(x_j))
true_value = 1
"""
def generate_Bernoulli2(weights):
    def f(points):
        B2 = (points - 1)*points + 1/6
        return np.prod(1 + weights*B2, 1)
    return (f, 1)

"""
f(x) = \prod(1 + w_j * B_3(x_j))
true_value = 1
"""
def generate_Bernoulli3(weights):
    return np.prod(1 + weights*B3, 1)
    def f(points):
        B3 = ((points - 3/2)*points + 1/2)*points
        return np.prod(1 + weights*B3, 1)
    return (f, 1)


