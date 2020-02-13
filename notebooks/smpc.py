from crypt import isPrime, InverseFermat
from random import randrange
from typing import List, Tuple
import numpy as np

# A class of points defining a polynomial [(x0, y0), (x1, y1), ... (xn, yn)]
Shares = List[Tuple[int, int]]


def ShamirRandomPolynomial(s: int, n: int, p: int) -> List[int]:
    """
    Calculates the coeficients living on a field Zp for a polynomial of degree n
    Input:
        s: secret, independent term of the polynomial
        n: degree of the polynomial, we need t=n+1 to reconstruct
        p: prime number to generate the coeffients
    Output:
        A list of integers with the coefficients [a0, a1, a2,...,an]
    """
    # s is the secret (independent term), n the degree of polynomial and p the prime number for the field
    if not isPrime(p, 40):
        raise ValueError(f"p={p} is not prime, please feed a prime")

    coef = [s]
    while len(coef)<n+1:
        coef.append(randrange(1, p))
    return coef

def PolyEvaluate(coef: List[int], x: int) -> int:
    """
    Horner's method to evaluate polynomial
    coef is a list of the coefficients [a0, a1,..., an], note this is not modulo

    Input:
        coef: list of coefficients of the polynomial
        x: input where we want to evaluate the polynomial
    Returns:
        f(x), the evaluation of the polynomial at x
    """
    n = len(coef)
    coef = coef[::-1]
    r = coef[0]
    for i in range(1, n): 
        r = (r*x + coef[i])
    return r 

def PolyEvaluateModulo(coef: List[int], x: int, p: int) -> int:
    """
    Horner's method to evaluate polynomial
    coef is a list of the coefficients [a0, a1,..., an], note this is 
    modulo p evaluation. Not real valued polynomials

    Input:
        coef: list of coefficients of the polynomial
        x: input where we want to evaluate the polynomial (x<p)
        p: the prime number over which the ring of polynomials is defined
    Returns:
        f(x), the evaluation of the polynomial at x in modulo 
    """
    if not isPrime(p, 40):
        raise ValueError(f"{p} should be a prime number")
    n = len(coef)
    #apply modulo just in case we lie out of the range of the field
    x = x%p 
    coef = coef[::-1]
    r = coef[0]
    for i in range(1, n): 
        r = (r*x + coef[i])%p
    return r 

def DrawPolynomialPoints(shares:int, coef: List[int], p: int,  rng: Tuple[int, int]=(0, 20))-> Shares:
    """
    Draw distinct points from a coefficient of certain degree and coefficients in a list
    Inputs:
        shares: number of points to be drawn randomly
        coef: the list of coefficients of the polynomial. Degree is len(coef)-1
        p: the prime over which the ring is generated
        rng: range of x to be selected randomly.

    returns:
        A list of tuples [(x0, f(x0)), (x1, f(x1)),...,(x_npoints-1, f(x_npoints-1))]
        those are random points (not repeated) sampled from the polynomial with x in the range
        of the parameter rng

    """
    if not isPrime(p, 40):
        raise ValueError(f"{p} is not a prime, please enter the prime that generates the coefficients")
    if len(rng)!=2:
        raise ValueError(f"rng value is incorrect, must be a tuple (min, max) in the range of x to draw")

    degree = len(coef)-1
    x = []
    while len(x)<shares:
        i = randrange(rng[0], rng[1])
        
        # avoid sharing the secret (i=0) and repeat shares
        if i not in x and i!=0:
            x.append(i)
    
    # apply modulo operation to the evaluation, otherwise we evaluate regular polynomial
    y = [PolyEvaluateModulo(coef, point, p) for point in x]
    
    return list(zip(x, y))

def LagrangeInterpolation(x: int, shares: Shares, p: int) -> int:
    """
    Given t shares calculate the interpolation of the polynomial of degree
    t-1 at point x.
    Input: 
        x: the point where we want to calculate f(x)
        shares: a set of values (xi, yi) randomly sampled from the polynomial
    Returns:
        the value f(x) for the polynomial with degree len(points)-1 (fully determined)
    """
    if not isPrime(p, 40):
        raise ValueError(f"{p} is not a prime, please enter the prime with which you have generated the random points")

    #Check all points are distinct
    distinct_xi = set([x for x,y in shares])
    if len(shares)!=len(distinct_xi):
        raise ValueError("Points at which we evaluate the polynomial must be all different")
    
    summands = []
    for i, (xi, yi) in enumerate(shares):
        val = 1
        for j, (xj, yj) in enumerate(shares):
            if i!=j:
                val*=(x-xj)*InverseFermat(xi-xj,p)%p
                
        val*=yi%p
        summands.append(val)
    
    # cast it to integer (has to be since x and y are all integers and f is polynomial)
    # same time to avoid numerical errors we round the number
    return int(round(sum(summands)))%p

def RevealSecret(shares: Shares, p: int) -> int:
    """
    Returns the coefficient a0 from the polynomial constructed
    using the shares
    Input: 
        shares: a set of points (x, y) drawn from the polynomial
        p: prime number that generates the field
    Output:
        The intercept of the polynomial a0 or equivalently, the secret
    """
    return LagrangeInterpolation(0, shares, p)
