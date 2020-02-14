from crypt import isPrime, InverseFermat, RandomPrime, fastPowering
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

import numpy as np
# TODO: Implement an interface for all LSSS
# TODO: Implement printing the object __str__ method

class AdditiveLSSS():
    """
    Full Threshold linear secret sharing scheme a.k.a additive secret sharing

    Variables:
        self._size: The size in bits of the generated prime
        self._p: randomly generated prime of size size in bits
        self._n: number of parties to split the secert to
        self._m: M matrix of LSSS
        self._v: V vector of LSSS
        self._k: n values randomly generated s.t their sum modulo p is the secret
        self._secret: the secret to split, a number between 0 and p-1

    Methods:
        self.GenerateShares: Generate a random k vector s.t the sum modulo p
                             of its components is equal to the secret

        self.S: Returns the vector of shares, in this scheme is the same as
                returning the vector K directly since S is the unit matrix M
                multiplied by vector K.

        self.RevealSecret: Reveals the secret. Calculated using the scalar
                           multplication of v (which is ones) and k. This
                           modulated by p is the secret.


    WARNING: This class uses numpy to do matrix operations and maximum precision is
             uint64, however it is not recommended to work with prime numbers of 64
             bits as sometimes causes overflow. This is why by default we work on 
             32 bits. In a better implementation we have to improve this.
    """
    def __init__(self, n, size=32):
        self._size = size
        self._p = RandomPrime(size, 100)
        self._n = n
        self._m = self._init_m()
        self._v = self._init_v()

        # to be initialized
        self._k = None
        self._secret = None

    @property
    def p(self):
        return self._p

    @property
    def n(self):
        return self._n

    @property
    def k(self):
        return self._k
    
    @property
    def secret(self):
        return self._secret

    # Setter for secret, atomatically generates 
    # the random K to get the Shares S
    @secret.setter
    def secret(self, secret):
        self._secret = secret
        self.GenerateShares() 

    def _init_m(self):
        return np.identity(self.n, dtype=np.uint64)

    def _init_v(self):
        return np.ones(self.n, dtype=np.uint64)

    def GenerateShares(self):
        if not self._secret:
            raise ValueError("No secret on the scheme, please feed with secret")
        # generate n random additive shares
        shares = [randrange(self._p) for _ in range(self._n-1)]
        shares.append((self._secret-sum(shares))%self._p)
        self._k = np.array(shares, dtype=np.uint64)

    def S(self):
        return self._m.dot(self._k)

    def RevealSecret(self):
        return int(np.dot(self._v, self._k)%self._p)


class ShamirLSSS():
    def __init__(self, n, t, size=32):
        self._size = size
        self._p = RandomPrime(size, 100)
        self._n = n
        self._t = t

        # Initialize in this class
        self._m = self._init_m()
        self._v = self._init_v()

        # to be initialized
        self._k = None
        self._secret = None

    @property
    def p(self):
        return self._p

    @property
    def n(self):
        return self._n

    @property
    def k(self):
        return self._k

    @property
    def t(self):
        return self._t

    @property
    def m(self):
        return self._m

    @property
    def v(self):
        return self._v

    @property
    def secret(self):
        return self._secret

    @secret.setter
    def secret(self, secret):
        self._secret = secret
        self.GenerateShares() 

    def _init_m(self):
        m = np.ones((self._n, self._t+1), dtype=np.uint64)
        for row in range(self._n):
            for col in range(self._t+1):
                m[row][col] = fastPowering(row+1, col, self._p)
        return m

    def _init_v(self):
        v = np.zeros(self._t+1, dtype=np.uint64)
        v[0] = 1
        return v

    def GenerateShares(self):
        if not self._secret:
            raise ValueError("No secret on the scheme, please feed with secret")
        # TODO: should modify function ShamnirRandomPolynomial to generate t+1
        shares = ShamirRandomPolynomial(self._secret, self._t, self._p)
        self._k = np.array(shares, dtype=np.uint64)

    def S(self):
        return self._m.dot(self._k)

    def RevealSecret(self):
        return int(np.dot(self._v, self._k)%self._p)














