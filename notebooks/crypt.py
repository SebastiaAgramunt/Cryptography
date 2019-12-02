from random import randrange


def xgcd(a, b):
    """
    Extended greatest common divisor.
    ttps://en.wikipedia.org/wiki/Extended_Euclidean_algorithm

    Solving equation au+bv=gcd(a,b)
    result is: (g,u,v) where g is greatest common divisor and u,v the solution of the eq above

    e.g.
    a = 16261
    b = 85652
    g, u, v = xgcd(a,b)
    print("the solution to the equation au+bv=gcd(a,b) where (a,b)=({},{})".format(a,b))
    print("is g, u, v = ({}, {}, {}) where g is the gcd\n".format(g,u,v))
    print("therefore the equation holds:")
    print("{}u+{}v={} for the above mentioned values of u and v\n".format(a,b,g))
    """
    u0, u1, v0, v1 = 0, 1, 1, 0
    while a != 0:
        q, b, a = b // a, a, b % a
        v0, v1 = v1, v0 - q * v1
        u0, u1 = u1, u0 - q * u1
    return b, u0, v0

def primes_sieve_eratosthenes(limit):
    """
    Calculate list of primes until an integer "limit"
    https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes

    e.g. Calculate primes up to 100
    n = 100
    primes = list(primes_sieve_eratosthenes(n))
    print("All prime numbers bellow {} are:\n{}".format(n,primes))

    """
    a = [True] * limit                          
    a[0] = a[1] = False
    for i, isprime in enumerate(a):
        if isprime:
            yield i
            for n in range(i*i, limit, i):
                a[n] = False

def InvertibleNumbers(m):
    """
    Given a natural number m it returns the invertible numbers modulo(m)

    e.g.
    m = 24
    inv = InvertibleNumbers(m)
    print("The invertible numbers (product wise) of Z/mZ m={} are\n{}".format(m,inv))

    m = 11
    inv = InvertibleNumbers(m)
    print("The invertible numbers (product wise) of Z/mZ m={} are\n{}".format(m,inv))
    """
    inv = []
    for i in range(m):
        g, _, _ = xgcd(i,m)
        if g==1:
            inv.append(i)
    return inv

def InverseMod(a, m):
    """
    Given natural number a and m returns a^{-1}(mod m). The inverse modulo m of a

    e.g.
    a = 10
    m = 7
    inv_a = InverseMod(a,m)

    print("inverse of {} in modulo {} is {}\na*inv_a = {}".format(a,m,inv_a,a*inv_a%m))
    """
    g, u, _ = xgcd(a,m)
    if g!=1:
        print("{} has no inverse on modulo {}".format(a,m))
        return None
    return u

def fastPowering(a, k, m):
    """
    Fast powering algorithm
    a in $Z_m$ and integers 0 <= k <= m
    
    returns:
    a^k (mod m)
    
    Note: exactly the same as pow(a,k,p) implemented in python
    e.g.
    p = 7
    for k in range(3*p):
        print("3^{} (mod {}) = {}".format(k, p, fastPowering(3,k,p)))
    """
    b = 1
    if k == 0:
        return b
    A = a
    # If the least significant bit is 1, $a^1 = a$
    if 1 & k:
        b = a
    k = k >> 1
    while k:
        A = (A**2) % m
        if 1 & k:
            b = (b * A) % m
        k = k >> 1
    return b

def InverseFermat(a, p):
    """
    Given natural number a prime number p returns a^{-1}(mod p). The inverse
    modulo p of a.
    This way of getting the inverse of a number modulo a prime number is found
    using Fermat's little theorem
    
    Important! p must be a prime number.

    e.g.
    p = 17
    a = 10
    print("Inverse modulo of {} with p={} is {}. a*a^-1={}".format(a, p, InverseFermat(a,p), a*InverseFermat(a,p)%p))
    """
    # p has to be a prime number
    if not isPrime(p, 40):
        return None
    
    return fastPowering(a, p-2, p)

def EulerTotient(m):
    """
    Given a number m, gives the number of generators of the group, i.e. 
    the nubmer of elements that can generate all the group by powering
    """
    c = 0
    for i in range(m):
        g, _, _ = xgcd(i, m)
        if g == 1:
            c+=1
    return c 

def isPrime(n, m=5):
    """
    Check primality of a number n using Miller-Rabin algorithm:
    Book An introduction to mathemahical cryptography (Hoffstein, Pipher and Silverman)
    page 127

    n: The number we want to check primality
    m: Number of checks to perform

    Returns: False if the number is composite
             True if the number is probable to be prime
    """
    
    # Check if the number is in between the first 4 primes
    if n in [2, 3, 5, 7]:
        return True
    
    # Take out half of the options, divisible by two is composite
    if n%2 == 0:
        return False
    
    # Write "n-1 = 2^k * q"
    q = n-1
    k = 0
    while q%2 == 0:
        q = q//2
        k += 1
    
    # Loop over number of trials (the more, the better for primality test)
    for _ in range(m):
        
        # a is the potential whithess
        a = randrange(2, n-1)
        gcd, _, _ = xgcd(a, n)
        
        # If gcd of a random number<n and n is not 1 then n is composite.
        if gcd != 1:
            return False
        
        a = fastPowering(a, q, n)
        if a!=1:
            i = 0
            while a != (n-1):
                if i == (k-1):
                    return False
                else:
                    i += 1
                    a = fastPowering(a, 2, n) 
    return True

def generateLargePrime(size=256, m=40):
    """
    Return a random prime number of size bits
    Use Miller-Rabin algorithm to check primality

    m: Number of checks to perform in Miller-Rabin
    Returns: A very probable prime number
    """
    while True:
        n = randrange(2**(size-1), 2**(size))
        if isPrime(n, m):
            return n

# TODO: this function works only when p is prime, gerenalise for n
def findGeneratorPrime(p):
    """
    brute force to find a generator of a group integer modulo p
    with multiplicative operation being p a prime number.
    If p is not prime we would have to take out non invertible elements
    """
   
    while True:
        i = randrange(2, p-1)
        l = len(set([pow(i, j, p) for j in range(p)]))
        if l == p-1:
            return i

def _g(x, c, n):
    # Function x^2+c mod n
    return (fastPowering(x, 2, n) + c)%n

def PollardsRhoFactorisation(n):
    """Finds a non trivial factor for number n
    """

    # If n is prime return 1, the only factor
    if isPrime(n, m=40):
        return 1

    x, c = randrange(2, n), randrange(1, n)
    y, d = x, 1

    while d==1:
        x = _g(x, c, n)
        y = _g(_g(y, c, n), c, n)
        d, _, _ = xgcd(abs(x-y), n)


        if d==n:
            return PollardsRhoFactorisation(n)

    return d












