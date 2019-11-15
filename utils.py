
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
    using Fermat's little hteorem
    
    Important! p must be a prime number.

    e.g.
    p = 17
    a = 10
    print("Inverse modulo of {} with p={} is {}. a*a^-1={}".format(a, p, InverseFermat(a,p), a*InverseFermat(a,p)%p))
    """
    return fastPowering(a, p-2, p)


