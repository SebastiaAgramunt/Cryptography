from random import randrange


SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,\
 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139,\
 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,\
 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311,\
 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401,\
 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491,\
 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599,\
 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683,\
 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797,\
 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887,\
 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997]

def xgcd(a, b):
    """
    Extended greatest common divisor.
    ttps://en.wikipedia.org/wiki/Extended_Euclidean_algorithm

    Solving equation au+bv=gcd(a,b)
    result is: (g,u,v) where g is greatest common divisor and u,v the solution 
    of the eq above

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

def PrimesSieveEratosthenes(limit):
    """
    Calculate list of primes until an integer "limit"
    https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes

    e.g. Calculate primes up to 100
    n = 100
    primes = list(primes_sieve_eratosthenes(n))
    print("All prime numbers bellow {} are:\n{}".format(n,primes))

    e.g.2 Calculate all prime numbers smaller than 2^16
    n = 1 << 16
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
    return u%m

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
    Book An introduction to mathemahical cryptography (Hoffstein, 
    Pipher and Silverman)
    page 127

    n: The number we want to check primality
    m: Number of checks to perform

    Returns: False if the number is composite
             True if the number is probable to be prime (the larger the m, 
             the more sure we are)
    """
    # Take out half of the options, divisible by two is composite
    if n%2 == 0:
        return False
    
    # Check if the number is in among the list of small primes
    if n in SMALL_PRIMES:
        return True
    
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

def RandomPrime(size=256, m=40):
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

def _g(x, c, n):
    # Function x^2+c mod n
    return (fastPowering(x, 2, n) + c)%n

def PollardsRhoFactorisation(n):
    """Finds a non trivial factor for number n, 
    found factor is not necessarily a prime number(check!)
    """

    # If n is prime return 1, the only factor
    if isPrime(n, m=40):
        return 1

    if n%2==0:
        return 2

    x, c = randrange(2, n), randrange(1, n)
    y, d = x, 1

    while d==1:
        x = _g(x, c, n)
        y = _g(_g(y, c, n), c, n)
        d, _, _ = xgcd(abs(x-y), n)

        if d==n:
            return PollardsRhoFactorisation(n)

    return d

def BruteForceFactorisation(n, primes=SMALL_PRIMES):
    """
    Finds the list of factors of n using a set of prime numbers. And also
    the number of times each number has to be multiplied.
    If the number is too large we may need a larger list of primes.
    you can create large primes (up to for instance 2^16) by imputing 
    primes = list(primes_sieve_eratosthenes(1<<16))
    """
    if isPrime(n,40):
        return [1], [1]
    factors = []
    reps = []
    for prime in primes:
        #cannot factor further
        if prime>n: 
            break
        if n%prime == 0:
            factors.append(prime)
            reps.append(0)
            while n%prime == 0:
                n //= prime
                reps[-1]+=1
    assert n==1, "Cannot factor, primes list is too short"
    return factors, reps

def PrimeFieldGenerator(p, primes=SMALL_PRIMES):
    '''
     input p is a prime number
     a list of primes to do the factorisation of the order
     Finds a generator of the field Fp*, i.e. an element such that
     multiplied by itself (1..,k) times can generate all elements of the group
     ALGORITHM IN BOOK: A Computational Introduction to Number Theory
     and Algebra by Shoup. Page 269.
    if not isPrime(p, 40):
        return None
    '''
    order = p-1
    factors, reps = BruteForceFactorisation(order, primes)

    g=1
    for factor, rep in zip(factors, reps):
        while True:
            a = randrange(1, p) #generator candidate
            b = fastPowering(a, order*InverseFermat(factor, p)%p, p)
            if b!=1:
                exp = InverseFermat(fastPowering(factor, rep, p), p)*order%p
                gamma = fastPowering(a, exp, p)
                g=g*gamma%p
                break
    return g


def isGeneratorBruteForce(g, p):
    l = len(set([pow(g, j, p) for j in range(p)]))
    return True if l == p-1 else False

def LCM(a, b):
    '''
    Computes the least common multiple of two numbers
    '''
    g, _, _ = xgcd(a, b)

    return a*b//g


def RSAKeyGenerator(n=16):
    '''
    RSA key generation
    input bitsize of the key
    output: n, pubK, privK
    '''

    # Generate two random primes of n bits
    p = RandomPrime(n, m=40)
    q = RandomPrime(n, m=40)

    # p and q must be different primes
    while p==q:
      q = RandomPrime(n, m=40)  

    N = p*q
    lm = LCM(p-1, q-1)
    
    while True:
        e = randrange(2, lm)
        g, _, _ = xgcd(e, lm)
        if g==1:
            d = InverseMod(e, lm)
            # return public and private keys
            return (N, e), (N, d)

def RSAEncrypt(m, PublicKey):
    '''
    Input:
        m: message (An integer message)
        PublicKey: A tuple (N, e)
    Returns:
        c: Encrypted message
    '''
    N = PublicKey[0]
    e = PublicKey[1]

    return fastPowering(m, e, N)

def RSADecrypt(c, PrivateKey):
    '''
    Input:
        c: Encrypted message
        PrivateKey: A tuple (N, d)
    Returns:
        m: Decrypted message
    '''
    N = PrivateKey[0]
    d = PrivateKey[1]
    return fastPowering(c, d, N)


