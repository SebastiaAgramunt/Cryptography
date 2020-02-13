from random import randrange
from typing import List, Tuple


# A handy list of small primes
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

def xgcd(a: int, b: int) -> int:
    """
    Extended greatest common divisor.
    ttps://en.wikipedia.org/wiki/Extended_Euclidean_algorithm

    Function that solves the equation au+bv=gcd(a,b) for a, b, 
    u, v integers.

    Input:
        a: an integer
        b: an integer

    Output:
        g: g(a, b) greatest common divisor of a and b
        u: the factor multiplying a in eq above
        v: the factor multiplying b in eq above

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

def PrimesSieveEratosthenes(limit: int):
    """
    Calculate list of primes until an integer "limit"
    https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes

    The function is a generator, so you have to cast it
    to a list if you want a full list of all the primes.

    Input:
        limit: an max integer up to where we want to
               calculate prime numbers
    Output:
        A generator of prime numbers up until limit.

    e.g. Calculate primes up to 100
    n = 100
    primes = list(primes_sieve_eratosthenes(n))
    print("All prime numbers bellow {} are:\n{}".format(n,primes))

    e.g.2 (bits) Calculate all prime numbers smaller than 2^16
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

def InvertibleNumbers(m: int) -> List[int]:
    """
    Given a natural number m it returns all the invertible numbers modulo(m) 
    for the product operation.

    Input:
        m: An integer

    Output:
        inv: A list of all the invertible numbers

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

def InverseMod(a: int, m: int) -> int:
    """
    Given natural number a and m returns b = a^{-1}(mod m). The inverse modulo m of a.
    This is b*a (mod p )=1

    Input:
        a: an integer element of the field (1 < a < m)
        m: an integer
    Output:
        b: The inverse modulo m of a

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

def fastPowering(a: int, k: int, m: int) -> int:
    """
    Fast powering algorithm in Zm. Calculates
    b = a^{k} (mod m)

    Input:
        a: the base (1 < a < m)
        k: the exponent
        m: the moudlar factor

    Output:
        b: a^k (mod m)
    
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

def InverseFermat(a: int, p: int) -> int:
    """
    Given a natural number "a" over a field generated by "p", the function
    generates "ai" the inverse of "a", this is ai*a(mod p)=1.

    ALGORITHM: Using Fermat's little theorem. If p is prime then a^{p-1}(mod p)=1
    so the inverse of a is got by multiplying both sides by a^{-1}, therefore:
    a^{-1}=a^{p-2} (mod p). The inverse is then a^{p-2}(mod p)

    Input: 
        a: an integer number smaller than p
        p: a prime number

    Output:
        a^{-1}: The inverse of a

    e.g.
    p = 17
    a = 10
    print("Inverse modulo of {} with p={} is {}. a*a^-1={}".format(a, p, InverseFermat(a,p), a*InverseFermat(a,p)%p))
    """
    # p has to be a prime number
    if not isPrime(p, 40):
        raise ValueError(f"{p} is not prime number, cannot calculate inverse using Fermat's theorem")

    if a>p:
        raise ValueError(f"a={a} cannot be larger than p={p}, the element is not in the group, try input a%p instead")
    
    return fastPowering(a, p-2, p)

def EulerTotient(m: int) -> int:
    """
    Given a number m, gives the number of generators of the group, i.e. 
    the nubmer of elements that can generate all the group by powering

    Input:
        m: an integer that generates a multiplicative group Gm
    Output:
        c: the Euler totient

    Observation: If m is prime the Euler totient is m-1
    """
    c = 0
    for i in range(m):
        g, _, _ = xgcd(i, m)
        if g == 1:
            c+=1
    return c 

def isPrime(n: int, m: int =20) -> bool:
    """
    Check primality of a number n using Miller-Rabin algorithm:
    Book An introduction to mathemahical cryptography (Hoffstein, 
    Pipher and Silverman) page 127.

    n: The number we want to check primality
    m: Number of checks to perform

    Returns: False if the number is composite with 100% certainty
             True if the number is probable to be prime (the larger the m, 
             the more sure we are)
    """
    # Check if the number is in among the list of small primes
    if n in SMALL_PRIMES:
        return True

    if n < 2:
        return False

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

def RandomPrime(size: int=256, m: int=40) -> int:
    """
    Returns a random prime number of size bits
    Use Miller-Rabin algorithm to check primality

    Inputs:
        size: the size in bits of the desired prime
        m: number of times to run Miller-Rabin primality
           check.
    Outputs:
        p: a random prime of the expected size

    """
    while True:
        p = randrange(2**(size-1), 2**(size))
        if isPrime(p, m):
            return p

def _g(x: int, c: int, n: int) -> int:
    # Function x^2+c mod n
    return (fastPowering(x, 2, n) + c)%n

def PollardsRhoFactorisation(n: int) -> int:
    """
    Finds a non trivial factor for number n, 
    found factor is not necessarily a prime number(check!).

    ALGORITHM: Pollard's rho factorisation from book
    An introduction to mathematical cryptography by Hoffstein,
    Pipher and Silverman page 133

    Inputs:
        n: a natural number
    Outputs:
        d: a factor of n. i.e. n%d==0
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

def BruteForceFactorisation(n: int, primes: List[int]=SMALL_PRIMES):
    """
    Finds the list of factors of n using a set of prime numbers. And also
    the number of times each number has to be multiplied.

    Input;
        n: a composite integer (not prime)
        primes: a list of candiate primes to factorise n

    Output:
        a tuple, (factors, reps)
        factors: the prime numbers that factor n
        reps: the repetitions of those prime numbers

    ##calculate back n:
    n = 1
    for factor, rep in zip(factors, reps):
        n*=factor**rep
    print(f"Recovered n from factors and reps {n}")

    WARNING: if n is too large compared to the list of primes we may not be
    able to factor n. You can then increase the list of primes by inputing
    as parameter list(primes_sieve_eratosthenes(1<<16)).
    We have assured convergence if primes[-1]>n

    """
    if isPrime(n,40):
        return [n], [1]
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

def PrimeFieldGenerator(p: int, primes: List[int]=SMALL_PRIMES):
    '''
    Gives the generator of a prime field defined by p

    Input:
        p: prime number
        primes: a list of prime numbers to be used to factorise

    Output:
        g: a generator of the field

    ALGORITHM IN BOOK: A Computational Introduction to Number Theory
    and Algebra by Shoup. Page 269.
    WARNING: If primes list is small or p is very large we may not be able
    to find the generator, since this method requires factorisation of the
    order = p-1.
    '''
    if not isPrime(p, 100):
        raise ValueError(f"{p} is not prime, please select a prime number")

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

def _FindSafePrimes(size: int=256, krange: int=500):
    """
    Function that finds two primes p and q related by
    p = kq +1, so order(p)=p-1=kq
    it is used to factor p-1 number

    input:
        size: size of the prime numbers
    output: 
        the triplet (p, q, k) following the equation
        p = kq + 1
            p: prime number
            q: prime number
            k: factorisation
    """
    # Start with non prime numbers
    while True:
        p = RandomPrime(size, 100)
        for k in range(2, krange):
            # p-1=kq, check if p-1 is divisible by q
            q = (p-1)//k
            if (p-1)%k==0 and isPrime(q, 100):
                # if divisible, assign q = (p-1)/k
                # we've found the values
                return p, q, k


def GeneratePrimeGeneratorPair(size: int=256):
    """
    Generates a random prime number p and a
    random generator g for the field Fp.

    Input: 
        size: Size in bits of the random prime p
    Output:
        tuple (p, g)
        p is prime number of specified size
        g is generator for the field Fp

    """

    #find primes p, q such that p = kq + 1
    p, q, k = _FindSafePrimes(size)
    assert (p-1)==k*q, "FindSafePrimes is not well implemented"
    #p, q, k = 38327, 19163, 2

    order = p-1
    # order of p is p-1. Has factorisation k*q
    # but k is composed (smaller than krange in FindSafePrimes,
    # we have to find its prime factors 
    factors, reps = BruteForceFactorisation(k)

    #add q as a factor of p-1
    factors+=[q]
    reps+=[1]

    # A value 1 < g < p will be a generator to a prime p iff
    # for every prime factor q of p-1 (the order) we have
    # b = g^((p-1)/q) != 1 (mod p)
    #
    # ALGORITHM IN BOOK: A Computational Introduction to Number Theory
    # and Algebra by Shoup. Page 269.
    while True:
        g = randrange(1, p)
        Found = True
        for factor, rep in zip(factors, reps):
            b = fastPowering(g, order*InverseFermat(factor, p)%p, p)
            if b==1:
                Found = False
                break
        if Found:
            return p, g


def isGeneratorBruteForce(g: int, p: int):
    l = len(set([pow(g, j, p) for j in range(p)]))
    return True if l == p-1 else False

def LCM(a: int, b: int):
    '''
    Computes the least common multiple of two integers
    Input:
        a: integer
        b: integer
    Output:
        least common multiple
    '''
    g, _, _ = xgcd(a, b)

    return a*b//g


def RSAKeyGenerator(size: int=16):
    '''
    RSA key generation. Generates public and
    private keys in RSA protocol

    Input:
        size: size in bits of the field

    Output:
        PublicKey: (N, e)
        PrivateKey: (N, d)

    '''

    # Generate two random primes of n bits
    p = RandomPrime(size, m=40)
    q = RandomPrime(size, m=40)

    # p and q must be different primes
    while p==q:
      q = RandomPrime(size, m=40)  

    N = p*q
    lm = LCM(p-1, q-1)
    
    while True:
        e = randrange(2, lm)
        g, _, _ = xgcd(e, lm)
        if g==1:
            d = InverseMod(e, lm)
            # return public and private keys
            return (N, e), (N, d)

def RSAEncrypt(m: int, PublicKey: Tuple[int]):
    '''
    Encrypts a message m using the RSA public key

    Input:
        m: message (An integer message)
        PublicKey: A tuple (N, e)
    Returns:
        c: Encrypted message
    '''
    N = PublicKey[0]
    e = PublicKey[1]

    return fastPowering(m, e, N)

def RSADecrypt(c: int, PrivateKey: Tuple[int]):
    '''
    Decrcypts the ciphertext c using the private key
    Input:
        c: Encrypted message
        PrivateKey: A tuple (N, d)
    Returns:
        m: Decrypted message
    '''
    N = PrivateKey[0]
    d = PrivateKey[1]
    return fastPowering(c, d, N)


