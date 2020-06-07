from crypt import isPrime

from random import randrange, shuffle, seed
from typing import List

# small prime number
field = 67

# bits for the representation of integers
l = 64

# number of parties
parties = 2

# seed
sem = 1

#
def to_binl(x: int, l: int) -> str:
    """
    Calculates the coeficients living on a field Zp for a polynomial of degree n
    Input:
        x: integer to encode (smaller than 2**l)
        l: number of bits to encode
    Output:
        x expressed in l bits
    """
    assert x<2**l, f"{x} has to be smaller than {2**l} and therefore cannot be represented as int64"
    
    b = bin(x)
    if len(b)==l+2:
        return b
    else:
        n = l+2-len(b)
        return '0b' + '0'*n + b.split("b")[1]


def share_bin(x: int, l: int, field: int, parties: int = 2) -> List[List[int]]:
    """
    Given an integer x, calculate the shares for n parties in the field.
    Input:
        x: integer (smaller than 2**l)
        l: number of bits to encode
        field: The field of the additive sharing
        parties: number of parties for secret sharing
    Output:
        A list of a list of integers with max range field. Each list is
        the share for each party.
    """

    assert isPrime(field, 100), f"{field} is not prime"

    x_bin = to_binl(x, l)
    shares = [[] for _ in range(parties)]

    for bit in x_bin[2:]:
        rdm = []
        for i in range(parties-1):
            rdm.append(randrange(field))
            shares[i].append(rdm[-1])

        shares[parties-1].append((int(bit)-sum(rdm))%field)

    return shares


def reconstruct_bin(shares: List[List[int]], field: int) -> str:
    """
    From the shares reconstruct the binary number
    Input:
        x: integer (smaller than 2**l)
        l: number of bits to encode
        field: The field of the additive sharing
        parties: number of parties for secret sharing
    Output:
        A list of a list of integers with max range field. Each list is
        the share for each party.
    """

    assert isPrime(field, 100), f"{field} is not prime"

    x = []
    l = len(shares[0])

    bit = 0
    while bit < l:
        s = 0
        for share in shares:
            s+=share[bit]

        x.append(s%field)
        bit+=1

    return '0b' + ''.join(map(str, x))


# def calc_w_c(x: List[int], r: List[int], p: int, party: int, beta: bool):

#     assert len(x)==len(r), "x and r have to be same lenth l"
#     l = len(x)

#     w, c = [], []
#     for i in range(l-1, -1, -1):
#         if beta==0:
            
#             for x_i, r_i in zip(x, r):
#                 w.append((x_i + party*r_i -2*r_i*x_i)%p)

#                 k = i+1
#                 tmp = 0
#                 for k in range(i+1, l)
#                 c.append(party*r_i -x_i + party + )






if __name__ == "__main__":

    print("Initial parameters\n====================")
    print(f"Field={field}\nl_bits={l}\nparties={parties}\nseed={sem}")
    seed(sem)


    x = randrange(0, 2**l)
    print(f"\nPick a random value x:\n\tx={x} < {2**l}\n")
    print(f"\tx_bin={to_binl(x, l)}")

    shares = share_bin(x, l, field, parties=parties)

    print()
    for i, share in enumerate(shares):
        print(f"Party {i}:\n\t({', '.join(map(str, share))})")

    reconstruct = reconstruct_bin(shares, field) 

    

    print(f"\n\nReconstruction of x is: \n\t{reconstruct}")
    print(f"The original value of x is: \n\t{to_binl(x, l)}")

    # checking if reconstruction is done ok
    assert reconstruct == to_binl(x, l)
    print("which is the same!")



    # calculating random common parameters
    r = randrange(0, 2**l)
    r_bin = to_binl(r, l)
    t = (r+1)%2**l

    b = randrange(0, 2)

    # s and u are random in F*p
    s = [randrange(1, field) for _ in range(l)]
    u = [randrange(1, field) for _ in range(l)]


    # generate the permutation
    perm = list(range(l))
    shuffle(perm)


    print(f"\n\nr is:\n\t{r}\n\t{r_bin}")
    print(f"t is:\n\t{t}")
    print(f"beta is:\n\t{b}")
    print(f"s is:\n\t{s}")
    print(f"u is:\n\t{u}")
    print(f"permutation is:\n\t{perm}")
    



    
