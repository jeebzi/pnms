from sage.all import *
import random as r


def squaremult(a, b, m):
    c = 1
    while b:
        if b & 1:
            c = c * a % m
        a = a * a % m
        b = b >> 1
    return c % m

def make_matrice_gamma(gamma, n, p):
    """
    fait la lsite des polynomes que s'annule en gamma
    """
    matrice = []
    i = 0
    while i < n:
        j = 0
        tab = []
        while j < n:
            if i == 0 and j == 0:
                tab += [p]
            elif i != 0 and j == 0:
                tab += [-gamma**i]
            elif i == j:
                tab += [1]
            else:
                tab += [0]
            j += 1
        matrice += [tab]
        i += 1
    return matrice


def gen_pnms():

    p = random_prime(2**256)
    n = 5
    lambd = 2
    phi = 2**64
    K = GF(p)
    pol = PolynomialRing(K, "X")
    X = pol("X")
    found_lambda = False
    while not found_lambda:
        # E = ZZ["X"](list(X**n-lambd))
        E = ZZ["X"](f"X^{n} - {lambd}")
        d, v, u = xgcd((squaremult(X, p, E) - X) % E, E)

        # gamma = -d[0] % p
        gamma = d.roots()[0][0]
        B = matrix(ZZ, make_matrice_gamma(gamma, n, p)).LLL()
        if 2 * n * abs(lambd) * B.norm(1) < phi:
            found_lambda = True
        elif lambd > phi / 2 * n:
            lambd = 2
            n += 1
        lambd += 2
    rho = B.norm(1) - 1
    #M
    i = 0
    while B[i][0] % 2 == 0:
        i += 1
    M = B[i]
    M = ZZ["X"](list(M))
    #inverse M
    d, Mu, Ev = xgcd(M, E)
    d = int(d)
    d_inv = pow(d, -1, phi)
    M_inv = d_inv * Mu % phi

    pnms = [p, n, gamma, rho, E]
    return pnms, M, M_inv, phi, lambd

def mult_montg_pnms_rand(pnms, M, M_inv, phi):
    """
    prend en parametre un pnms génère 2 poly aléatoire
    puis effetctue la mutliplication de montgomery sur ces dernier
    pnms: [p, n, gamma, rho, E]
    """
    p, n, gamma, rho, E = pnms
    A = gen_poly_random(n - 1, int(-rho) + 1, int(rho))
    B = gen_poly_random(n - 1, int(-rho) + 1, int(rho))
    A = ZZ["X"](A)
    B = ZZ["X"](B)

    C = A * B % E
    Q = (C * M_inv % E) % phi
    C_prime = (C - (Q * M % E)) / phi
    C_prime = ZZ["X"](C_prime)
    print("C\'", C_prime)

    return (C_prime(gamma) == A(gamma) * B(gamma) * pow(phi, -1, p) % p)


def gen_poly_random(n, inf, sup):
    """
    génère un polynôme amélatoire de degrè n avec ces facteurs compris entre inf et sup
    """
    i = 0
    pol = []
    while i <= n:
        pol += [r.randint(inf, sup)]
        i += 1
    return pol
    
if __name__ == "__main__":
    pnms, M, M_inv, phi, lambd = gen_pnms()
    print(mult_montg_pnms_rand(pnms, M, M_inv, phi))
