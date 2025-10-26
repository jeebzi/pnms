#include "small_param.h"
int degree_trim(const int64_t *p, int len) {
    int i;
    for (i = len - 1; i >= 0; --i) if (p[i] != 0) return i;
    return 0;
}

/* addition: res = a + b (taille = n) */
void poly_add(int64_t *res, const int64_t *a, const int64_t *b, int n) {
    int i;
    for (i = 0; i < n; ++i) res[i] = a[i] + b[i];
}

/* soustraction: res = a - b (taille = n) */
void poly_sub(int64_t *res, const int64_t *a, const int64_t *b, int n) {
    int i;
    for (i = 0; i < n; ++i) res[i] = a[i] - b[i];
}

/* multiplication brute: a,b de degré < n  => prod de taille 2n-1 */
int64_t *poly_mul_raw(const int64_t *a, const int64_t *b, int n) {
    int outlen = 2 * n - 1;
    int64_t *out = calloc((size_t)outlen, sizeof(int64_t));
    if (!out) return NULL;
    int i, j;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            out[i + j] += a[i] * b[j];
        }
    }
    return out;
}

/* Division euclidienne: remainder = A mod B
   A: array a of length la (degA = la-1)
   B: array b of length lb (degB = lb-1), on suppose b[degB] != 0
   Retourne un tableau remainder de longueur lb (caller free).
*/
int64_t *poly_mod(const int64_t *a, int la, const int64_t *b, int lb) {
    int64_t *R = calloc((size_t)la, sizeof(int64_t));
    if (!R) return NULL;
    int i;
    for (i = 0; i < la; ++i) R[i] = a[i];

    int degA = la - 1;
    int degB = lb - 1;
    if (degB < 0) { free(R); return NULL; }

    for (i = degA; i >= degB; --i) {
        int64_t rcoeff = R[i];
        if (rcoeff == 0) continue;
        int64_t blead = b[degB];
        if (blead == 0) { free(R); return NULL; }
        int64_t factor = rcoeff / blead;
        int j;
        for (j = 0; j <= degB; ++j) {
            R[i - degB + j] -= factor * b[j];
        }
    }

    int64_t *rem = calloc((size_t)lb, sizeof(int64_t));
    if (!rem) { free(R); return NULL; }
    for (i = 0; i < lb; ++i) rem[i] = R[i];
    free(R);
    return rem;
}

/* réduction par E: prend poly de longueur up to 2n-1 et renvoie remainder de longueur n (deg < n) */
int64_t *reduce_mod_E(const int64_t *poly, int len, const int64_t *E, int degE, int n) {
    int lb = degE + 1;
    int la = len;
    int64_t *rem_full = poly_mod(poly, la, E, lb);
    if (!rem_full) return NULL;
    int64_t *rem = calloc((size_t)n, sizeof(int64_t));
    if (!rem) { free(rem_full); return NULL; }
    int copy = (lb < n) ? lb : n;
    int i;
    for (i = 0; i < copy; ++i) rem[i] = rem_full[i];
    free(rem_full);
    return rem;
}

/* coefficient-wise modulo phi (donne valeurs dans [0, phi-1]) */
void poly_coeffs_mod_phi(int64_t *out, const int64_t *in, int n, uint64_t phi) {
    int i;
    for (i = 0; i < n; ++i) {
        int64_t v = in[i] % (int64_t)phi;
        if (v < 0) v += (int64_t)phi;
        out[i] = v;
    }
}

/* division par phi (vérifie la divisibilité exacte des coefficients) */
int poly_divide_by_phi(int64_t *out, const int64_t *in, int n, uint64_t phi) {
    int i;
    for (i = 0; i < n; ++i) {
        if (in[i] % (int64_t)phi != 0) return -1;
        out[i] = in[i] / (int64_t)phi;
    }
    return 0;
}

/* produit mod E: C = (A*B) mod E ; A,B length n ; renvoie tableau length n (malloc) */
int64_t *poly_mul_mod(const int64_t *A, const int64_t *B, int n, const int64_t *E, int degE) {
    int64_t *prod = poly_mul_raw(A, B, n);
    if (!prod) return NULL;
    int prollen = 2 * n - 1;
    int64_t *rem = reduce_mod_E(prod, prollen, E, degE, n);
    free(prod);
    return rem;
}

/* multiplication de Montgomery version PNMS */
int64_t *montgomery_mul_pnms(const int64_t *A, const int64_t *B,
                             int n,
                             const int64_t *E, int degE,
                             const int64_t *M, const int64_t *Minv,
                             uint64_t phi)
{
    int64_t *C = poly_mul_mod(A, B, n, E, degE);
    if (!C) return NULL;

    int64_t *t1 = poly_mul_mod(C, Minv, n, E, degE);
    if (!t1) { free(C); return NULL; }

    int64_t *Q = calloc((size_t)n, sizeof(int64_t));
    if (!Q) { free(C); free(t1); return NULL; }
    poly_coeffs_mod_phi(Q, t1, n, phi);
    free(t1);

    int64_t *S = poly_mul_mod(Q, M, n, E, degE);
    if (!S) { free(C); free(Q); return NULL; }

    int64_t *U = calloc((size_t)n, sizeof(int64_t));
    if (!U) { free(C); free(Q); free(S); return NULL; }
    int i;
    for (i = 0; i < n; ++i) U[i] = C[i] - S[i];

    int64_t *Cprime = calloc((size_t)n, sizeof(int64_t));
    if (!Cprime) { free(C); free(Q); free(S); free(U); return NULL; }
    int ok = poly_divide_by_phi(Cprime, U, n, phi);
    if (ok != 0) {
        free(C); free(Q); free(S); free(U); free(Cprime);
        return NULL;
    }

    int64_t *final = reduce_mod_E(Cprime, n, E, degE, n);
    free(C); free(Q); free(S); free(U); free(Cprime);
    return final;
}


int64_t rand_int64(int64_t inf, int64_t sup) {
	uint64_t range = (uint64_t)sup - (uint64_t)inf + 1;
	uint64_t limit = UINT64_MAX - (UINT64_MAX % range);
	uint64_t r;
	do {
		r = (uint64_t)rand();
	} while (r >= limit);

	return (int64_t)(r % range) + inf;
}

int64_t * gen_random_poly(uint64_t n, int64_t inf, int64_t sup) {
	/* genère aléatoirement un polynôme
	 * de degrès n avec ces coefs entre inf et sup
	 */
	int64_t *poly;
	poly = (int64_t*)malloc((n + 1) * sizeof(int64_t));

	int i = 0;
	while (i <= n) {
		poly[i] = rand_int64(inf, sup);
		i += 1;
	}
	return poly;
}

void print_poly(int64_t *poly, int deg) {
	int i = 0;
	printf("[");
	while (i <= deg) {
		printf("%ld", poly[i]);
		if (i < deg) printf(", ");
		i += 1;
	}
	printf("]\n");
}

int main() {
	srand(time(NULL));
	int64_t *A = gen_random_poly(n - 1, -rho + 1, rho);
	int64_t *B = gen_random_poly(n - 1, -rho + 1, rho);
	int64_t *C_prime = montgomery_mul_pnms(A, B, n, E, n, M, M_inv, phi);
	print_poly(A, n - 1);
	print_poly(B, n - 1);
	print_poly(C_prime, n - 1);
	return 0;
}
	
