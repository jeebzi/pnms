#include "param.h"

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
	return 0;
}
	
