// First lead changes in prime races, by direct counting (segmented sieve of Eratosthenes).
// For every modulus 3 <= q <= QMAX with q != 2 mod 4 and every ordered pair of units (a, b), report the first prime x
// with pi(x;q,a) > pi(x;q,b) (0 if none <= X). For every q also report the first x with
// avg_{residues} pi(x;q,.) > avg_{non-residues} pi(x;q,.) (team race), residues = squares of units.
// usage: count_races X [QMAX]           cc -O3 -march=native count_races.c -o count_races -lm
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define QM 64
static int qs[QM], nq;
static int units[QM][QM], nunits[QM], isres[QM][QM], isunit[QM][QM];
static int64_t cnt[QM][QM];                 // cnt[iq][a]
static int64_t first[QM][QM][QM];           // first[iq][a][b]
static int pend[QM][QM][QM], npend[QM][QM]; // pending b's for (iq, a)
static int64_t sumR[QM], sumN[QM], nR[QM], nN[QM], team_first[QM];

static int gcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }

static inline void add_prime(int64_t p) {
    for (int i = 0; i < nq; i++) {
        int q = qs[i];
        int a = (int)(p % q);
        if (!isunit[i][a]) continue;
        int64_t ca = ++cnt[i][a];
        if (isres[i][a]) sumR[i]++; else sumN[i]++;
        int *pl = pend[i][a], n = npend[i][a];
        for (int k = 0; k < n; ) {
            int b = pl[k];
            if (ca > cnt[i][b]) { first[i][a][b] = p; pl[k] = pl[--n]; } else k++;
        }
        npend[i][a] = n;
        if (!team_first[i] && nN[i] * sumR[i] > nR[i] * sumN[i]) team_first[i] = p;
    }
}

int main(int argc, char **argv) {
    int64_t X = (int64_t)atof(argv[1]);
    int QMAX = argc > 2 ? atoi(argv[2]) : 30;
    for (int q = 3; q <= QMAX; q++) {
        if (q % 4 == 2) continue;
        int i = nq++; qs[i] = q;
        for (int a = 1; a < q; a++) if (gcd(a, q) == 1) { units[i][nunits[i]++] = a; isunit[i][a] = 1; }
        for (int k = 0; k < nunits[i]; k++) { int y = units[i][k]; isres[i][(y * y) % q] = 1; }
        for (int k = 0; k < nunits[i]; k++) { if (isres[i][units[i][k]]) nR[i]++; else nN[i]++; }
        for (int k = 0; k < nunits[i]; k++) {
            int a = units[i][k];
            for (int l = 0; l < nunits[i]; l++) if (l != k) pend[i][a][npend[i][a]++] = units[i][l];
        }
    }
    // base primes up to sqrt(X)
    int64_t R = (int64_t)sqrt((double)X) + 2;
    char *bs = calloc(R + 1, 1);
    int64_t *bp = malloc(sizeof(int64_t) * (R / 2 + 10)); int64_t nbp = 0;
    for (int64_t i = 2; i <= R; i++) if (!bs[i]) { bp[nbp++] = i; for (int64_t j = i * i; j <= R; j += i) bs[j] = 1; }
    // segmented sieve over odd numbers; prime 2 is not a unit for even q and is counted for odd q
    add_prime(2);
    const int64_t S = 1 << 21;                       // odd numbers per segment
    char *seg = malloc(S);
    int64_t *nxt = malloc(sizeof(int64_t) * nbp);   // next odd multiple index for each odd base prime
    for (int64_t k = 1; k < nbp; k++) { int64_t p = bp[k]; nxt[k] = (p * p - 1) / 2; } // index of p^2 (odd n = 2i+1)
    int64_t imax = (X - 1) / 2;                       // largest index with 2i+1 <= X
    int64_t next_report = 1000000;
    for (int64_t lo = 1; lo <= imax; lo += S) {        // index 0 is n = 1 (not prime)
        int64_t hi = lo + S; if (hi > imax + 1) hi = imax + 1;
        memset(seg, 0, hi - lo);
        for (int64_t k = 1; k < nbp; k++) {
            int64_t p = bp[k], j = nxt[k];
            if (j >= hi) continue;
            if (j < lo) j = lo + ((p - (lo - j) % p) % p);
            for (; j < hi; j += p) seg[j - lo] = 1;
            nxt[k] = j;
        }
        for (int64_t j = lo; j < hi; j++) if (!seg[j - lo]) add_prime(2 * j + 1);
        if (2 * hi + 1 >= next_report) {
            fprintf(stderr, "reached %.3e\n", (double)(2 * hi + 1)); next_report *= 10;
        }
    }
    printf("# X=%lld\n", (long long)X);
    for (int i = 0; i < nq; i++) {
        int q = qs[i];
        printf("team %d first %lld nR %lld nN %lld\n", q, (long long)team_first[i], (long long)nR[i], (long long)nN[i]);
        for (int k = 0; k < nunits[i]; k++) for (int l = 0; l < nunits[i]; l++) if (k != l) {
            int a = units[i][k], b = units[i][l];
            printf("pair %d %d %d first %lld res %d %d\n", q, a, b, (long long)first[i][a][b], isres[i][a], isres[i][b]);
        }
    }
    return 0;
}
