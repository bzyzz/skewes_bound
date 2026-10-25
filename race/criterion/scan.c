// E(u) = c0*Phi(u) + sum_j Re(A_j e^{i g_j u}) + (1/u) sum_j Re(B_j e^{i g_j u}) on a uniform grid
// by phase rotation; each block of `blk` points is re-anchored with cexp.
// Output: per-block min and argmin; optionally the full grid.
// gcc -O3 -march=native -fopenmp -shared -fPIC scan.c -o libscan.so -lm
#include <math.h>
#include <stdlib.h>
#include <complex.h>

void scan(const double *g, const double *Ar, const double *Ai, const double *Br, const double *Bi,
          long nz, double c0, double u0, double h, long npts, long blk,
          double *bmin, double *bu, double *full) {
    long nblk = (npts + blk - 1) / blk;
    #pragma omp parallel
    {
        double *sa = malloc(sizeof(double) * blk);
        double *sb = malloc(sizeof(double) * blk);
        double complex *w = malloc(sizeof(double complex) * nz);
        for (long j = 0; j < nz; j++) w[j] = cexp(I * g[j] * h);
        #pragma omp for schedule(dynamic)
        for (long b = 0; b < nblk; b++) {
            long i0 = b * blk, i1 = i0 + blk; if (i1 > npts) i1 = npts;
            long m = i1 - i0;
            double ub = u0 + h * i0;
            for (long k = 0; k < m; k++) { sa[k] = 0; sb[k] = 0; }
            for (long j = 0; j < nz; j++) {
                double complex zz = cexp(I * g[j] * ub), wj = w[j];
                double ar = Ar[j], ai = Ai[j], br = Br[j], bi = Bi[j];
                for (long k = 0; k < m; k++) {
                    double zr = creal(zz), zi = cimag(zz);
                    sa[k] += ar * zr - ai * zi;
                    sb[k] += br * zr - bi * zi;
                    zz *= wj;
                }
            }
            double mn = 1e300, mu = 0;
            for (long k = 0; k < m; k++) {
                double u = u0 + h * (i0 + k);
                double phi = 1 + 2 / u + 8 / (u * u) + 48 / (u * u * u);
                double e = c0 * phi + sa[k] + sb[k] / u;
                if (full) full[i0 + k] = e;
                if (e < mn) { mn = e; mu = u; }
            }
            bmin[b] = mn; bu[b] = mu;
        }
        free(sa); free(sb); free(w);
    }
}
