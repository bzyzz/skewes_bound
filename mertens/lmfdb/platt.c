// Decoder for D. J. Platt's zeta-zero files as distributed by the LMFDB.
// File: uint64 number_of_blocks; then blocks: double t0, double t1, uint64 Nt0, uint64 Nt1, followed by
// (Nt1 - Nt0) records of 13 bytes (uint64 z1, uint32 z2, uint8 z3). Within a block Z accumulates
// Z += z1 + z2*2^64 + z3*2^96, and the zero is gamma = t0 + Z * 2^-101 (stored to +-2^-102).
// gcc -O2 -shared -fPIC platt.c -o libplatt.so
#include <stdio.h>
#include <stdint.h>
#include <string.h>

// Count zeros and blocks in a file. Returns number of zeros, or -1 on error.
long platt_count(const char *path, long *nblocks_out) {
    FILE *f = fopen(path, "rb"); if (!f) return -1;
    uint64_t nb; if (fread(&nb, 8, 1, f) != 1) { fclose(f); return -1; }
    long total = 0;
    for (uint64_t b = 0; b < nb; b++) {
        double t0, t1; uint64_t N0, N1;
        if (fread(&t0, 8, 1, f) != 1 || fread(&t1, 8, 1, f) != 1 || fread(&N0, 8, 1, f) != 1 || fread(&N1, 8, 1, f) != 1) { fclose(f); return -1; }
        long k = (long)(N1 - N0);
        total += k;
        if (fseek(f, 13L * k, SEEK_CUR)) { fclose(f); return -1; }
    }
    *nblocks_out = (long)nb; fclose(f); return total;
}

// Decode all zeros. Per zero: t0 of its block, Z = hi*2^64 + lo, global index N (1-based: N = Nt0 + j + 1).
// Per block: t0, t1, Nt0, Nt1. Returns number of zeros written, or -1 on error / overflow.
long platt_decode(const char *path, double *t0z, uint64_t *hi, uint64_t *lo, uint64_t *idx, long cap,
                  double *bt0, double *bt1, uint64_t *bN0, uint64_t *bN1, long bcap) {
    FILE *f = fopen(path, "rb"); if (!f) return -1;
    uint64_t nb; if (fread(&nb, 8, 1, f) != 1) { fclose(f); return -1; }
    if ((long)nb > bcap) { fclose(f); return -1; }
    long n = 0;
    unsigned char rec[13];
    for (uint64_t b = 0; b < nb; b++) {
        double t0, t1; uint64_t N0, N1;
        if (fread(&t0, 8, 1, f) != 1 || fread(&t1, 8, 1, f) != 1 || fread(&N0, 8, 1, f) != 1 || fread(&N1, 8, 1, f) != 1) { fclose(f); return -1; }
        bt0[b] = t0; bt1[b] = t1; bN0[b] = N0; bN1[b] = N1;
        unsigned __int128 Z = 0;
        for (uint64_t j = N0; j < N1; j++) {
            if (fread(rec, 13, 1, f) != 1 || n >= cap) { fclose(f); return -1; }
            uint64_t z1; uint32_t z2; uint8_t z3;
            memcpy(&z1, rec, 8); memcpy(&z2, rec + 8, 4); z3 = rec[12];
            Z += (unsigned __int128)z1 + ((unsigned __int128)z2 << 64) + ((unsigned __int128)z3 << 96);
            t0z[n] = t0; hi[n] = (uint64_t)(Z >> 64); lo[n] = (uint64_t)Z; idx[n] = j + 1;
            n++;
        }
    }
    fclose(f); return n;
}
