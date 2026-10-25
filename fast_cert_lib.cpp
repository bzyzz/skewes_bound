#include "fast_cert_zeros.hpp"
#include <complex>
#include <cstring>

static inline long double kappa_2m_hat_ld(long double xi, int m) {
    constexpr long double PI = 3.14159265358979323846264338327950288L;
    if (fabsl(xi) < 1e-12L) return 1.0L;
    long double xi2 = xi * xi;
    for (int k = 1; k <= m; ++k) {
        long double kpi = k * PI;
        if (fabsl(fabsl(xi) - kpi) < 1e-7L) {
            long double u = (fabsl(xi) - kpi);
            long double sinc_u = (fabsl(u) < 1e-12L) ? 1.0L : (sinl(u) / u);
            long double sign = (k & 1) ? 1.0L : -1.0L;
            long double res = sign * sinc_u * (kpi * kpi) / (xi * (kpi + fabsl(xi)));
            for (int j = 1; j <= m; ++j) {
                if (j != k) {
                    long double jp2 = (j * PI) * (j * PI);
                    res *= jp2 / (jp2 - xi2);
                }
            }
            return res;
        }
    }
    long double res = sinl(xi) / xi;
    for (int k = 1; k <= m; ++k) {
        long double kp2 = (k * PI) * (k * PI);
        res *= kp2 / (kp2 - xi2);
    }
    return res;
}

extern "C" {

void certify_chunks_omp(
    const double *T_edges,
    const int *K_expected,
    int N_chunks,
    double u_ref,
    const char *out_dir,
    int *status_out,
    double *max_rad_out,
    double *bin_gamma_out,
    double *Z0_re, double *Z0_im,
    double *Z1_re, double *Z1_im,
    double *Z2_re, double *Z2_im,
    int n_threads
) {
    double T_max = T_edges[N_chunks];
    int max_M = (int)std::floor(std::sqrt(T_max / (2.0 * M_PI))) + 4;
    std::vector<long double> log_n_ld(max_M + 1);
    std::vector<double> log_n(max_M + 1), inv_sqrt_n(max_M + 1);
    for (int n = 1; n <= max_M; ++n) {
        log_n_ld[n] = std::log((long double)n);
        log_n[n] = (double)log_n_ld[n];
        inv_sqrt_n[n] = 1.0 / std::sqrt((double)n);
    }

    constexpr long double TWO_PI_LD = 6.28318530717958647692528676655900577L;

    #pragma omp parallel for schedule(dynamic, 1) num_threads(n_threads)
    for (int c = 0; c < N_chunks; ++c) {
        double Ta = T_edges[c];
        double Tb = T_edges[c + 1];
        double chunk_len = Tb - Ta;
        int K_exp = K_expected[c];
        double gb = 0.5 * (Ta + Tb);
        bin_gamma_out[c] = gb;

        char fpath[512];
        bool have_cached = false;
        if (out_dir && out_dir[0] != '\0') {
            snprintf(fpath, sizeof(fpath), "%s/chunk_%08lld_%08lld.bin",
                     out_dir, (long long)std::llround(Ta), (long long)std::llround(Tb));
            FILE *fc = fopen(fpath, "rb");
            if (fc) {
                int32_t k_hdr = 0;
                double ta_hdr = 0.0;
                if (fread(&k_hdr, sizeof(int32_t), 1, fc) == 1 &&
                    fread(&ta_hdr, sizeof(double), 1, fc) == 1 &&
                    k_hdr == K_exp && std::abs(ta_hdr - Ta) < 1e-6) {
                    std::vector<double> omids(K_exp);
                    std::vector<float> rads(K_exp);
                    bool read_ok = true;
                    double mrad = 0.0;
                    for (int i = 0; i < K_exp; ++i) {
                        if (fread(&omids[i], sizeof(double), 1, fc) != 1 ||
                            fread(&rads[i], sizeof(float), 1, fc) != 1) {
                            read_ok = false; break;
                        }
                        if (rads[i] > mrad) mrad = rads[i];
                    }
                    fclose(fc);
                    if (read_ok) {
                        status_out[c] = 1;
                        max_rad_out[c] = mrad;
                        double z0r = 0.0, z0i = 0.0, z1r = 0.0, z1i = 0.0, z2r = 0.0, z2i = 0.0;
                        for (int i = 0; i < K_exp; ++i) {
                            long double g_ld = (long double)Ta + (long double)omids[i];
                            double g = (double)g_ld;
                            double dg = (double)(g_ld - (long double)gb);
                            double den = 0.25 + g * g;
                            double inv_rho_r = 0.5 / den;
                            double inv_rho_i = -g / den;
                            double inv_rho2_r = (0.25 - g * g) / (den * den * u_ref);
                            double inv_rho2_i = (-g) / (den * den * u_ref);
                            double cr = -2.0 * (inv_rho_r + inv_rho2_r);
                            double ci = -2.0 * (inv_rho_i + inv_rho2_i);
                            long double ph_ld = g_ld * (long double)u_ref;
                            double ph = (double)(ph_ld - std::floor(ph_ld / TWO_PI_LD) * TWO_PI_LD);
                            double c_ph = std::cos(ph), s_ph = std::sin(ph);
                            double zr = cr * c_ph - ci * s_ph;
                            double zi = cr * s_ph + ci * c_ph;
                            z0r += zr; z0i += zi;
                            z1r += -dg * zi; z1i += dg * zr;
                            double hdg2 = -0.5 * dg * dg;
                            z2r += hdg2 * zr; z2i += hdg2 * zi;
                        }
                        Z0_re[c] = z0r; Z0_im[c] = z0i;
                        Z1_re[c] = z1r; Z1_im[c] = z1i;
                        Z2_re[c] = z2r; Z2_im[c] = z2i;
                        continue;
                    }
                } else {
                    fclose(fc);
                }
            }
        }

        std::vector<ZeroOffsetInterval> certs = certify_chunk(Ta, Tb, K_exp, log_n_ld, log_n, inv_sqrt_n, 48);
        if ((int)certs.size() != K_exp) {
            certs = certify_chunk(Ta, Tb, K_exp, log_n_ld, log_n, inv_sqrt_n, 512);
        }

        bool ok = ((int)certs.size() == K_exp);
        double mrad = 0.0;
        if (ok && K_exp > 0) {
            if (!(0.0 < certs[0].oa && certs.back().ob < chunk_len)) ok = false;
            for (int i = 0; i < K_exp; ++i) {
                if (!(certs[i].oa < certs[i].ob)) { ok = false; break; }
                if (i > 0 && !(certs[i - 1].ob < certs[i].oa)) { ok = false; break; }
                if (certs[i].rad > mrad) mrad = certs[i].rad;
            }
        }

        if (!ok) {
            status_out[c] = 0;
            max_rad_out[c] = -1.0;
            Z0_re[c] = Z0_im[c] = Z1_re[c] = Z1_im[c] = Z2_re[c] = Z2_im[c] = 0.0;
            continue;
        }

        status_out[c] = 1;
        max_rad_out[c] = mrad;

        if (out_dir && out_dir[0] != '\0') {
            char tmp_path[540];
            snprintf(tmp_path, sizeof(tmp_path), "%s.tmp.%d", fpath, omp_get_thread_num());
            FILE *f = fopen(tmp_path, "wb");
            if (f) {
                int32_t k_hdr = K_exp;
                fwrite(&k_hdr, sizeof(int32_t), 1, f);
                fwrite(&Ta, sizeof(double), 1, f);
                for (int i = 0; i < K_exp; ++i) {
                    fwrite(&certs[i].omid, sizeof(double), 1, f);
                    fwrite(&certs[i].rad, sizeof(float), 1, f);
                }
                fclose(f);
                rename(tmp_path, fpath);
            }
        }

        double z0r = 0.0, z0i = 0.0, z1r = 0.0, z1i = 0.0, z2r = 0.0, z2i = 0.0;
        for (int i = 0; i < K_exp; ++i) {
            long double g_ld = (long double)Ta + (long double)certs[i].omid;
            double g = (double)g_ld;
            double dg = (double)(g_ld - (long double)gb);
            double den = 0.25 + g * g;
            double inv_rho_r = 0.5 / den;
            double inv_rho_i = -g / den;
            double inv_rho2_r = (0.25 - g * g) / (den * den * u_ref);
            double inv_rho2_i = (-g) / (den * den * u_ref);
            double cr = -2.0 * (inv_rho_r + inv_rho2_r);
            double ci = -2.0 * (inv_rho_i + inv_rho2_i);

            long double ph_ld = g_ld * (long double)u_ref;
            double ph = (double)(ph_ld - std::floor(ph_ld / TWO_PI_LD) * TWO_PI_LD);
            double c_ph = std::cos(ph), s_ph = std::sin(ph);
            double zr = cr * c_ph - ci * s_ph;
            double zi = cr * s_ph + ci * c_ph;

            z0r += zr; z0i += zi;
            z1r += -dg * zi; z1i += dg * zr;
            double hdg2 = -0.5 * dg * dg;
            z2r += hdg2 * zr; z2i += hdg2 * zi;
        }
        Z0_re[c] = z0r; Z0_im[c] = z0i;
        Z1_re[c] = z1r; Z1_im[c] = z1i;
        Z2_re[c] = z2r; Z2_im[c] = z2i;
    }
}

// Rigorously evaluate S_K(omega) and its total uncertainty radius over cached binary chunks
// [T_edges[c], T_edges[c+1]] using 80-bit long double + exact Mean Value Theorem Lipschitz bound.
int eval_certified_S_chunks_omp(
    const double *T_edges,
    const int *K_expected,
    int N_chunks,
    const char *out_dir,
    long double omega,
    long double eta1,
    long double delta,
    long double lam1,
    long double lam2,
    long double lam3,
    int M_moll,
    double *chunk_S_mid,
    double *chunk_S_rad,
    int n_threads
) {
    constexpr long double TWO_PI_LD = 6.28318530717958647692528676655900577L;
    long double eta0 = eta1 + delta;
    int any_err = 0;

    #pragma omp parallel for schedule(dynamic, 1) num_threads(n_threads)
    for (int c = 0; c < N_chunks; ++c) {
        if (any_err) continue;
        double Ta = T_edges[c];
        double Tb = T_edges[c + 1];
        int K_exp = K_expected[c];
        char fpath[512];
        snprintf(fpath, sizeof(fpath), "%s/chunk_%08lld_%08lld.bin",
                 out_dir, (long long)std::llround(Ta), (long long)std::llround(Tb));
        FILE *f = fopen(fpath, "rb");
        if (!f) { any_err = 1; continue; }
        int32_t k_hdr = 0;
        double ta_hdr = 0.0;
        if (fread(&k_hdr, sizeof(int32_t), 1, f) != 1 ||
            fread(&ta_hdr, sizeof(double), 1, f) != 1 ||
            k_hdr != K_exp || std::abs(ta_hdr - Ta) > 1e-6) {
            fclose(f); any_err = 1; continue;
        }
        std::vector<double> omids(K_exp);
        std::vector<float> rads(K_exp);
        for (int i = 0; i < K_exp; ++i) {
            if (fread(&omids[i], sizeof(double), 1, f) != 1 ||
                fread(&rads[i], sizeof(float), 1, f) != 1) {
                any_err = 1; break;
            }
        }
        fclose(f);
        if (any_err) continue;

        long double s_sum = 0.0L;
        long double r_sum = 0.0L;
        for (int i = 0; i < K_exp; ++i) {
            long double g = (long double)Ta + (long double)omids[i];
            long double rad = (long double)rads[i];
            long double xi = eta1 * g;
            long double kh_core = lam1 * kappa_2m_hat_ld(xi, 1)
                                + lam2 * kappa_2m_hat_ld(xi, 2)
                                + lam3 * kappa_2m_hat_ld(xi, 3);
            long double kh_moll = kappa_2m_hat_ld(delta * g, M_moll);
            long double w = kh_core * kh_moll;

            long double den = 0.25L + g * g;
            long double inv_rho_r = 0.5L / den;
            long double inv_rho_i = -g / den;
            long double inv_rho2_r = (0.25L - g * g) / (den * den * omega);
            long double inv_rho2_i = (-g) / (den * den * omega);
            long double cr = inv_rho_r + inv_rho2_r;
            long double ci = inv_rho_i + inv_rho2_i;

            long double ph_ld = g * omega;
            long double ph = ph_ld - floorl(ph_ld / TWO_PI_LD) * TWO_PI_LD;
            long double re_part = cr * cosl(ph) - ci * sinl(ph);
            long double term = -2.0L * w * re_part;
            s_sum += term;

            // Rigorous MVT + roundoff bound for term uncertainty:
            long double g_lo = g - rad;
            long double w_ub = fabsl(w) + eta0 * rad;
            long double lip = (w_ub * (2.0L * omega + 3.0L) + eta0 * (1.0L + 1.0L / (g_lo * omega))) / g_lo;
            r_sum += rad * lip + 1e-15L / g_lo;
        }
        chunk_S_mid[c] = (double)s_sum;
        chunk_S_rad[c] = (double)(r_sum + (long double)K_exp * 1e-18L);
    }
    return any_err ? -1 : 0;
}

}
