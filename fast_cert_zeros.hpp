#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <algorithm>
#include <omp.h>

static inline long double rs_theta(long double t) {
    const long double pi = 3.14159265358979323846264338327950288L;
    long double half_t = 0.5L * t;
    long double log_t_2pi = std::log(t / (2.0L * pi));
    long double inv_t = 1.0L / t;
    long double inv_t2 = inv_t * inv_t;
    return half_t * log_t_2pi - half_t - 0.125L * pi
           + inv_t / 48.0L
           + 7.0L * inv_t * inv_t2 / 5760.0L;
}

static inline double rs_remainder_01(long double t_ld, int &M_out) {
    static constexpr double P_PSI0[13] = {
        1.18052391178158452e-09, 1.44720063736079051e-08, -1.07379560253753367e-07,
        -1.42845983385324725e-06, 4.62498401730884284e-07, 7.94344157397609076e-05,
        2.97053155469926026e-04, -1.62372528016937068e-03, -1.35676219644003601e-02,
        -1.36050260500901795e-02, 1.32376575480621400e-01, 4.37240468077510658e-01,
        3.82683432365090559e-01
    };
    static constexpr double P_C1_Z[13] = {
        2.58702233523163268e-09, 5.44246895386571258e-09, -1.85497573774654085e-07,
        -9.61489139738042587e-07, 5.97817295473849217e-06, 5.92285661512381844e-05,
        -1.32083184567283487e-05, -1.46478077033620861e-03, -3.31075977373197685e-03,
        9.87106630018547669e-03, 3.84912504821840543e-02, 1.37847734263526738e-02,
        -2.68251026283753587e-02
    };
    long double a_ld = std::sqrt(t_ld / (2.0L * 3.14159265358979323846264338327950288L));
    int M = (int)std::floor(a_ld);
    M_out = M;
    double a = (double)a_ld;
    double p = (double)(a_ld - (long double)M);
    double z = 2.0 * p - 1.0;
    double u = z * z;
    double psi0 = P_PSI0[0], c1_z = P_C1_Z[0];
    for (int k = 1; k < 13; ++k) {
        psi0 = psi0 * u + P_PSI0[k];
        c1_z = c1_z * u + P_C1_Z[k];
    }
    double c1 = (z * c1_z) / a;
    double sign = (M & 1) ? 1.0 : -1.0;
    return sign * (psi0 + c1) / std::sqrt(a);
}

static inline double rigorous_err_Z(double t) {
    return 0.053 * std::pow(t, -1.25) + 2.2e-11;
}

// Offset interval relative to chunk start T_a: gamma in [T_a + oa, T_a + ob]
struct ZeroOffsetInterval {
    double oa, ob;
    double omid;
    float rad;
};

std::vector<ZeroOffsetInterval> certify_chunk(
    double T_a, double T_b, int expected_count,
    const std::vector<long double> &log_n_ld,
    const std::vector<double> &log_n,
    const std::vector<double> &inv_sqrt_n,
    int n_scan_base
) {
    constexpr int K = 22;
    constexpr double H = 0.25;
    constexpr long double TWO_PI_LD = 6.28318530717958647692528676655900577L;
    std::vector<ZeroOffsetInterval> certs;
    certs.reserve(expected_count + 16);

    double last_off = 0.0; // offset from T_a
    double last_z = 0.0;
    int last_s = 0;

    int cur_M = 0;
    std::vector<double> cn, sn, dcn, dsn;
    int step_count = 0;

    int n_blocks = (int)std::llround((T_b - T_a) / (2.0 * H));
    for (int blk = 0; blk < n_blocks; ++blk, ++step_count) {
        double off_c = (blk + 0.5) * (2.0 * H); // exact multiple of 0.25 in (0, T_b - T_a)
        double t_c = T_a + off_c;
        long double tc_ld = (long double)T_a + (long double)off_c;
        int M = (int)std::floor(std::sqrt((tc_ld - (long double)H) / (2.0L * 3.14159265358979323846264338327950288L)));
        if (step_count % 64 == 0 || M != cur_M) {
            cur_M = M;
            cn.resize(M + 1);
            sn.resize(M + 1);
            dcn.resize(M + 1);
            dsn.resize(M + 1);
            for (int n = 1; n <= M; ++n) {
                long double ang_ld = tc_ld * log_n_ld[n];
                double ang = (double)(ang_ld - std::floor(ang_ld / TWO_PI_LD) * TWO_PI_LD);
                cn[n] = inv_sqrt_n[n] * std::cos(ang);
                sn[n] = -inv_sqrt_n[n] * std::sin(ang);
                double step_ang = 2.0 * H * log_n[n];
                dcn[n] = std::cos(step_ang);
                dsn[n] = -std::sin(step_ang);
            }
        } else {
            for (int n = 1; n <= M; ++n) {
                double c = cn[n], s = sn[n];
                double dc = dcn[n], ds = dsn[n];
                cn[n] = c * dc - s * ds;
                sn[n] = s * dc + c * ds;
            }
        }

        double Ar[K] = {0.0}, Ai[K] = {0.0};
        for (int n = 1; n <= M; ++n) {
            double ln = log_n[n];
            double c = cn[n], s = sn[n];
            double p_ln = 1.0;
            for (int k = 0; k < K; ++k) {
                if (k > 0) p_ln *= ln / (double)k;
                int mod4 = k & 3;
                if (mod4 == 0)      { Ar[k] += c * p_ln; Ai[k] += s * p_ln; }
                else if (mod4 == 1) { Ar[k] += s * p_ln; Ai[k] -= c * p_ln; }
                else if (mod4 == 2) { Ar[k] -= c * p_ln; Ai[k] -= s * p_ln; }
                else                { Ar[k] -= s * p_ln; Ai[k] += c * p_ln; }
            }
        }

        long double th0_ld = rs_theta(tc_ld);
        double th0 = (double)(th0_ld - std::floor(th0_ld / TWO_PI_LD) * TWO_PI_LD);
        double th1 = (double)(0.5L * std::log(tc_ld / (2.0L * 3.14159265358979323846264338327950288L)) - 1.0L / (48.0L * tc_ld * tc_ld));
        double th2 = 0.25 / t_c;

        auto eval_Z_offset = [&](double d) -> double {
            long double t_ld = tc_ld + (long double)d;
            int M_pt = 0;
            double rem = rs_remainder_01(t_ld, M_pt);
            double fr = Ar[K - 1], fi = Ai[K - 1];
            for (int k = K - 2; k >= 0; --k) {
                fr = fr * d + Ar[k];
                fi = fi * d + Ai[k];
            }
            for (int n = M + 1; n <= M_pt; ++n) {
                long double ang_ld = t_ld * log_n_ld[n];
                double ang = (double)(ang_ld - std::floor(ang_ld / TWO_PI_LD) * TWO_PI_LD);
                fr += inv_sqrt_n[n] * std::cos(ang);
                fi -= inv_sqrt_n[n] * std::sin(ang);
            }
            double th = th0 + d * (th1 + d * th2);
            double c_th = std::cos(th), s_th = std::sin(th);
            return 2.0 * (c_th * fr - s_th * fi) + rem;
        };

        auto cert_sign_offset = [&](double d, double &z_out) -> int {
            double z = eval_Z_offset(d);
            z_out = z;
            double err = rigorous_err_Z(t_c + d);
            if (z > err) return +1;
            if (z < -err) return -1;
            return 0;
        };

        auto process_point = [&](double cur_d) {
            double cur_off = off_c + cur_d;
            double cur_z = 0.0;
            int cur_s = cert_sign_offset(cur_d, cur_z);
            if (cur_s == 0) return;
            if (last_s != 0 && cur_s != last_s) {
                double da = last_off - off_c, db = cur_d;
                double za = last_z, zb = cur_z;
                int sa = last_s, sb = cur_s;
                double ra = da, rb = db, rza = za, rzb = zb;
                for (int it = 0; it < 24; ++it) {
                    if (rb - ra < 1e-13) break;
                    double mid = 0.5 * (ra + rb);
                    double sec = ra - rza * (rb - ra) / (rzb - rza);
                    if (sec <= ra + 0.05 * (rb - ra) || sec >= rb - 0.05 * (rb - ra) || (it % 3 == 2)) sec = mid;
                    double zs = eval_Z_offset(sec);
                    if ((rza < 0.0 && zs > 0.0) || (rza > 0.0 && zs < 0.0)) { rb = sec; rzb = zs; }
                    else { ra = sec; rza = zs; }
                }
                double root_d = 0.5 * (ra + rb);
                for (double eps : {3e-10, 6e-10, 1.5e-9, 4e-9, 1e-8, 5e-8, 2e-7}) {
                    double cand_a = std::max(da, root_d - eps);
                    double cand_b = std::min(db, root_d + eps);
                    if (!((off_c + cand_a) < (off_c + cand_b))) continue;
                    double cza, czb;
                    if (cert_sign_offset(cand_a, cza) == sa && cert_sign_offset(cand_b, czb) == sb) {
                        da = cand_a; db = cand_b;
                        break;
                    }
                }
                double oa = off_c + da, ob = off_c + db;
                double omid = 0.5 * (oa + ob);
                float rad = (float)(0.5 * (ob - oa) + 1e-12);
                certs.push_back({oa, ob, omid, rad});
            }
            last_off = cur_off;
            last_z = cur_z;
            last_s = cur_s;
        };

        if (blk == 0) {
            last_s = cert_sign_offset(-H, last_z);
            last_off = 0.0;
        }

        double dt = (2.0 * H) / n_scan_base;
        double prev_d = -H;
        double prev_z_abs = std::abs(last_z);
        for (int i = 1; i <= n_scan_base; ++i) {
            double cur_d = (i == n_scan_base) ? H : (-H + i * dt);
            double z_probe = eval_Z_offset(cur_d);
            if (prev_z_abs < 0.25 || std::abs(z_probe) < 0.25) {
                constexpr int N_SUB = 64;
                double sub_dt = (cur_d - prev_d) / N_SUB;
                for (int s = 1; s <= N_SUB; ++s) {
                    double sd = (s == N_SUB) ? cur_d : (prev_d + s * sub_dt);
                    process_point(sd);
                }
            } else {
                process_point(cur_d);
            }
            prev_d = cur_d;
            prev_z_abs = std::abs(z_probe);
        }
    }
    return certs;
}
