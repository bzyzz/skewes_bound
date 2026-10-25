#include <cmath>
#include <cstdint>
#include <vector>
#include <omp.h>

extern "C" {

// Fast evaluation on a uniform grid u_k = u0 + k * du (k = 0 .. M-1)
// using phase stepping with periodic exact resync every 256 steps.
void eval_skewes_uniform(
    const double *gamma,
    const double *cc,
    const double *cs,
    int N,
    double u0,
    double du,
    int M,
    double *out
) {
    for (int k = 0; k < M; ++k) out[k] = 0.0;

    #pragma omp parallel num_threads(32)
    {
        std::vector<double> local_out(M, 0.0);
        #pragma omp for schedule(static, 4096)
        for (int j = 0; j < N; ++j) {
            double g = gamma[j];
            double a = cc[j];
            double b = cs[j];
            double dc = std::cos(g * du);
            double ds = std::sin(g * du);
            for (int k_base = 0; k_base < M; k_base += 256) {
                int k_end = std::min(M, k_base + 256);
                double ang0 = g * (u0 + k_base * du);
                double c = std::cos(ang0);
                double s = std::sin(ang0);
                for (int k = k_base; k < k_end; ++k) {
                    local_out[k] += a * c + b * s;
                    double nc = c * dc - s * ds;
                    s = s * dc + c * ds;
                    c = nc;
                }
            }
        }
        #pragma omp critical
        {
            for (int k = 0; k < M; ++k) out[k] += local_out[k];
        }
    }
}

// Evaluate at arbitrary points u_grid[0..M-1], parallelizing over zeros when M is small
void eval_skewes_pts(
    const double *gamma,
    const double *cc,
    const double *cs,
    int N,
    const double *u_grid,
    int M,
    double *out
) {
    for (int k = 0; k < M; ++k) out[k] = 0.0;
    #pragma omp parallel num_threads(32)
    {
        std::vector<double> local_out(M, 0.0);
        #pragma omp for schedule(static, 4096)
        for (int j = 0; j < N; ++j) {
            double g = gamma[j];
            double a = cc[j];
            double b = cs[j];
            for (int k = 0; k < M; ++k) {
                double ang = g * u_grid[k];
                local_out[k] += a * std::cos(ang) + b * std::sin(ang);
            }
        }
        #pragma omp critical
        {
            for (int k = 0; k < M; ++k) out[k] += local_out[k];
        }
    }
}

}
