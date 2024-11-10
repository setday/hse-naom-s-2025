#pragma once

#include <cstdlib>

namespace ADAAI::LSU_SOLVERS {

  const int K = 1000; // number of iterations

  /**
   * @brief Solves a system of linear equations using interactive scheme.
   *
   * This function takes the coefficient matrix, the x vector and the right-hand side vector
   * and solves for the unknowns using the interactive method. 
   * The solution is then stored back in the x vector.
   *
   * @param N The size of the system of equations.
   * @param A The coefficient matrix of the system.
   * @param x The solution vector.
   * @param b The right-hand side vector of the system.
   * @param eps The epsilon.
   */
  void solve_linear_system_interactive(size_t N, double **A, double *x, const double *b, const double eps = 1e-4) {
    if (N == 0) {
      return;
    }

    if (N == 1) {
      x[0] = b[0] / A[0][0];
      return;
    }

    for (size_t i = 0; i < N; ++i)
      x[i] = rand() + rand() + rand();
    // memcpy(x, b, N * sizeof(double)); // maybe ok too

    auto *x_old = new double[N];
    memcpy(x_old, x, N * sizeof(double));

    double prev_diff = -1;

    auto **B = new double *[N];// B = I - A
    for (size_t i = 0; i < N; ++i) {
      B[i] = new double[N];
      for (size_t j = 0; j < N; ++j) {
        B[i][j] = (i == j) - A[i][j];
      }
    }

    for (int it = 0; it < K; ++it) { // x = B @ x_old + b
        for (size_t i = 0; i < N; ++i) {
            x[i] = 0;
            for (size_t j = 0; j < N; ++j) {
                x[i] += B[i][j] * x_old[j];
            }
            x[i] += b[i];
        }


        double max_diff = 0.0; // l_inf norm, maybe use l_2 (also, "rule of thumb" mentioned?)
        for (size_t i = 0; i < N; ++i) {
            max_diff = std::max(max_diff, std::abs(x[i] - x_old[i]));
        }
        if (max_diff < eps) {
            break;
        }

        if (i > K / 2 && max_diff > prev_diff) { 
            // we require monotonic decrease after K / 2 iterations
            throw std::runtime_error("Iterative solver diverged");
        }

        std::memcpy(x_old, x, N * sizeof(double));
        prev_diff = max_diff;
    }
    if (prev_diff >= eps) {
        throw std::runtime_error("Iterative solver diverged");
    }
  }
} // namespace ADAAI::LSU_SOLVERS
