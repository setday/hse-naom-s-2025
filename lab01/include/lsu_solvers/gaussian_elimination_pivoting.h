#pragma once

#include <cstddef>
#include <cstring>
#include <cassert>

namespace ADAAI::LSU_SOLVERS {
  /**
   * @brief Solves a system of linear equations using GaussElimPivoting without
   * GSL.
   *
   * This function takes a reference to a SLAE object, extracts the coefficient
   * matrix and the right-hand side vector, and solves for the unknowns using the
   * GaussElimPivoting method. The solution is then stored back in the SLAE
   * object.
   *
   * @param N The size of the system of equations.
   * @param A The coefficient matrix of the system.
   * @param x The solution vector (also the right-hand side vector).
   * @param b The right-hand side vector of the system.
   */
  void solve_linear_system_GEP(size_t N, double **A, double *x, const double *b) {
    if (N == 0) {
      return;
    }

    memcpy(x, b, N * sizeof(double));

    if (N == 1) {
      x[0] /= A[0][0];
      return;
    }

    auto **A_copy = new double *[N];
    for (size_t i = 0; i < N; i++) {
      A_copy[i] = new double[N];
      memcpy(A_copy[i], A[i], N * sizeof(double));
    }

    auto *p = new size_t[N];
    for (size_t i = 0; i < N; i++) {
      p[i] = i;
    }

    // Perform Gaussian elimination with pivoting
    for (size_t i = 0; i < N; i++) {
      // Find pivot
      size_t pivot = i;
      for (size_t j = i + 1; j < N; j++) {
        if (std::abs(A_copy[j][i]) > std::abs(A_copy[pivot][i])) {
          pivot = j;
        }
      }

      if (A_copy[pivot][i] == 0) {
        throw std::runtime_error("GEP solver: Singular matrix encountered");
      }

      // Swap rows
      std::swap(A_copy[i], A_copy[pivot]);
      std::swap(x[i], x[pivot]);
      std::swap(p[i], p[pivot]);

      // Eliminate
      for (size_t j = i + 1; j < N; j++) {
        double factor = A_copy[j][i] / A_copy[i][i];
        A_copy[j][i] = 0;
        for (size_t k = i + 1; k < N; k++) {
          A_copy[j][k] -= factor * A_copy[i][k];
        }
        x[j] -= factor * x[i];
      }
    }

    // Back substitution
    for (int i = N - 1; i >= 0; i--) {
      for (size_t j = i + 1; j < N; j++) {
        x[i] -= A_copy[i][j] * x[j];
      }
      x[i] /= A_copy[i][i];
    }

    // Reorder the solution
    for (size_t i = 0; i < N; i++) {
      while (p[i] != i) {
        std::swap(x[i], x[p[i]]);
        std::swap(p[i], p[p[i]]);
      }
    }

    for (size_t i = 0; i < N; i++) {
      delete[] A_copy[i];
    }
    delete[] A_copy;
    delete[] p;
  }
} // namespace ADAAI::LSU_SOLVERS
