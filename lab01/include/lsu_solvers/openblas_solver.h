#pragma once

#include <cstdlib>
#include <cstddef>
#include <gsl/gsl_blas.h>
#include <lapacke.h>

namespace ADAAI::LSU_SOLVERS {
  /**
   * @brief Solves a system of linear equations using LU decomposition with OpenBLAS.
   *
   * This function takes a reference to a SLAE object, extracts the coefficient
   * matrix and the right-hand side vector, and solves for the unknowns using the
   * OpenBLAS library. The solution is then stored back in the SLAE object.
   *
   * @param N The size of the system of equations.
   * @param A The coefficient matrix of the system.
   * @param x The solution vector.
   * @param b The right-hand side vector of the system.
   */
  void solve_linear_system_openblas(std::size_t N, double **A, double *x, const double *b) {
    int *ipiv = (int *)std::malloc(N * sizeof(*ipiv));
    std::memcpy(x, b, N * sizeof(double));
    double *C = (double *)std::malloc(N * N * sizeof(*C));
    for(std::size_t i = 0; i < N; ++i) {
        std::memcpy(C + i * N, A[i], N * sizeof(double));
    }
    LAPACKE_dgetrf(LAPACK_COL_MAJOR, N, N, C, N, ipiv);
    LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', N, 1, C, N, ipiv, x, N);
    free(C);
    free(ipiv);
  }
} // namespace ADAAI::LSU_SOLVERS
