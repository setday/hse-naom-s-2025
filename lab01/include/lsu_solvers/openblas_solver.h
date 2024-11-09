#pragma once

#include <cstddef>

//#include <cblas.h>

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
  void solve_linear_system_openblas(size_t N, double **A, double *x, const double *b) {
    throw std::runtime_error("Not implemented");
  }
} // namespace ADAAI::LSU_SOLVERS
