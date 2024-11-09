#pragma once

#include <cstddef>
#include <stdexcept>

#include "gsl_solver.h"
#include "gaussian_elimination_pivoting.h"
#include "openblas_solver.h"

namespace ADAAI::LSU_SOLVERS {
  enum class LSSolveMethod {
    GSL,
    GEP,
    OPENBLAS,
  };

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
   * @param method The method to use for solving the system.
   */
  void solve_linear_system(size_t N, double **A, double *x, const double *b, LSSolveMethod method) {
    switch (method) {
      case LSSolveMethod::GSL:
        solve_linear_system_gsl(N, A, x, b);
        break;
      case LSSolveMethod::GEP:
        solve_linear_system_GEP(N, A, x, b);
        break;
      case LSSolveMethod::OPENBLAS:
        solve_linear_system_openblas(N, A, x, b);
        break;
      default:
        throw std::runtime_error("Unknown method");
        break;
    }
  }
} // namespace ADAAI::LSU_SOLVERS
