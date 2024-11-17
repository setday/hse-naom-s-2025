#pragma once

#include <cstddef>
#include <stdexcept>

#include "../matrix_work/matrix.h"

#include "gaussian_elimination_pivoting.h"
#include "gsl_solver.h"
#include "iterative_solver.h"
#include "openblas_solver.h"

namespace ADAAI::LSE_SOLVERS
{
enum class LSSolveMethod
{
  GSL,
  GEP,
  OPENBLAS,
  ITERATIVE,
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
void solve_linear_system( size_t                             N,
                          const ADAAI::MATH::Matrix<double>& A,
                          double*                            x,
                          const double*                      b,
                          LSSolveMethod                      method )
{
  const ADAAI::MATH::Matrix<double>& A_copy = A;
  switch ( method )
  {
    case LSSolveMethod::GSL:
      solve_linear_system_gsl( N, A_copy.data, x, b );
      break;
    case LSSolveMethod::GEP:
      solve_linear_system_GEP( N, A_copy.data, x, b );
      break;
    case LSSolveMethod::OPENBLAS:
      solve_linear_system_openblas( N, A_copy.data, x, b );
      break;
    case LSSolveMethod::ITERATIVE:
      solve_linear_system_iterative<double>( N, A_copy.data, x, b );
      break;
    default:
      throw std::runtime_error( "Unknown method" );
  }
}
} // namespace ADAAI::LSE_SOLVERS
