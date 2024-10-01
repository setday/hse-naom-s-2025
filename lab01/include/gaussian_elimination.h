#pragma once

#include "SLAE.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <stdio.h>

namespace ADAII {

/**
 * @brief Solves a system of linear equations using LU decomposition.
 * 
 * This function takes a reference to a SLAE object, extracts the coefficient matrix
 * and the right-hand side vector, and solves for the unknowns using the GNU Scientific
 * Library (GSL). The solution is then stored back in the SLAE object.
 * 
 * @param slae A reference to the SLAE object containing the system of equations.
 */
void solve_linear_system(SLAE &slae) {
  int N = slae.N;
  // Create GSL matrix and vectors
  gsl_matrix *gsl_A = gsl_matrix_alloc(N, N);
  gsl_vector *gsl_b = gsl_vector_alloc(N);
  gsl_vector *gsl_x = gsl_vector_alloc(N);

  // Fill GSL matrix A
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) {
      gsl_matrix_set(gsl_A, i, j, slae.A[i][j]);
    }
  }

  // Fill GSL vector b
  for (size_t i = 0; i < N; i++) {
    gsl_vector_set(gsl_b, i, slae.b[i]);
  }

  // Decompose the matrix
  gsl_permutation *p = gsl_permutation_alloc(N);
  int signum;
  gsl_linalg_LU_decomp(gsl_A, p, &signum);

  // Solve the system
  gsl_linalg_LU_solve(gsl_A, p, gsl_b, gsl_x);

  // Fill in the solution vector x
  for (size_t i = 0; i < N; i++) {
    slae.x[i] = gsl_vector_get(gsl_x, i);
  }

  // Free allocated memory
  gsl_permutation_free(p);
  gsl_vector_free(gsl_b);
  gsl_vector_free(gsl_x);
  gsl_matrix_free(gsl_A);
}
} // namespace ADAII
