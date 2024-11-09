#pragma once

#include <cstring>
#include <iostream>

//#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_vector.h>

#include "SLAE.h"

namespace ADAAI {
enum class LSSolveMethod {
  GSL,
  GEP,
};

/**
 * @brief Solves a system of linear equations using LU decomposition with GSL.
 *
 * This function takes a reference to a SLAE object, extracts the coefficient
 * matrix and the right-hand side vector, and solves for the unknowns using the
 * GNU Scientific Library (GSL). The solution is then stored back in the SLAE
 * object.
 *
 * @param N The size of the system of equations.
 * @param A The coefficient matrix of the system.
 * @param x The solution vector.
 * @param b The right-hand side vector of the system.
 */
void solve_linear_system_gsl(size_t N, double **A, double *x, const double *b) {
  // Create GSL matrix and vectors
  gsl_matrix *gsl_A = gsl_matrix_alloc(N, N);
  gsl_vector *gsl_b = gsl_vector_alloc(N);
  gsl_vector *gsl_x = gsl_vector_alloc(N);

  // Fill GSL matrix A
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) {
      gsl_matrix_set(gsl_A, i, j, A[i][j]);
    }
  }

  // Fill GSL vector b
  for (size_t i = 0; i < N; i++) {
    gsl_vector_set(gsl_b, i, b[i]);
  }

  // Decompose the matrix
  gsl_permutation *p = gsl_permutation_alloc((int)N);
  int signum;
  gsl_linalg_LU_decomp(gsl_A, p, &signum);

  // Solve the system
  gsl_linalg_LU_solve(gsl_A, p, gsl_b, gsl_x);

  // Fill in the solution vector x
  for (size_t i = 0; i < N; i++) {
    x[i] = gsl_vector_get(gsl_x, i);
  }

  // Free allocated memory
  gsl_permutation_free(p);
  gsl_vector_free(gsl_b);
  gsl_vector_free(gsl_x);
  gsl_matrix_free(gsl_A);
}

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
    default:
      std::cerr << "Unknown method" << std::endl;
      break;
  }
}
} // namespace ADAAI
