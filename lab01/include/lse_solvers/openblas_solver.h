#pragma once

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <openblas_config.h>

extern "C"
{
  void dgetrf_( const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info );
  void dgetrs_( const char* trans, const int* n, const int* nrhs, const double* a, const int* lda, const int* ipiv, double* b, const int* ldb, int* info );
}

namespace ADAAI::LSE_SOLVERS
{
/**
 * @brief Solves a system of linear equations using LU decomposition with
 * OpenBLAS.
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
void solve_linear_system_openblas( std::size_t N, double** A, double* x, const double* b )
{
  int* ipiv = reinterpret_cast<int*>( std::malloc( N * sizeof( *ipiv ) ) );
  if ( ipiv == nullptr )
  {
    throw std::bad_alloc();
  }
  double* C = reinterpret_cast<double*>( std::malloc( N * N * sizeof( *C ) ) );
  if ( C == nullptr )
  {
    free( ipiv );
    throw std::bad_alloc();
  }

  std::memcpy( x, b, N * sizeof( double ) );
  for ( std::size_t i = 0; i < N; ++i )
  {
    std::memcpy( C + i * N, A[i], N * sizeof( double ) );
  }

  int info;

  int n_int = static_cast<int>( N );
  dgetrf_( &n_int, &n_int, C, &n_int, ipiv, &info );
  if ( info != 0 )
  {
    free( C );
    free( ipiv );
    throw std::runtime_error( "LU decomposition failed" );
  }

  char trans = 'N';
  int  nrhs  = 1;
  dgetrs_( &trans, &n_int, &nrhs, C, &n_int, ipiv, x, &n_int, &info );
  if ( info != 0 )
  {
    free( C );
    free( ipiv );
    throw std::runtime_error( "Solving the linear system failed" );
  }

  free( C );
  free( ipiv );
}
} // namespace ADAAI::LSE_SOLVERS
