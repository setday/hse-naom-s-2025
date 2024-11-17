#pragma once

#include <cstdlib>
#include <random>

namespace ADAAI::LSE_SOLVERS
{

inline std::mt19937& generator()
{
  static thread_local std::mt19937 gen( std::random_device {}() );
  return gen;
}

template<typename T>
T random_float( T min = 0.0, T max = 1.0 )
{
  std::uniform_real_distribution<T> dist( min, max );
  return dist( generator() );
}

const int    K   = 1000; // number of iterations
const double EPS = 1e-4; // epsilon

/**
 * @brief Solves a system of linear equations using iterative scheme.
 *
 * This function takes the coefficient matrix, the x vector and the right-hand
 * side vector and solves for the unknowns using the iterative method. The
 * solution is then stored back in the x vector.
 *
 * @param N The size of the system of equations.
 * @param A The coefficient matrix of the system.
 * @param x The solution vector.
 * @param b The right-hand side vector of the system.
 * @param eps The epsilon.
 * @tparam T The type of the elements in the matrix and vectors (must be a floating point type).
 */
template<typename T>
  requires std::is_floating_point_v<T>
void solve_linear_system_iterative( size_t N, T** A, T* x, const T* b )
{
  if ( N == 0 )
  {
    return;
  }

  if ( N == 1 )
  {
    x[0] = b[0] / A[0][0];
    return;
  }

  for ( size_t i = 0; i < N; ++i )
    x[i] = random_float<T>() + random_float<T>() + random_float<T>();
  //   memset(x, 0, N * sizeof(T)); // maybe ok too

  auto* x_new = new T[N];

  T prev_diff = -1;

  auto** B = new T*[N]; // B = I - A
  for ( size_t i = 0; i < N; ++i )
  {
    B[i] = new T[N];
    for ( size_t j = 0; j < N; ++j )
    {
      B[i][j] = ( i == j ) - A[i][j];
    }
  }

  for ( int it = 0; it < K; ++it )
  { // x = B @ x_new + b
    T max_diff = 0.0;

    for ( size_t i = 0; i < N; ++i )
    {
      x_new[i] = b[i];
      for ( size_t j = 0; j < N; ++j )
      {
        x_new[i] += B[i][j] * x[j];
      }

      // l_inf norm, maybe use l_2 (also, "rule of thumb" mentioned?)
      max_diff = std::max( max_diff, std::abs( x[i] - x_new[i] ) );
    }

    if ( max_diff < EPS )
    {
      prev_diff = max_diff;
      break;
    }

    if ( it > K / 2 && max_diff > prev_diff )
    {
      // we require monotonic decrease after K / 2 iterations
      throw std::runtime_error( "Iterative solver diverged (the difference increased)" );
    }

    std::swap( x_new, x );
    prev_diff = max_diff;
  }
  if ( prev_diff >= EPS )
  {
    throw std::runtime_error( "Iterative solver diverged (the difference at the end is too big)" );
  }
}
} // namespace ADAAI::LSE_SOLVERS
