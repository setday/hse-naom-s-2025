#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>

#include "../include/matrix_work/matrix.h"

#include "../include/lse_solvers/solve_lse.h"
#include "include/random_management.h"

template<typename F, ADAAI::LSE_SOLVERS::LSSolveMethod M>
void find_inverse( ADAAI::MATH::Matrix<F> A, F* b, F* x, size_t N )
{
  ADAAI::LSE_SOLVERS::solve_linear_system( N, A, x, b, M );
}

template<typename F>
F compute_eps( F** A, F* b, F* x, size_t N )
{
  F eps = 0.0;

  for ( size_t i = 0; i < N; ++i )
  {
    F AXb = -b[i];

    for ( size_t j = 0; j < N; ++j )
    {
      AXb += A[i][j] * x[j];
    }

    eps += AXb * AXb;
  }

  return eps;
}

template<typename F, ADAAI::LSE_SOLVERS::LSSolveMethod M>
std::pair<F, F> test_lse_solvers( size_t N )
{
  ADAAI::MATH::Matrix<F> A( N, N );
  ADAAI::MATH::Vector<F> b( N );
  ADAAI::MATH::Vector<F> x( N );

  size_t         num_tests = 10;
  std::vector<F> eps( num_tests );

  for ( size_t i = 0; i < num_tests; ++i )
  {
    ADAAI::LAB01::TEST::RAND::generate_normally_distributed_matrix( A );
    ADAAI::LAB01::TEST::RAND::generate_normally_distributed_vector( b );

    find_inverse<F, M>( A, b.data, x.data, N );

    eps[i] = compute_eps( A.data, b.data, x.data, N );
  }

  F sum_eps = 0.0;
  F max_eps = 0.0;

  for ( size_t i = 0; i < num_tests; ++i )
  {
    sum_eps += eps[i];
    max_eps = std::max( max_eps, eps[i] );
  }

  return { sum_eps / num_tests, max_eps };
}

void dump_results(
    const std::vector<std::tuple<size_t, double, double>>& results )
{
  for ( const auto& [N, avg_eps, max_eps] : results )
  {
    std::cout << N << " " << avg_eps << " " << max_eps << std::endl;
  }

  // save to file
  std::ofstream file( "GEP_eps_results.json" );
  file << "[\n";
  for ( const auto& [N, avg_eps, max_eps] : results )
  {
    file << "  {\"N\": " << N << ", \"avg_eps\": " << avg_eps
         << ", \"max_eps\": " << max_eps << "},\n";
  }
  file << "]\n";
}

template<typename F = double, ADAAI::LSE_SOLVERS::LSSolveMethod M = ADAAI::LSE_SOLVERS::LSSolveMethod::GEP>
void test_lse_solvers()
{
  size_t min_N = 10;
  size_t max_N = 300;

  std::vector<std::tuple<size_t, double, double>> results( max_N - min_N + 1 );

  for ( size_t N = min_N; N <= max_N; ++N )
  {
    auto [avg_eps, max_eps] = test_lse_solvers<F, M>( N );
    results[N - min_N]      = { N, avg_eps, max_eps };

    std::cout << static_cast<double>( N - min_N ) /
                     static_cast<double>( max_N - min_N ) * 100
              << "%\n";
  }

  dump_results( results );
}
