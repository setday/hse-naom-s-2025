#include <iostream>
#include <vector>
#include <tuple>
#include <fstream>

#include "../include/lsu_solvers/solve_lsu.h"
#include "include/memory_management.h"
#include "include/random_management.h"

template <typename F>
void find_inverse(F **A, F *b, F *x, size_t N) {
  ADAAI::LSU_SOLVERS::solve_linear_system(N, A, x, b, ADAAI::LSU_SOLVERS::LSSolveMethod::OPENBLAS);
}

template <typename F>
F compute_eps(F **A, F *b, F *x, size_t N) {
  F eps = 0.0;

  for (size_t i = 0; i < N; ++i) {
    F AXb = -b[i];

    for (size_t j = 0; j < N; ++j) {
      AXb += A[i][j] * x[j];
    }

    eps += AXb * AXb;
  }

  return eps;
}

template <typename F>
std::pair<F, F> test_lsu_solvers(size_t N) {
  F **A = ADAAI::LAB01::TEST::MEM::allocate_matrix<F>(N);
  F *b = ADAAI::LAB01::TEST::MEM::allocate_vector<F>(N);
  F *x = ADAAI::LAB01::TEST::MEM::allocate_vector<F>(N);

  size_t num_tests = 10;
  std::vector<F> eps(num_tests);

  for (size_t i = 0; i < num_tests; ++i) {
    ADAAI::LAB01::TEST::RAND::generate_normally_distributed_matrix(A, N);
    ADAAI::LAB01::TEST::RAND::generate_normally_distributed_vector(b, N);

    find_inverse(A, b, x, N);

    eps[i] = compute_eps(A, b, x, N);
  }

  ADAAI::LAB01::TEST::MEM::deallocate_matrix(A, N);
  ADAAI::LAB01::TEST::MEM::deallocate_vector(b);
  ADAAI::LAB01::TEST::MEM::deallocate_vector(x);

  F sum_eps = 0.0;
  F max_eps = 0.0;

  for (size_t i = 0; i < num_tests; ++i) {
    sum_eps += eps[i];
    max_eps = std::max(max_eps, eps[i]);
  }

  return {sum_eps / num_tests, max_eps};
}

void dump_results(const std::vector<std::tuple<size_t, double, double>>& results) {
  for (const auto& [N, avg_eps, max_eps] : results) {
    std::cout << N << " " << avg_eps << " " << max_eps << std::endl;
  }

  // save to file
  std::ofstream file("GEP_eps_results.json");
  file << "[\n";
  for (const auto& [N, avg_eps, max_eps] : results) {
    file << "  {\"N\": " << N << ", \"avg_eps\": " << avg_eps << ", \"max_eps\": " << max_eps << "},\n";
  }
  file << "]\n";
}

void test_lsu_solvers() {
  size_t min_N = 10;
  size_t max_N = 300;

  std::vector<std::tuple<size_t, double, double>> results(max_N - min_N + 1);

  for (size_t N = min_N; N <= max_N; ++N) {
    auto [avg_eps, max_eps] = test_lsu_solvers<double>(N);
    results[N - min_N] = {N, avg_eps, max_eps};

    std::cout << static_cast<double>(N - min_N) / static_cast<double>(max_N - min_N) * 100 << "%\n";
  }

  dump_results(results);
}
