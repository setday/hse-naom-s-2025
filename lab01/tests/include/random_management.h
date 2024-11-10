#pragma once

#include <cstddef>
#include <random>
#include <ctime>

namespace ADAAI::LAB01::TEST::RAND {
  template <typename F>
  F generate_normally_distributed_double() {
    static std::default_random_engine generator(42);
    static std::normal_distribution<F> distribution(0.0, 1.0);

    return distribution(generator);
  }

  template <typename F>
  void generate_normally_distributed_vector(F *v, size_t N) {
    for (size_t i = 0; i < N; ++i) {
      v[i] = generate_normally_distributed_double<F>();
    }
  }

  template <typename F>
  void generate_normally_distributed_matrix(F **A, size_t N) {
    for (size_t i = 0; i < N; ++i) {
      generate_normally_distributed_vector(A[i], N);
    }
  }
}
