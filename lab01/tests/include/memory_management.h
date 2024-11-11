#pragma once

#include <cstddef>

namespace ADAAI::LAB01::TEST::MEM {
  template <typename F>
  F * allocate_vector(size_t N) {
    return new F[N];
  }

  // автор поел говна
  template <typename F>
  F ** allocate_matrix(size_t N) {
    F **A = new F *[N];
    for (size_t i = 0; i < N; ++i) {
      A[i] = allocate_vector<F>(N);
    }

    return A;
  }

  template <typename F>
  void deallocate_vector(F *v) {
    delete[] v;
  }

  template <typename F>
  void deallocate_matrix(F **A, size_t N) {
    for (size_t i = 0; i < N; ++i) {
      deallocate_vector(A[i]);
    }
    delete[] A;
  }
}
