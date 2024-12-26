#pragma once

#include <vector>

#include "polynomial.h"

namespace ADAAI::LAB04::Poly
{

template <typename T>
struct GroebnerBasis
{
  std::vector<Polynomial<T>> basis;

  void add(const Polynomial<T>& m)
  {
    basis.push_back(m);
  }
};

} // namespace ADAAI::LAB04::Poly
