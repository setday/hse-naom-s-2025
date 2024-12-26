#include <iostream>

#include "include/buchberger.h"
#include "include/GMPwrapper.hpp"

using namespace ADAAI::LAB04::Poly;

int main()
{
  Monomial<Rational> ms[] = {
      { { 1, 1 }, 2, new size_t[2]{0, 1}},
      { { 3, 1 }, 2, new size_t[2]{1, 0}},
      { { 1, 2 }, 2, new size_t[2]{1, 1}},
      { { 1, 1 }, 2, new size_t[2]{2, 0}},
  }; // 3x, y, 0.5xy, x^2

  Polynomial<Rational> ps[] = {
      {2, ms},
      {2, ms + 2},
  }; // 3x + y, 0.5xy + x^2

//  GroebnerBasis<double> gb = buchberger({ms, ms + 2});

  std::cout << "Monomials:\n";
  for (const auto& m : ms)
  {
    std::cout << m << " | ";
  }

  std::cout << "\n\nPolynomials:\n";

  for (const auto& p : ps)
  {
    std::cout << p << '\n';
  }

  return 0;
}
