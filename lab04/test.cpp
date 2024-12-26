#include <iostream>

#include "include/buchberger.h"

using namespace ADAAI::LAB04::Poly;

int main()
{
  Monomial<double> ms[] = {
      {3, 2, new size_t[2]{1, 0}},
      {1, 2, new size_t[2]{0, 1}},
      {0.5, 2, new size_t[2]{1, 1}},
      {1, 2, new size_t[2]{2, 0}},
  }; // 3x, y, 0.5xy, x^2

  Polynomial<double> ps[] = {
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
