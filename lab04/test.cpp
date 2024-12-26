#include <iostream>

#include "include/buchberger.h"
//#include "include/GMPwrapper.hpp"

using namespace ADAAI::LAB04::Poly;

int main()
{
  Monomial<double> ms[] = {
      { 1, 2, new size_t[2] { 0, 1 } },
      { 3, 2, new size_t[2] { 1, 0 } },
      { 0.5, 2, new size_t[2] { 1, 1 } },
      { 1, 2, new size_t[2] { 2, 0 } },
  }; // 3x, y, 0.5xy, x^2

  std::vector<Polynomial<double>> ps = {
      { 2, ms },
      { 2, ms + 1 },
  }; // 3x + y, 0.5xy + 3x

  Polynomial<double> p = {
      4, ms
  }; // x^2 + 0.5xy + 3x + y

  bool RedFlag = false;

  ReduceWrtNSet(p, ps, RedFlag);

//  GroebnerBasis<double> gb = buchberger({ms, ms + 2});

  std::cout << "Monomials:\n";
  for ( const auto& m : ms )
  {
    std::cout << m << " | ";
  }

  std::cout << "\n\nPolynomials:\n";

  for ( const auto& pol : ps )
  {
    std::cout << pol << '\n';
  }

  std::cout << "\nReduced polynomial:\n" << p << '\n';

  return 0;
}
