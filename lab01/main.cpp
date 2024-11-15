#include <iostream>

#include "include/solver.h"

#include "tests/test_lse_solvers.cpp"

int main()
{
  //  double px =
  //  ADAAI::LAB01::Solver<ADAAI::LAB01::SolveMethod::IMPLICIT_GSL>::solve();
  //  std::cout << "Px of the European Call Option: " << px << '\n';

  test_lse_solvers<double, ADAAI::LSE_SOLVERS::LSSolveMethod::GEP>();
}
