#include <iostream>

#include "include/solver.h"

//#include "tests/test_gaussian_elimination.cpp"
#include "tests/test_openblas.cpp"

int main() {
//  double px = ADAAI::LAB01::Solver<ADAAI::LAB01::SolveMethod::IMPLICIT_GSL>::solve();
//  std::cout << "Px of the European Call Option: " << px << '\n';

  test_lsu_solvers();
}
