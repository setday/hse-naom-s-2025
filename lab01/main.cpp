#include "include/solver.h"
#include <iostream>

int main() {
  double px = ADAAI::LAB01::Solver<ADAAI::LAB01::SolveMethod::IMPLICIT_GSL>::solve();
  std::cout << "Px of the European Call Option: " << px << '\n';
}