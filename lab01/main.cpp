#include "include/solver.h"
#include <iostream>

int main() {
  double px = ADAAI::LAB01::Solver::solve();
  std::cout << "Px of the European Call Option: " << px << '\n';
}