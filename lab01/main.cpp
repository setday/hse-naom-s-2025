#include "include/solver.h"
#include <iostream>

int main() {
  ADAII::Solver solver;
  double px = solver.solve();
  std::cout << "Px of the European Call Option: " << px << '\n';
}