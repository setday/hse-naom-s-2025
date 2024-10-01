#pragma once

#include "SLAE.h"
#include "gaussian_elimination.h"

namespace ADAII {

/**
 * @brief A class that represents a solver for a system of linear algebraic equations (SLAE).
 * 
 * This class is responsible for constructing the SLAE, setting initial conditions,
 * and iteratively solving the linear system using a specified number of time steps.
 */
class Solver {
public:
  /**
   * @brief Constructs a Solver instance with a specified number of time steps.
   * 
   * @param p The number of time steps (default is 30).
   */
  Solver(int p = 30) : p(p) {}

  double solve() {
    SLAE slae = SLAE(p);
    slae.construct_SLAE();
    slae.set_initial_RHS();
    slae.update_left_and_right_boundaries(p);

    for (int step = p; step >= 0; step--) {
      std::cout << "step: " << p - step << "/" << p << '\n';
      solve_linear_system(slae);
      slae.update_RHS();
      slae.update_left_and_right_boundaries(step);
    }
    return slae.get_the_answer();
  }

private:
  int p; ///< The number of time steps for the solver.
};

} // namespace ADAII
