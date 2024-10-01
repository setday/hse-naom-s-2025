#pragma once

#include "SLAE.h"
#include "gaussian_elimination.h"
#include "../../previous_labs/intergartor/Interator.hpp"

using namespace ADAAI::Integration::Integrator;

namespace ADAAI
{

  /**
   * @brief A class that represents a solver for a system of linear algebraic
   * equations (SLAE).
   *
   * This class is responsible for constructing the SLAE, setting initial
   * conditions, and iteratively solving the linear system using a specified
   * number of time steps.
   */
  class Solver
  {
  public:
    /**
     * @brief Constructs a Solver instance with a specified number of time steps.
     *
     * @param p The number of time steps (default is 30).
     */
    Solver(int p = 30) : p(p) {}

    double solve()
    {
      SLAE slae = SLAE(p);
      auto *b = new double[slae.N];
      slae.set_initial_RHS(b);

      auto stepper = Stepper::OverrideTimeStepper(&slae);
      auto observer = BlindObserver<SLAE>();

      auto integrator = ODE_Integrator<SLAE, Stepper::OverrideTimeStepper<SLAE>, BlindObserver<SLAE>>(&stepper, &observer);
      integrator(b, b, 0, 30, 1);
      return slae.get_the_answer(b);
    }

  private:
    int p; ///< The number of time steps for the solver.
  };

} // namespace ADAAI
