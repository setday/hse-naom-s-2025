#pragma once

#include "SLAE.h"
#include "gaussian_elimination.h"
#include "../../previous_labs/intergartor/Interator.hpp"

using namespace ADAAI::Integration::Integrator;

namespace ADAAI::LAB01
{

  /**
   * @brief A class that represents a solver for a system of linear algebraic
   * equations (SLAE).
   *
   * This class is responsible for constructing the SLAE, setting initial
   * conditions, and iteratively solving the linear system using a specified
   * number of time steps.
   */
  struct Solver
  {
    static double solve(size_t with_p_steps = 30)
    {
      double t_start = 0.0;
      double t_end = 1.0;
      double t_tau = (t_end - t_start) / (double) with_p_steps;

      OptionEnvironment oe = OptionEnvironment(t_end, t_tau);
      SLAE slae = SLAE(&oe);

      auto *b = new double[slae.N];
      slae.set_initial_RHS(b);

      auto stepper = Stepper::OverrideTimeStepper(&slae);
      auto observer = BlindObserver<SLAE>();

      auto integrator = ODE_Integrator<SLAE, Stepper::OverrideTimeStepper<SLAE>, BlindObserver<SLAE>>(&stepper, &observer);
      integrator(b, b, t_start, t_end, t_tau);
      return slae.get_the_answer(b);
    }
  };

} // namespace ADAAI
