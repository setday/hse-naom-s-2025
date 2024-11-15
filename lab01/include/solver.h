#pragma once

#include "../../HSE_NaOM_S2024/integration/intergartor/Interator.hpp"
#include "../../HSE_NaOM_S2024/integration/intergartor/steppers/BasicTimeStepper.hpp"

#include "extra_classes/BlindObserver.hpp"
#include "extra_classes/ImplicitStepper.hpp"
#include "lse_solvers/solve_lse.h"
#include "SLAE.h"

using namespace ADAAI::Integration::Integrator;

namespace ADAAI::LAB01
{

enum class SolveMethod
{
  EXPLICIT,
  IMPLICIT_GSL,
  IMPLICIT_OPENBLAS,
  IMPLICIT_GAUSSIAN_ELIMINATION,
};

/**
 * @brief A class that represents a solver for a system of linear algebraic
 * equations (SLAE).
 *
 * This class is responsible for constructing the SLAE, setting initial
 * conditions, and iteratively solving the linear system using a specified
 * number of time steps.
 */
template<SolveMethod method = SolveMethod::IMPLICIT_GAUSSIAN_ELIMINATION>
struct Solver
{
  static double solve( size_t with_p_steps = 30 )
  {
    double t_start = 0.0;
    double t_end   = 1.0;
    double t_tau   = ( t_end - t_start ) / ( double ) with_p_steps;

    OptionEnvironment oe   = OptionEnvironment( t_end, t_tau );
    SLAE              slae = SLAE( &oe );

    auto* b = new double[SLAE::N];
    slae.set_initial_RHS( b );

    if constexpr ( method == SolveMethod::EXPLICIT )
    {
      auto stepper =
          ADAAI::Integration::Integrator::Stepper::RFK45_TimeStepper<SLAE>(
              &slae );
      auto observer = BlindObserver<SLAE>();

      auto integrator = ODE_Integrator<
          SLAE,
          ADAAI::Integration::Integrator::Stepper::RFK45_TimeStepper<SLAE>,
          BlindObserver<SLAE>>( &stepper, &observer );
      integrator( b, b, t_start, t_end, t_tau );
    }
    else if constexpr ( method == SolveMethod::IMPLICIT_GSL )
    {
      auto stepper =
          Stepper::ImplicitStepper<SLAE, LSE_SOLVERS::LSSolveMethod::GSL>(
              &slae );
      auto observer = BlindObserver<SLAE>();

      auto integrator = ODE_Integrator<
          SLAE, Stepper::ImplicitStepper<SLAE, LSE_SOLVERS::LSSolveMethod::GSL>,
          BlindObserver<SLAE>>( &stepper, &observer );
      integrator( b, b, t_start, t_end, t_tau );
    }
    else if constexpr ( method == SolveMethod::IMPLICIT_OPENBLAS )
    {
      auto stepper =
          Stepper::ImplicitStepper<SLAE, LSE_SOLVERS::LSSolveMethod::OPENBLAS>(
              &slae );
      auto observer = BlindObserver<SLAE>();

      auto integrator = ODE_Integrator<
          SLAE,
          Stepper::ImplicitStepper<SLAE, LSE_SOLVERS::LSSolveMethod::OPENBLAS>,
          BlindObserver<SLAE>>( &stepper, &observer );
      integrator( b, b, t_start, t_end, t_tau );
    }
    else if constexpr ( method == SolveMethod::IMPLICIT_GAUSSIAN_ELIMINATION )
    {
      auto stepper =
          Stepper::ImplicitStepper<SLAE, LSE_SOLVERS::LSSolveMethod::GEP>(
              &slae );
      auto observer = BlindObserver<SLAE>();

      auto integrator = ODE_Integrator<
          SLAE, Stepper::ImplicitStepper<SLAE, LSE_SOLVERS::LSSolveMethod::GEP>,
          BlindObserver<SLAE>>( &stepper, &observer );
      integrator( b, b, t_start, t_end, t_tau );
    }

    return slae.get_the_answer( b );
  }
};

} // namespace ADAAI::LAB01
