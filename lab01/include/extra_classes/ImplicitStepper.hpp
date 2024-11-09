#pragma once

#include <utility>

#include "../../../HSE_NaOM_S2024/integration/intergartor/steppers/BasicTimeStepper.hpp"

#include "../lsu_solvers/solve_lsu.h"

namespace ADAAI::Integration::Integrator::Stepper
{
template<typename RHS, LSU_SOLVERS::LSSolveMethod method = LSU_SOLVERS::LSSolveMethod::GEP>
  requires std::is_base_of_v<ADAAI::Integration::Integrator::RHS, RHS>
class ImplicitStepper : public TimeStepper<RHS>
  {
  double **implicit_matrix;

  public:
    explicit ImplicitStepper(const RHS* rhs )
        : TimeStepper<RHS>( rhs )
    {
      implicit_matrix = new double*[RHS::N];
      for (int i = 0; i < RHS::N; ++i)
      {
        implicit_matrix[i] = new double[RHS::N];
      }
    }

    ~ImplicitStepper()
    {
      for (int i = 0; i < RHS::N; ++i)
      {
        delete[] implicit_matrix[i];
      }
      delete[] implicit_matrix;
    }

    std::pair<double, double>
    operator()( double *current_state, double *next_state, double current_time, double suggested_d_time = 1e-2 ) const override
    {
      this->m_rhs->get_implicit_matrix(implicit_matrix);

      solve_linear_system(RHS::N, implicit_matrix, next_state, current_state, method);

       this->m_rhs->update_RHS_boundary(next_state, current_time);

      return { current_time + suggested_d_time, suggested_d_time };
    }
  }; // class ImplicitStepper
} // namespace ADAAI::Integration::Integrator::Stepper
