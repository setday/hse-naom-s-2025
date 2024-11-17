#pragma once

#include <utility>

#include "../../../HSE_NaOM_S2024/integration/intergartor/steppers/BasicTimeStepper.hpp"

#include "../lse_solvers/solve_lse.h"
#include "../matrix_work/matrix.h"

namespace ADAAI::Integration::Integrator::Stepper
{
template<typename RHS,
         LSE_SOLVERS::LSSolveMethod method = LSE_SOLVERS::LSSolveMethod::GEP>
  requires std::is_base_of_v<ADAAI::Integration::Integrator::RHS, RHS>
class ImplicitStepper : public TimeStepper<RHS>
{
  ADAAI::MATH::Matrix<double> implicit_matrix;

public:
  explicit ImplicitStepper( const RHS* rhs )
      : TimeStepper<RHS>( rhs ), implicit_matrix( RHS::N, RHS::N )
  {
  }

  ~ImplicitStepper() = default;

  std::pair<double, double>
  operator()( double* current_state, double* next_state, double current_time, double suggested_d_time = 1e-2 ) const override
  {
    solve_linear_system( RHS::N, this->m_rhs->A, next_state, current_state, method );

    this->m_rhs->update_RHS_boundary( next_state, current_time );

    return { current_time + suggested_d_time, suggested_d_time };
  }
}; // class ImplicitStepper
} // namespace ADAAI::Integration::Integrator::Stepper
