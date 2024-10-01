#pragma once

namespace ADAAI::Integration::Integrator
{
  struct RHS
  {
    size_t N = 0; // The number of equations

    /// \brief The right-hand side of the system of equations
    virtual void operator()( double current_time, const double* current_state, double* rhs ) const = 0;
  };
} // namespace ADAAI::Integration::Integrator
