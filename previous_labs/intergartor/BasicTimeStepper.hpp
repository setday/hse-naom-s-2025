#pragma once

#include <utility>

#include "RHS.hpp"

namespace ADAAI::Integration::Integrator::Stepper
{
  template<typename RHS>
  class TimeStepper
  {
  public:
    const RHS* m_rhs;

    explicit TimeStepper( const RHS* rhs )
        : m_rhs( rhs )
    {
    }

    /// \brief The stepper function
    /// \param current_time The current time
    /// \param current_state The current state of the system
    /// \param next_state The next state of the system
    /// \return The next time (current_time + dt) and the delta time

    virtual std::pair<double, double>
    operator()( double *current_state, double *next_state, double current_time, double suggested_d_time ) const = 0;
  }; // class Stepper

  template<typename RHS>
  class OverrideTimeStepper : public TimeStepper<RHS>
  {
  public:
    explicit OverrideTimeStepper( const RHS* rhs )
        : TimeStepper<RHS>( rhs )
    {
    }

    std::pair<double, double>
    operator()( double *current_state, double *next_state, double current_time, double suggested_d_time = 1e-2 ) const override
    {
      ( *this->m_rhs )( current_time, current_state, next_state );

      return { current_time + suggested_d_time, suggested_d_time };
    }
  }; // class DiscreteTimeStepper
} // namespace ADAAI::Integration::Integrator::Stepper
