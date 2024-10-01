#pragma once

#include <iostream>

#include "BasicTimeStepper.hpp"
#include "Observer.hpp"

namespace ADAAI::Integration::Integrator
{
  /// \brief The ODE integrator
  /// \tparam RHS The right-hand side of the system of equations
  /// \tparam TS The time stepper
  /// \tparam RHS_O The observer

  template<typename RHS_I = RHS, typename TS = Stepper::TimeStepper<RHS_I>, typename RHS_O = BlindObserver<RHS_I>>
    requires std::is_base_of_v<RHS, RHS_I> && std::is_base_of_v<Observer<RHS_I>, RHS_O> && std::is_base_of_v<Stepper::TimeStepper<RHS_I>, TS>
  class ODE_Integrator
  {
    const TS*    m_stepper;
    const RHS_O* m_observer;

  public:
    ODE_Integrator( const TS* stepper, const RHS_O* observer )
        : m_stepper( stepper ), m_observer( observer )
    {
    }

    /// \brief The integrator function
    /// \param state_start The initial state of the system
    /// \param state_end The final state of the system
    /// \param t_start The initial time
    /// \param t_end The final time
    /// \return The time of the final state
    double operator()( const double *state_start, double *state_end, double t_start = 0.0, double t_end = 2e3, double suggested_dt = 1e-2 ) const
    {
      size_t N = m_stepper->m_rhs->N;

      double current_time = t_start;
      auto *current_state = new double[N];
      auto *next_state = new double[N];

      for ( size_t i = 0; i < N; ++i )
      {
        current_state[i] = state_start[i];
      }

      int percent_count = 0;

      while ( current_time < t_end )
      {
        if ( !( *m_observer )( current_time, current_state ) )
        {
          break;
        }

        auto [next_time, dt] = ( *m_stepper )( current_state, next_state, current_time, suggested_dt );

        if ( next_time > t_end )
        {
          dt        = t_end - current_time;
          next_time = t_end;
        }

        current_time = next_time;
        for ( size_t i = 0; i < N; ++i )
        {
          current_state[i] = next_state[i];
        }

        if ( current_time / t_end * 100.0 > percent_count + 5 )
        {
          percent_count += 5;
          std::cout << percent_count << "%" << std::endl;
        }
      }

      for ( size_t i = 0; i < N; ++i )
      {
        state_end[i] = current_state[i];
      }

      return current_time;
    }
  };
} // namespace ADAAI::Integration::Integrator
