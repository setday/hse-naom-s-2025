#pragma once

#include "RHS.hpp"

namespace ADAAI::Integration::Integrator
{
  template<typename RHS>
  struct Observer
  {
    virtual bool operator()( double current_time, const double *current_state ) const = 0;
  };

  template<typename RHS>
  struct BlindObserver : Observer<RHS>
  {
    bool operator()( double current_time, const double *current_state ) const override {
      return true;
    }
  };
} // namespace ADAAI::Integration::Integrator
