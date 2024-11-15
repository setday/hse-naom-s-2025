#pragma once

#include "../../../HSE_NaOM_S2024/integration/intergartor/Observer.hpp"

namespace ADAAI::Integration::Integrator
{
template<typename RHS>
struct BlindObserver : Observer<RHS>
{
  bool operator()( double        current_time,
                   const double* current_state ) const override
  {
    return true;
  }
};
} // namespace ADAAI::Integration::Integrator
