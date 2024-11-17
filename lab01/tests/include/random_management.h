#pragma once

#include <cstddef>
#include <ctime>
#include <random>

#include "../../include/matrix_work/matrix.h"

namespace ADAAI::LAB01::TEST::RAND
{
template<typename F>
F generate_normally_distributed_double()
{
  static std::default_random_engine  generator( 42 );
  static std::normal_distribution<F> distribution( 0.0, 1.0 );

  return distribution( generator );
}

template<typename F>
void generate_normally_distributed_vector( ADAAI::MATH::Vector<F>& v )
{
  for ( size_t i = 0; i < v.size; ++i )
  {
    v[i] = generate_normally_distributed_double<F>();
  }
}

template<typename F>
void generate_normally_distributed_matrix( ADAAI::MATH::Matrix<F>& m )
{
  for ( size_t i = 0; i < m.rows; ++i )
  {
    for ( size_t j = 0; j < m.cols; ++j )
    {
      m[i][j] = generate_normally_distributed_double<F>();
    }
  }
}
} // namespace ADAAI::LAB01::TEST::RAND
