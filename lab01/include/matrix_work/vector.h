#pragma once

#include <cmath>
#include <cstdlib>

namespace ADAAI::MATH
{

template<typename T>
struct Vector
{
  T*     data;
  size_t size;

  Vector( size_t _size )
      : size( _size )
  {
    allocate( size );
  }

  Vector( const Vector& )            = delete;
  Vector& operator=( const Vector& ) = delete;

  ~Vector()
  {
    deallocate();
  }

  /**
   * @brief Cast to the pointer to the data
   *
   * @return Pointer to the data
   *
   * @example double* x = A;
   */
  explicit operator T*() const
  {
    return data;
  }

  /**
   * @brief Access element of the vector
   * 
   * @param i Index
   * @return Reference to the element
   * 
   * @example A[0] = 5;
   */
  T& operator[]( size_t i )
  {
    return data[i];
  }

  /**
   * @brief Access element of the vector
   * 
   * @param i Index
   * @return Reference to the element
   * 
   * @example double x = A[0];
   */
  T operator[]( size_t i ) const
  {
    return data[i];
  }

  /**
   * @brief Add vector to the current vector
   *
   * @param v Vector to add
   * @return Reference to the current vector
   *
   * @example A += B;
   */
  Vector& operator+=( Vector& v )
  {
    for ( size_t i = 0; i < size; ++i )
    {
      data[i] += v[i];
    }
    return *this;
  }

  /**
   * @brief Calculate the norm of the vector
   * 
   * @return Norm of the vector
   * 
   * @example double len_squared = v % v;
   */
  T operator%( Vector& v )
  {
    T result = 0;
    for ( size_t i = 0; i < size; ++i )
    {
      result += data[i] * v[i];
    }
    return result;
  }

private:
  void allocate( size_t _size )
  {
    data = new T[_size];
  }

  void deallocate()
  {
    delete[] data;
  }
};

} // namespace ADAAI::MATH