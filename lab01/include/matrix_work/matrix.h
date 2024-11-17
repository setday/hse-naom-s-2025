#pragma once

#include "vector.h"

namespace ADAAI::MATH
{

template<typename T>
struct Matrix
{
  T**    data;
  size_t rows;
  size_t cols;

  Matrix( size_t _rows, size_t _cols )
      : rows( _rows ), cols( _cols )
  {
    allocate( rows, cols );
  }

  Matrix( const Matrix& other )
      : rows( other.rows ), cols( other.cols )
  {
    allocate( rows, cols );
    for ( size_t i = 0; i < rows; ++i )
    {
      for ( size_t j = 0; j < cols; ++j )
      {
        data[i][j] = other.data[i][j];
      }
    }
  }

  Matrix& operator=( const Matrix& ) = delete;

  ~Matrix()
  {
    deallocate();
  }

  /**
   * @brief Cast to the pointer to the data
   *
   * @return Pointer to the data
   *
   * @example double** A = matrix;
   */
  explicit operator T**() const
  {
    return data;
  }

  /**
   * @brief Access element of the matrix
   * 
   * @param i Row index
   * @return Reference to the element
   * 
   * @example A[0][1] = 5;
   */
  T* operator[]( size_t i ) const
  {
    return data[i];
  }

  /**
   * @brief Transform ~b = b
   *
   * @note This function is used for compatibility with the SparseMatrix class
   */
  void transform_direct( Vector<T>& b )
  {
    /// Do nothing
  }

  /**
   * @brief Transform A = ~B = A - I
   */
  void transform_iterative()
  {
    for ( size_t i = 0; i < rows; ++i )
    {
      data[i][i] -= 1;
    }
  }

private:
  void allocate( size_t _rows, size_t _cols )
  {
    data = new T*[_rows];
    for ( size_t i = 0; i < _rows; ++i )
    {
      data[i] = new T[_cols];
    }
  }

  void deallocate()
  {
    for ( size_t i = 0; i < rows; ++i )
    {
      delete[] data[i];
    }
    delete[] data;
  }
};

/**
 * @brief Multiply matrix by vector (result = A @ v)
 *
 * @param A Matrix to multiply by
 * @param v Vector to multiply by
 * @param result Result of multiplication
 * @param keep_result If true, the result vector will be added to the existing values in the result vector
 */
template<typename T>
void vec_mul( Matrix<T>& A, Vector<T>& v, Vector<T>& result, bool keep_result = true )
{
  for ( size_t i = 0; i < A.rows; ++i )
  {
    if ( !keep_result )
    {
      result[i] = 0;
    }
    for ( size_t j = 0; j < A.cols; ++j )
    {
      result[i] += A.data[i][j] * v[j];
    }
  }
}

/**
 * @brief Multiply vector by matrix (result = v @ A)
 *
 * @param v Vector to multiply by
 * @param A Matrix to multiply by
 * @param result Result of multiplication
 * @param keep_result If true, the result vector will be added to the existing values in the result vector
 */
template<typename T>
void vec_mul( Vector<T>& v, Matrix<T>& A, Vector<T>& result, bool keep_result = true )
{
  for ( size_t i = 0; i < A.cols; ++i )
  {
    if ( !keep_result )
    {
      result[i] = 0;
    }
    for ( size_t j = 0; j < A.rows; ++j )
    {
      result[i] += A.data[j][i] * v[j];
    }
  }
}
} // namespace ADAAI::MATH