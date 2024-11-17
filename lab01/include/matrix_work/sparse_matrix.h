#pragma once

#include <stdexcept>

#include "matrix.h"

namespace ADAAI::MATH
{

template<typename T>
struct SparceElement
{
  size_t row;
  T      value;
};

template<typename T>
struct SparceRowView
{
  SparceElement<T>* data = nullptr;
  size_t            size = 0;

  explicit SparceRowView( SparceElement<T>* _data, size_t _size )
      : data( _data ), size( _size )
  {
  }

  ~SparceRowView() = default;

  SparceElement<T>& operator[]( size_t i )
  {
    for ( size_t k = 0; k < size; ++k )
    {
      if ( data[k].row == i )
      {
        return data[k];
      }
    }
    throw std::out_of_range( "Element not found" );
  }

  T operator[]( size_t i ) const
  {
    for ( size_t k = 0; k < size; ++k )
    {
      if ( data[k].row == i )
      {
        return data[k].value;
      }
    }
    return 0;
  }
};

template<typename T>
struct SparseMatrix
{
  SparceElement<T>** data      = nullptr;
  size_t*            row_sizes = nullptr;
  size_t             rows      = 0;

  SparseMatrix() = default;

  explicit SparseMatrix( size_t _rows )
      : rows( _rows )
  {
    allocate( rows );
  }

  explicit SparseMatrix( const Matrix<T>& matrix )
  {
    rows = matrix.rows;
    allocate( rows );

    for ( size_t i = 0; i < rows; ++i )
    {
      row_sizes[i] = 0;
      for ( size_t j = 0; j < matrix.cols; ++j )
      {
        row_sizes[i] += matrix[i][j] != 0;
      }

      data[i]  = new SparceElement<T>[row_sizes[i]];
      size_t k = 0;

      for ( size_t j = 0; j < matrix.cols; ++j )
      {
        if ( matrix[i][j] != 0 )
        {
          data[i][k] = { j, matrix[i][j] };
          ++k;
        }
      }
    }
  }

  SparseMatrix( const SparseMatrix& )            = delete;
  SparseMatrix& operator=( const SparseMatrix& ) = delete;

  ~SparseMatrix()
  {
    deallocate();
  }

  /**
   * @brief Fill the matrix with data from a dense matrix
   *
   * @param matrix Dense matrix to fill the sparse matrix with
   */
  void from_dense( const Matrix<T>& matrix )
  {
    if ( rows != 0 )
    {
      throw std::runtime_error( "Matrix already have data. Create a clean matrix to fill it with new data" );
    }

    rows = matrix.rows;
    allocate( rows );

    for ( size_t i = 0; i < rows; ++i )
    {
      row_sizes[i] = 0;
      for ( size_t j = 0; j < matrix.cols; ++j )
      {
        row_sizes[i] += matrix[i][j] != 0;
      }

      data[i]  = new SparceElement<T>[row_sizes[i]];
      size_t k = 0;

      for ( size_t j = 0; j < matrix.cols; ++j )
      {
        if ( matrix[i][j] != 0 )
        {
          data[i][k] = { j, matrix[i][j] };
          ++k;
        }
      }
    }
  }

  /**
   * @brief Access element of the matrix
   *
   * @param i Row index
   * @return Reference to the element
   *
   * @example A[0][0] = 5;
   */
  SparceRowView<T>& operator[]( size_t i )
  {
    return SparceRowView<T>( data[i], row_sizes[i] );
  }

  /**
   * @brief Transform ~b = D^-1 @ b, where D = diag(A)
   */
  void transform_direct( Vector<T>& b )
  {
    for ( size_t i = 0; i < rows; ++i )
    {
      for ( size_t j = 0; j < row_sizes[i]; ++j )
      {
        if ( data[i][j].row == i )
        {
          b[i] /= data[i][j].value;
          break;
        }
      }
    }
  }

  /**
   * @brief Transform A = ~B = D^-1 @ B, where B = A - D, D = diag(A)
   */
  void transform_iterative()
  {
    for ( size_t i = 0; i < rows; ++i )
    {
      T diagonal_element = 1;
      for ( size_t j = 0; j < row_sizes[i]; ++j )
      {
        if ( data[i][j].row == i )
        {
          diagonal_element = data[i][j].value;
          data[i][j].value = 0; /// This is too expensive to erase an element so we just set it to 0
          break;
        }
      }
      for ( size_t j = 0; j < row_sizes[i]; ++j )
      {
        data[i][j].value /= diagonal_element;
      }
    }
  }

private:
  void allocate( size_t _rows )
  {
    row_sizes = new size_t[_rows];
    data      = new SparceElement<T>*[_rows];
  }

  void deallocate()
  {
    delete[] row_sizes;
    delete[] data;
  }
};

/**
 * @brief Multiply matrix by vector (result = A @ v)
 *
 * @param A Matrix to multiply by
 * @param v Vector to multiply by
 * @param result Result of multiplication
 * @param keep_result If set to false, the result vector will be zeroed before the multiplication
 */
template<typename T>
void vec_mul( SparseMatrix<T>& A, Vector<T>& v, Vector<T>& result, bool keep_result = true )
{
  for ( size_t i = 0; i < A.rows; ++i )
  {
    if ( !keep_result )
    {
      result[i] = 0;
    }
    for ( size_t j = 0; j < A.row_sizes[i]; ++j )
    {
      result[i] += A.data[i][j].value * v[A.data[i][j].row];
    }
  }
}

/**
 * @brief Multiply vector by matrix (result = v @ A)
 *
 * @param v Vector to multiply by
 * @param A Matrix to multiply by
 * @param result Result of multiplication
 * @param keep_result If set to false, the result vector will be zeroed before the multiplication
 */
template<typename T>
void vec_mul( Vector<T>& v, SparseMatrix<T>& A, Vector<T>& result, bool keep_result = true )
{
  for ( size_t i = 0; i < A.rows; ++i )
  {
    if ( !keep_result )
    {
      result[i] = 0;
    }
    for ( size_t j = 0; j < A.row_sizes[i]; ++j )
    {
      result[A.data[i][j].row] += A.data[i][j].value * v[i];
    }
  }
}
} // namespace ADAAI::MATH