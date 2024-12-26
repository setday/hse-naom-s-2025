#pragma once

#include <cstddef>
#include <cassert>
#include <iostream>

namespace ADAAI::LAB04::Poly
{

constexpr size_t MAX_VARS = 32;

template <typename T>
struct Monomial
{
  T coeff;
  size_t powers[MAX_VARS]{};

  Monomial() : coeff() { }

  Monomial(T c, size_t n, const size_t* p) : coeff(c)
  {
    for (size_t i = 0; i < n; i++)
    {
      powers[i] = p[i];
    }
  }

  T get_coeff() const
  {
    return coeff;
  }

  [[nodiscard]] bool is_zero() const
  {
    return coeff == 0;
  }

  Monomial<T> operator+(const Monomial<T>& m) const
  {
    for (size_t i = 0; i < MAX_VARS; i++)
    {
      assert(powers[i] == m.powers[i]);
    }

    return Monomial<T>(coeff + m.coeff, MAX_VARS, powers);
  }

  Monomial<T> operator-(const Monomial<T>& m) const
  {
    for (size_t i = 0; i < MAX_VARS; i++)
    {
      assert(powers[i] == m.powers[i]);
    }

    return Monomial<T>(coeff - m.coeff, MAX_VARS, powers);
  }

  Monomial<T> operator*(const Monomial<T>& m) const
  {
    T c = coeff * m.coeff;
    size_t p[MAX_VARS] = {0};

    for (size_t i = 0; i < MAX_VARS; i++)
    {
      p[i] = powers[i] + m.powers[i];
    }

    return Monomial<T>(c, MAX_VARS, p);
  }

  Monomial<T> operator/(const Monomial<T>& m) const
  {
    T c = coeff / m.coeff;
    size_t p[MAX_VARS] = {0};

    for (size_t i = 0; i < MAX_VARS; i++)
    {
      assert(powers[i] >= m.powers[i]);
      p[i] = powers[i] - m.powers[i];
    }

    return Monomial<T>(c, MAX_VARS, p);
  }

  Monomial<T> negate() const
  {
    return Monomial<T>(-coeff, MAX_VARS, powers);
  }

  bool divides(const Monomial<T>& m) const
  {
    for (size_t i = 0; i < MAX_VARS; i++)
    {
      if (powers[i] < m.powers[i])
      {
        return false;
      }
    }
    return true;
  }
};

template <typename T>
Monomial<T> operator*(T c, const Monomial<T>& m)
{
  return Monomial<T>(c * m.coeff, MAX_VARS, m.powers);
}

template <typename T>
bool operator>(const Monomial<T>& zeta, const Monomial<T>& nu)
{
  // L (Lexicographic order)
  for (size_t i = 0; i < MAX_VARS; i++)
  {
    if (zeta.powers[i] == nu.powers[i])
    {
      continue;
    }
    return zeta.powers[i] > nu.powers[i];
  }
  return false;
}

template <typename T>
bool operator<(const Monomial<T>& zeta, const Monomial<T>& nu)
{
  return nu > zeta;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Monomial<T>& m)
{
  os << m.coeff;
  for (size_t i = 0; i < MAX_VARS; i++)
  {
    if (m.powers[i] == 0)
    {
      continue;
    }
    os << "x_" << i << "^" << m.powers[i];
  }
  return os;
}

} // namespace ADAAI::LAB04::Poly
