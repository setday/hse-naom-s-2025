#pragma once

#include <cstddef>
#include <cassert>
#include <algorithm>
#include <functional>
#include <iostream>

#include "monomial.h"

namespace ADAAI::LAB04::Poly
{

constexpr size_t MAX_TERMS = 32;

template <typename T>
struct MonomialPointer {
  Monomial<T>* m;
  size_t n;

  size_t current;

  MonomialPointer(Monomial<T>* m_, size_t n_, size_t current_) : m(m_), n(n_), current(current_) {
    assert(current < n);
    assert(current >= 0);
  }

  Monomial<T> operator*() const {
    return m[current];
  }

  MonomialPointer<T> next() const {
    return MonomialPointer<T>(m, n, current + 1);
  }

  MonomialPointer<T> prev() const {
    return MonomialPointer<T>(m, n, current - 1);
  }
};

template <typename T>
struct Polynomial
{
  size_t n_terms;
  Monomial<T> terms[MAX_TERMS];

  Polynomial() : n_terms(0)
  {
    for (auto & term : terms)
    {
      term = Monomial<T>();
    }
  }

  Polynomial(size_t n, Monomial<T>* t) : n_terms(n)
  {
    for (size_t i = 0; i < n; i++)
    {
      terms[i] = t[i];
    }

    sort();
  }

  MonomialPointer<T> leading_term() const {
    return MonomialPointer<T>(terms, n_terms, 0);
  }

  void append(const Monomial<T>& m)
  {
    terms[n_terms++] = m;
  }

  void sort()
  {
    std::sort(terms, terms + n_terms, std::greater<Monomial<T>>());
  }

  [[nodiscard]] bool is_zero() const
  {
    return n_terms == 0;
  }

  Polynomial<T> operator+(const Polynomial<T>& other) const
  {
    Polynomial<T> p;
    size_t i = 0;
    size_t j = 0;

    while (i < n_terms && j < other.n_terms)
    {
      if (terms[i] > other.terms[j])
      {
        p.append(terms[i++]);
      }
      else if (terms[i] < other.terms[j])
      {
        p.append(other.terms[j++]);
      }
      else
      {
        Monomial<T> m = terms[i] + other.terms[j];
        if (m.coeff != 0)
        {
          p.append(m);
        }
        i++;
        j++;
      }
    }

    while (i < n_terms)
    {
      p.append(terms[i++]);
    }

    while (j < other.n_terms)
    {
      p.append(other.terms[j++]);
    }

    return p;
  }

  Polynomial<T> operator-(const Polynomial<T>& other) const
  {
    Polynomial<T> p;
    size_t i = 0;
    size_t j = 0;

    while (i < n_terms && j < other.n_terms)
    {
      if (terms[i] > other.terms[j])
      {
        p.append(terms[i++]);
      }
      else if (terms[i] < other.terms[j])
      {
        p.append(other.terms[j++].negate());
      }
      else
      {
        Monomial<T> m = terms[i] - other.terms[j];
        if (m.coeff != 0)
        {
          p.append(m);
        }
        i++;
        j++;
      }
    }

    while (i < n_terms)
    {
      p.append(terms[i++]);
    }

    while (j < other.n_terms)
    {
      p.append(other.terms[j++].negate());
    }

    return p;
  }

  Polynomial<T> operator/(const Polynomial<T>& other) const
  {
    Polynomial<T> a = *this;
    Polynomial<T> b = other;
    Polynomial<T> q;

    while (a.n_terms >= b.n_terms)
    {
      size_t n = a.n_terms - b.n_terms;
      size_t p[MAX_VARS] = {0};
      p[0] = n;
      Monomial<T> m(1, 1, p);
      q += m;
      a -= m * b;
    }

    return q;
  }

  Polynomial<T> operator%(const Polynomial<T>& other) const
  {
    Polynomial<T> a = *this;
    Polynomial<T> b = other;

    while (!b.is_zero())
    {
      Polynomial<T> r = a - (a / b) * b;
      a = b;
      b = r;
    }

    return a;
  }

  Polynomial<T> gcd(const Polynomial<T>& other) const
  {
    Polynomial<T> a = *this;
    Polynomial<T> b = other;

    while (!b.is_zero())
    {
      Polynomial<T> r = a % b;
      a = b;
      b = r;
    }

    return a;
  }
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const Polynomial<T>& p)
{
  for (size_t i = 0; i < p.n_terms; i++)
  {
    os << p.terms[i];
    if (i < p.n_terms - 1)
    {
      os << " + ";
    }
  }
  return os;
}

} // namespace ADAAI::LAB04::Poly
