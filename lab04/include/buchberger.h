#pragma once

#include <vector>

#include "groebner_basis.h"

namespace ADAAI::LAB04::Poly
{

template <typename T>
void ReduceWrtNSet(Polynomial<T>& p, const std::vector<Polynomial<T>>& N, bool& RedFlag) {
  Monomial<T> * Term = p.leading_term();
  RedFlag = false;

  while (!Term->is_zero()) {
    Polynomial<T> r;
    bool found = false;

    // Find r in N such that LT(r) divides T
    for (const auto& poly : N) {
      if (poly.leading_term()->divides(Term)) {
        r = poly;
        found = true;
        break;
      }
    }

    if (!found) {
      Term = p->next_term(Term);
      continue;
    }

    // Eliminating the term T from p
    Monomial<T> U = p.prev_term(Term);
    p = p - (Term->coeff / r.leading_term().coeff) * (Term / r.leading_term()) * r;
    RedFlag = true;

    if (!U.is_zero()) {
      Term = p.next_term(U);
    } else {
      Term = p.leading_term();
    }
  }

  if (!p.is_zero()) {
    p.normalize();
  }
}

template <typename T>
GroebnerBasis<T> buchberger(const std::vector<Monomial<T>>& polys) {
  GroebnerBasis<T> basis;



  return basis;
}

} // namespace ADAAI::LAB04::Poly
