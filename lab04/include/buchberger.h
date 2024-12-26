#pragma once

#include <vector>

#include "groebner_basis.h"

namespace ADAAI::LAB04::Poly
{

template <typename T>
void ReduceWrtNSet(Polynomial<T>& p, const std::vector<Polynomial<T>>& N, bool& RedFlag) {
  MonomialPointer<T> Term = p.leading_term();
  RedFlag = false;

  while (!Term->is_zero()) {
    Polynomial<T> r;
    bool found = false;

    // Find r in N such that LT(r) divides T
    for (const auto& poly : N) {
      if (poly.leading_term()->divides(Term.get_val())) {
        r = poly;
        found = true;
        break;
      }
    }

    if (!found) {
      Term = Term.next();
      continue;
    }

    // Eliminating the term T from p
    MonomialPointer<T> U = Term.prev();
    p = p - (Term->coeff / r.leading_term()->coeff) * (Term.get_val() / r.leading_term().get_val()) * r;
    RedFlag = true;

    if (!U->is_zero()) {
      Term = U.next();
    } else {
      Term = p.leading_term();
    }
  }

  if (!p.is_zero()) {
    p.sort();
  }
}

template <typename T>
GroebnerBasis<T> buchberger(const std::vector<Monomial<T>>& polys) {
  GroebnerBasis<T> basis;



  return basis;
}

} // namespace ADAAI::LAB04::Poly
