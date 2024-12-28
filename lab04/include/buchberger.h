#pragma once

#include <queue>
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
    p = p - r.leading_term()->coeff * (Term.get_val() / r.leading_term().get_val()) * r;
    Term.m = p.terms;
    Term.n = p.n_terms;
    RedFlag = true;

    if (Term.prev()->is_zero()) {
      Term = p.leading_term();
    }
  }

  if (!p.is_zero()) {
    p.normalize();
  }
}

template <typename T>
void ReduceWrtEachOther(std::vector<Polynomial<T>>& N) {
  bool NO_Reductions;

  do {
    NO_Reductions = true;

    for (auto it = N.begin(); it != N.end(); ) {
      Polynomial<T> p = *it;
      it = N.erase(it);

      bool RedFlag;
      ReduceWrtNSet(p, N, RedFlag);

      if (!p.is_zero()) {
        N.push_back(p);
      }

      if (RedFlag) {
        NO_Reductions = false;
      }
    }
  } while (!NO_Reductions);
}

template <typename T>
Polynomial<T> SPoly(const Polynomial<T>& f, const Polynomial<T>& g) {
  Monomial<T> f_LT_val = f.leading_term().get_val();
  Monomial<T> g_LT_val = g.leading_term().get_val();

  Monomial<T> LCM = f_LT_val.lcm(g_LT_val);

  Polynomial<T> s = (LCM / f_LT_val) * f - (LCM / g_LT_val) * g;

  return s;
}

template <typename T>
bool SkipDt1(const Polynomial<T>& f, const Polynomial<T>& g, const Monomial<T>& LCM) {
  Monomial<T> f_LT_val = f.leading_term().get_val();
  Monomial<T> g_LT_val = g.leading_term().get_val();

  return LCM == f_LT_val * g_LT_val;
}

template <typename T>
bool SkipDt2(const Polynomial<T>& f, const Polynomial<T>& g, const std::vector<Polynomial<T>>& BG, const Monomial<T>& LCM) {
  for (const auto& h : BG) {
    Monomial<T> h_LT_val = h.leading_term().get_val();
    if (!h_LT_val.divides(LCM))
    {
      return true;
    }
  }

  return false;
}

template <typename T>
void GComplete(const std::vector<Polynomial<T>>& G, std::vector<Polynomial<T>>& BG) {
  BG = G;

  // Normalize each polynomial in BG
  for (auto& f : BG) {
    f.normalize();
  }

  // Reduce BG with respect to each other
  ReduceWrtEachOther(BG);

  // Initialize the queue of critical pairs
  std::priority_queue<std::pair<Polynomial<T>, Polynomial<T>>> CP;
  for (const auto& f : BG) {
    for (const auto& g : BG) {
      if (g < f) {
        CP.push(std::make_pair(f, g));
      }
    }
  }

  // Main critical pairs completion loop
  while (!CP.empty()) {
    auto [f, g] = CP.top();
    CP.pop();

    Monomial<T> LCM = f.leading_term().get_val().lcm(g.leading_term().get_val());

    if (SkipDt2(f, g, BG, LCM)) {
      continue;
    }

    Polynomial<T> s = SPoly(f, g);
    bool RedFlag;
    ReduceWrtNSet(s, BG, RedFlag);

    if (s.is_zero()) {
      continue;
    }

    s.normalize();
    BG.push_back(s);

    for (const auto& h : BG) {
      LCM = s.leading_term().get_val().lcm(h.leading_term().get_val());
      if (!(h == s) && !SkipDt1(h, s, LCM) && !SkipDt2(h, s, BG, LCM)) {
        CP.push(std::make_pair(h, s));
      }
    }
  }

  // Reduce the constructed Gröbner basis
  ReduceWrtEachOther(BG);
}

template <typename T>
GroebnerBasis<T> buchberger(const std::vector<Polynomial<T>>& polys) {
  GroebnerBasis<T> basis;

  std::vector<Polynomial<T>> BG;
  GComplete(polys, BG);

  for (const auto& f : BG) {
    basis.add(f);
  }

  return basis;
}

} // namespace ADAAI::LAB04::Poly
