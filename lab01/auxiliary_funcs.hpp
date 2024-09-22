#pragma once

#include <cmath> // Include for std::exp and std::sqrt

namespace ADAII {

double get_s_max(int T, int M_s, int K, double S_0, double sigma) {
  return std::max(S_0, K) * std::exp(M_s * sigma * std::sqrt(T));
}

double get_v_max(int T, double theta, int M_v, double eps) {
  return std::max(1.0, theta) * std::exp(M_v * eps * std::sqrt(T));
}

std::pair<double, double> adjust_s_max_and_v_max(int T, int M_s, int M_v,
                                                 double theta, double eps,
                                                 int K, double S_0, double V_0,
                                                 int n = 100, int m = 100) {
  // Step 1: Compute the initial estimate for S_max and V_max.
  double S_max = get_s_max(T, M_s, K, S_0);
  double V_max = get_v_max(T, theta, M_v, eps);

  // Step 2: Adjust S_max and V_max s.t. S(0) and V(0) lie exactly on the grid.
  // Derivation:
  //            We want S(0) to lie exactly on the grid.
  //            This means that S(0) = S_max * (j / n) for some natural j <= n.
  //            Hence, we can compute j as [n * (S(0) / S_max)] and then
  //            use this value to compute S_max. Now, S(0) = S_max * (j / n).
  //            The same procedure applies to V_max.

  int j = static_cast<int>(n * S_0 / S_max);
  S_max = S_0 * n / j;

  int i = static_cast<int>(m * V_0 / V_max);
  V_max = V_0 * m / i;

  return {S_max, V_max}; // Do we need to return new steps (h_v, h_s) or (i, j) or something else???
}
} // namespace ADAII
