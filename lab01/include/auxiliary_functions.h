#pragma once

#include <cmath>

namespace ADAII {

/**
 * @brief Calculates the maximum underlying asset price (S_max).
 *
 * @param T The time period.
 * @param M_s The multiplier for volatility.
 * @param K The strike price.
 * @param S_0 The initial asset price.
 * @param sigma The volatility of the asset.
 * @return The calculated maximum asset price (S_max).
 */
double get_s_max(int T, int M_s, double K, double S_0, double sigma) {
  return std::max(S_0, K) * std::exp(M_s * sigma * std::sqrt(T));
}

/**
 * @brief Calculates the maximum volatility (V_max).
 *
 * @param T The time period.
 * @param theta Long-term mean level.
 * @param M_v The multiplier for volatility.
 * @param eps Volatility of volatility.
 * @return The calculated maximum volatility (V_max).
 */
double get_v_max(int T, double theta, int M_v, double eps) {
  return std::max(1.0, theta) * std::exp(M_v * eps * std::sqrt(T));
}

/**
 * @brief Adjusts S_max and V_max to ensure they lie on the grid.
 *
 * @param T The time period.
 * @param M_s The multiplier for S_max calculation.
 * @param M_v The multiplier for V_max calculation.
 * @param theta Long-term mean level.
 * @param eps Volatility of volatility.
 * @param K The strike price.
 * @param S_0 The initial asset price.
 * @param V_0 The initial volatility.
 * @param sigma The volatility of the asset.
 * @param n The number of divisions for S_max grid.
 * @param m The number of divisions for V_max grid.
 * @return A pair containing the adjusted S_max and V_max.
 */
std::pair<double, double> get_adjusted_s_max_and_v_max(int T, int M_s, int M_v,
                                                       double theta, double eps,
                                                       double K, double S_0,
                                                       double V_0, double sigma,
                                                       int n, int m) {
  // Step 1: Compute the initial estimate for S_max and V_max.
  double S_max = get_s_max(T, M_s, K, S_0, sigma);
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

  return {S_max, V_max};
}

/**
 * @brief Computes the payoff for a given asset price.
 *
 * @param S The underlying asset price.
 * @param K The strike price.
 * @return The calculated payoff, which is max(S - K, 0).
 */
double payoff_function(double S, double K) { return std::max(S - K, 0.0); }
} // namespace ADAII