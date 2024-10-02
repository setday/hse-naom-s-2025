#pragma once

#include <tuple>

#include "auxiliary_functions.h"

namespace ADAAI::LAB01 {
struct OptionEnvironment {
  int n; ///< The number of divisions for S_max grid
  int m; ///< The number of divisions for V_max grid

  double r;       ///< Risk-free rate
  double d;       ///< Dividend yield
  double sigma;   ///< Volatility
  double epsilon; ///< Volatility of volatility
  double kappa;   ///< Mean reversion speed
  double theta;   ///< Long-term mean level
  double beta;    ///< Sensitivity to market movements
  double tau;     ///< Time step size
  double rho;     ///< Correlation with the market
  int K;          ///< Strike price

  double S_0;     ///< Initial stock price
  double V_0;     ///< Initial volatility
  double S_max;   ///< Maximum asset price
  double V_max;   ///< Maximum volatility
  double h_s_max; ///< Step size in asset price
  double h_v_max; ///< Step size in volatility

  /**
   * @brief Constructs an OptionEnvironment object with specified parameters.
   *
   * @param T_max The total time for the simulation (default is 1.0).
   * @param T_tau The time step size
   * @param n The number of divisions for the asset price grid (default is 50).
   * @param m The number of divisions for the volatility grid (default is 50).
   * @param r The risk-free rate (default is 0.01).
   * @param d The dividend yield (default is 0.02).
   * @param sigma The volatility (default is 0.2).
   * @param epsilon The volatility of volatility (default is 0.2).
   * @param kappa The mean reversion speed (default is 0.5).
   * @param theta The long-term mean level (default is 0.2).
   * @param beta Sensitivity to market movements (default is 1).
   * @param rho Correlation with the market (default is 0.5).
   * @param K The strike price (default is 100).
   * @param M_s The multiplier for S_max (default is 5).
   * @param M_v The multiplier for V_max (default is 5).
   * @param S_0 The initial stock price (default is 100).
   * @param V_0 The initial volatility (default is 1).
   */
  explicit OptionEnvironment(double T_max, double T_tau, int n = 50, int m = 50,
                             double r = 0.01, double d = 0.02,
                             double sigma = 0.2, double epsilon = 0.2,
                             double kappa = 0.5, double theta = 0.2,
                             double beta = 1, double rho = 0.5, int K = 100,
                             int M_s = 5, int M_v = 5, double S_0 = 100,
                             double V_0 = 1)
      : n(n), m(m), r(r), d(d), sigma(sigma), epsilon(epsilon), kappa(kappa),
        theta(theta), beta(beta), tau(T_tau), rho(rho), K(K), S_0(S_0),
        V_0(V_0) {

    std::tie(S_max, V_max) = get_adjusted_s_max_and_v_max(
        T_max, M_s, M_v, theta, epsilon, K, S_0, V_0, sigma, n, m);

    h_s_max = S_max / n;
    h_v_max = V_max / m;
  }

  [[nodiscard]] int get_closest_S_index(double S) const {
    return std::lround(S / h_s_max);
  }

  [[nodiscard]] int get_closest_V_index(double V) const {
    return std::lround(V / h_v_max);
  }

  [[nodiscard]] double get_payoff_for_S_j(int j) const {
    return payoff_function(h_s_max * j, K);
  }
};
} // namespace ADAAI::LAB01
