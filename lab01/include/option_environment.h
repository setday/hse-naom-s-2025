#pragma once

#include <tuple>

#include "auxiliary_functions.h"

namespace ADAAI::LAB01 {
struct OptionEnvironment {
  static constexpr int n = 50; ///< The number of divisions for S_max grid
  static constexpr int m = 50; ///< The number of divisions for V_max grid

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

  double T_max;   ///< Time of option expiration

  /**
   * @brief Constructs an OptionEnvironment object with specified parameters.
   *
   * @param T_max The total time for the simulation (default is 1.0).
   * @param T_tau The time step size
   * @param r The risk-free rate (default is 0.19).
   * @param d The dividend yield (default is 0.01).
   * @param sigma The volatility (default is 0.2).
   * @param epsilon The volatility of volatility (default is 0.3).
   * @param kappa The mean reversion speed (default is 0.1).
   * @param theta The long-term mean level (default is 1).
   * @param beta Sensitivity to market movements (default is 0.9).
   * @param rho Correlation with the market (default is -0.15).
   * @param K The strike price (default is 105).
   * @param S_0 The initial stock price (default is 100).
   * @param V_0 The initial volatility (default is 1).
   *
   * @details The default values are taken from Gazprom stock options.
   */
  explicit OptionEnvironment(double T_max, double T_tau,
                             double r = 0.19, double d = 0.01,
                             double sigma = 0.2, double epsilon = 0.3,
                             double kappa = 0.1, double theta = 1,
                             double beta = 0.9, double rho = -0.15, int K = 105,
                             double S_0 = 100, double V_0 = 1)
      : r(r), d(d), sigma(sigma), epsilon(epsilon), kappa(kappa),
        theta(theta), beta(beta), tau(T_tau), rho(rho), K(K), S_0(S_0),
        V_0(V_0), S_max(-1), V_max(-1), T_max(T_max) {

    // The multipliers for S_max and V_max
    int M_s = 5, M_v = 5;

    std::tie(S_max, V_max) = AUX_FUNC::get_adjusted_s_max_and_v_max(
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
    return AUX_FUNC::payoff_function(h_s_max * j, K);
  }
};
} // namespace ADAAI::LAB01
