#pragma once

#include <cmath>
#include <iostream>

#include "../../previous_labs/intergartor/RHS.hpp"

#include "auxiliary_functions.h"
#include "gaussian_elimination.h"

namespace ADAAI {
class SLAE : Integration::Integrator::RHS {
  /*
  This class defines a system of linear algebraic equations (Ax = b) that will
  be solved using Gaussian elimination. The system arises from the numerical
  solution of a SDE (stochastic differential equation) using an implicit method
  (Euler).
  */

public:
  /**
   * @brief Constructs an SLAE object with specified parameters.
   *
   * Initializes the internal matrices and vectors based on the given
   * parameters, setting up the grid sizes and time step.
   *
   * @param p The number of time steps.
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
   * @param T The total time for the simulation (default is 1).
   */
  explicit SLAE(int p, int n = 50, int m = 50, double r = 0.01, double d = 0.02,
       double sigma = 0.2, double epsilon = 0.2, double kappa = 0.5,
       double theta = 0.2, double beta = 1, double rho = 0.5, int K = 100,
       int M_s = 5, int M_v = 5, double S_0 = 100, double V_0 = 1, int T = 1)
      : n(n), m(m), r(r), d(d), sigma(sigma), epsilon(epsilon), kappa(kappa),
        theta(theta), beta(beta), rho(rho), p(p), K(K), S_0(S_0), V_0(V_0) {
    N = (n + 1) * (m + 1);

    A = new double *[N];
    for (size_t i = 0; i < N; ++i) {
      A[i] = new double[N]();
    }

    tau = static_cast<double>(T) / p; // time step
    auto [S_max, V_max] = get_adjusted_s_max_and_v_max(
        T, M_s, M_v, theta, epsilon, K, S_0, V_0, sigma, n, m);
    h_s = S_max / n;
    h_v = V_max / m;

    construct_SLAE();
  }

  /**
   * @brief Sets the initial conditions for the right-hand side vector (b).
   */
  void set_initial_RHS(double* rhs) { // utilize initial condition @t=T
    for (int i = 0; i <= m; i++) {
      int i0 = get_Q_index(i, 0);
      rhs[i0] = 0;
    }
    for (int i = 1; i <= m - 1; i++) {
      for (int j = 1; j <= n - 1; j++) {
        double S_j = h_s * j;
        int ij = get_Q_index(i, j);
        rhs[ij] = payoff_function(S_j, K);
      }
    }
    for (int j = 1; j <= n - 1; j++) {
      int mj = get_Q_index(m, j);
      rhs[mj] = 0;
    }

    int k = p;
    // Lower Boundary
    double const7 = (p - k + 1) * tau;
    double const8 = K * std::exp((d - r) * const7) / h_s;
    double const9 = std::exp(-d * const7) * h_s;
    double const10 = std::exp(-r * const7) * K;
    for (int j = 0; j <= n; j++) {
      int oj = get_Q_index(0, j);
      if (j >= const8) {
        rhs[oj] = const9 * j - const10;
      } else {
        rhs[oj] = 0;
      }
    }

    // Right Boundary
    double const11 = const9 * n - const10;
    for (int i = 1; i <= m; i++) {
      int i_n = get_Q_index(i, n);
      rhs[i_n] = const11;
    }
  }

  /**
   * @brief Retrieves the solution for the given initial stock price and volatility.
   *
   * Computes the indices corresponding to the current stock price (S_0)
   * and volatility (V_0) in the solution vector and returns the corresponding
   * value from the solution.
   *
   * @return The solution corresponding to the initial conditions.
   */
  double get_the_answer(double const *rhs) {
    int i = 0;
    int j = 0;
    double S_i = 0;
    double V_j = 0;
    for (; i <= m; i++) {
      if (S_i >= S_0) {
        break;
      }
      S_i += h_s;
    }
    for (; j <= n; j++) {
      if (V_j >= V_0) {
        break;
      }
      V_j += h_v;
    }

    return rhs[get_Q_index(i, j)];
  }

  /**
   * @brief Destructor to clean up dynamically allocated memory.
   *
   * Frees the memory allocated for the solution vector, RHS vector, and the
   * coefficient matrix.
   */
  ~SLAE() {
    for (size_t i = 0; i < N; ++i) {
      delete[] A[i];
    }
    delete[] A;
  }

  void operator()( double current_time, const double* current_state, double* rhs ) const override {
    solve_linear_system(N, A, rhs, current_state);
  }

private:
  /**
   * @brief Constructs the coefficient matrix (A) and the right-hand side vector
   * (b), BUT only partly!
   *
   * Populates the internal matrix and vector based on the specified
   * parameters and the boundary conditions derived from the problem.
   */
  void construct_SLAE() {
    // Inner matrix construction

    // We define many constants which don't have any sense, but
    // they are used many times in different places. It is only a matter of
    // optimization.
    double const1 = tau * (r - d) / 2;
    double const2 = tau * sigma * sigma / 2 * h_v * std::pow(h_s, 2 * beta - 2);
    double const4 = tau * epsilon * epsilon / (2 * h_v);
    double const5 = tau * sigma * epsilon * rho * std::pow(h_s, beta - 1) / 4;

    for (int i = 1; i <= m - 1; i++) {

      double const3 = const2 * i;
      double const6 = const5 * i;

      for (int j = 1; j <= n - 1; j++) {
        int ij = get_Q_index(i, j);
        // The first term in LSH
        A[ij][ij] += 1;

        // The second term in LHS (operator L: convective + diffusive terms)
        int ipj = get_Q_index(i + 1, j);
        int ijp = get_Q_index(i, j + 1);
        int imj = get_Q_index(i - 1, j);
        int ijm = get_Q_index(i, j - 1);
//        int ipjp = get_Q_index(i + 1, j + 1);
//        int imjm = get_Q_index(i - 1, j - 1);
//        int ipjm = get_Q_index(i + 1, j - 1);
//        int imjp = get_Q_index(i - 1, j + 1);

        double temp_const = const1 * j;
        A[ij][ijp] -= temp_const;
        A[ij][ijm] += temp_const;

        temp_const = const3 * std::pow(j, 2 * beta - 2) * j * j;
        A[ij][ijp] -= temp_const;
        A[ij][ij] += 2 * temp_const;
        A[ij][ijm] -= temp_const;

        temp_const = tau * (kappa / h_v - i) / 2;
        A[ipj][ij] -= temp_const;
        A[imj][ij] += temp_const;

        temp_const = const4 * i;
        A[ipj][ij] -= temp_const;
        A[ij][ij] += 2 * temp_const;
        A[imj][ipj] -= temp_const;

        temp_const = const6 * std::pow(j, beta);
        A[ipj][ipj] -= temp_const;
        A[ipj][ijm] += temp_const;
        A[imj][ijp] += temp_const;
        A[imj][imj] -= temp_const;

        A[ij][ij] += tau * r;
      }
    }

    // Boundary conditions
    // Left Boundary
    for (int i = 0; i <= m; i++) {
      int i0 = get_Q_index(i, 0);
      A[i0][i0] += 1;
    }

    // Bottom Boundary
    for (int j = 0; j <= n; j++) {
      int oj = get_Q_index(0, j);
      A[oj][oj] = 1;
    }

    // Top Boundary
    for (int j = 1; j <= n - 1; j++) {
      int mj = get_Q_index(m, j);
      int mmj = get_Q_index(m - 1, j);
      A[mj][mj] = 1;
      A[mj][mmj] = -1;
    }


    // Right Boundary
    for (int i = 1; i <= m; i++) {
      int i_n = get_Q_index(i, n);
      A[i_n][i_n] = 1;
    }
  }

public:
  double **A; ///< Coefficient matrix. Size = NxN
  size_t N;      ///< Total number of equations (N = (n + 1) * (m + 1))

private:
  int n;          ///< The number of divisions for S_max grid
  int m;          ///< The number of divisions for V_max grid
  double r;       ///< Risk-free rate
  double d;       ///< Dividend yield
  double sigma;   ///< Volatility
  double epsilon; ///< Volatility of volatility
  double kappa;   ///< Mean reversion speed
  double theta;   ///< Long-term mean level
  double beta;    ///< Sensitivity to market movements
  double h_s;     ///< Step size in asset price
  double h_v;     ///< Step size in volatility
  double tau;     ///< Time step size
  double rho;     ///< Correlation with the market
  int p;          ///< Maximum time step index (p)
  int K;          ///< Strike price
  double S_0;     ///< Initial stock price
  double V_0;     ///< Initial volatility

  /**
   * @brief Computes the index in the linear system for given grid coordinates.
   *
   * @param i The index in the volatility dimension.
   * @param j The index in the asset price dimension.
   * @return The corresponding index in the solution vector.
   */
  [[nodiscard]] int get_Q_index(int i, int j) const { return i * (n + 1) + j; }
};
} // namespace ADAAI
