#pragma once

#include <cmath>
#include <cstring>
#include <iostream>

#include "../../HSE_NaOM_S2024/integration/intergartor/RHS.hpp"

#include "auxiliary_functions.h"
#include "gaussian_elimination.h"
#include "option_environment.h"

namespace ADAAI::LAB01 {
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
  explicit SLAE(OptionEnvironment *option_environment)
      : m_oe(option_environment), n(OptionEnvironment::n),
        m(OptionEnvironment::m) {
    A = new double *[N];
    for (int i = 0; i < N; ++i) {
      A[i] = new double[N];
      std::memset(A[i], 0, N * sizeof(double));
    }

    construct_SLAE();
  }

  /**
   * @brief Destructor for the SLAE object.
   *
   * Frees the memory allocated for the coefficient matrix.
   */
  ~SLAE() {
    for (int i = 0; i < N; ++i) {
      delete[] A[i];
    }
    delete[] A;
  }

  /**
   * @brief Sets the initial conditions for the right-hand side vector (b).
   * @param rhs The right-hand side vector to be updated.
   */
  void set_initial_RHS(double *rhs) { // utilize initial condition @t=T
    for (int i = 1; i <= m - 1; i++) {
      for (int j = 1; j <= n - 1; j++) {
        int ij = get_Q_index(i, j);
        rhs[ij] = m_oe->get_payoff_for_S_j(j);
      }
    }

    update_RHS_boundary(rhs, 0.0); // utilize boundary condition @tau=0 (option is expired)
  }

  /**
   * @brief Updates the right-hand side vector (b) with boundary conditions.
   * @param rhs The right-hand side vector to be updated.
   * @param current_tau Time to option expiration.
   */
  void update_RHS_boundary(double *rhs, double current_tau) const {
    double K_star = m_oe->K * std::exp(-(m_oe->r - m_oe->d) * current_tau);

    // Lower-Upper Boundary
    // double const8 = m_oe->K * std::exp((m_oe->d - m_oe->r) * current_tau) / m_oe->h_s_max;
    double const9 = std::exp(-m_oe->d * current_tau) * m_oe->h_s_max;
    double const10 = std::exp(-m_oe->r * current_tau) * K_star;

    for (int j = 0; j <= n; j++) {
      int oj = get_Q_index(0, j);
      rhs[oj] = std::max(const9 * j - const10, 0.0);  // relu(e^(-d * tau) * S - e^(-r * tau) * K)

      int mj = get_Q_index(m, j);
      rhs[mj] = const9 * j; // e^(-d * tau) * S
    }

    // Left-Right Boundary
    double const11 = const9 * n - const10;

    for (int i = 1; i <= m - 1; i++) {
      int in = get_Q_index(i, n);
      rhs[in] = const11; // e^(-d * tau) * S_max - e^(-r * tau) * K

      int i0 = get_Q_index(i, 0);
      rhs[i0] = 0;
    }
  }

  /**
   * @brief Retrieves the solution for the given initial stock price and
   * volatility.
   *
   * Computes the indices corresponding to the current stock price (S_0)
   * and volatility (V_0) in the solution vector and returns the corresponding
   * value from the solution.
   *
   * @return The solution corresponding to the initial conditions.
   */
  double get_the_answer(double const *rhs) {
    int i = m_oe->get_closest_S_index(m_oe->S_0);
    int j = m_oe->get_closest_V_index(m_oe->V_0);

    return rhs[get_Q_index(i, j)];
  }

  /**
   * @brief Retrieves the implicit matrix.
   *
   * Copies the internal matrix to the specified matrix.
   *
   * @param matrix The matrix to be populated with the implicit matrix.
   */
  void get_implicit_matrix(double **matrix) const {
    for (int i = 0; i < N; i++) {
      memcpy(matrix[i], A[i], N * sizeof(double));
    }
  }

  void operator()(double current_time, const double *current_state,
                  double *rhs) const override {
    double const_ds = (m_oe->r - m_oe->d) / 2;
    double const_ds_ds = m_oe->sigma * m_oe->sigma / 2 * m_oe->h_v_max *
                         std::pow(m_oe->h_s_max, 2 * m_oe->beta - 2);
    double const_dv_dv = m_oe->epsilon * m_oe->epsilon / (2 * m_oe->h_v_max);
    double const_dv_ds = m_oe->sigma * m_oe->epsilon * m_oe->rho *
                         std::pow(m_oe->h_s_max, m_oe->beta - 1) / 4;

    for (int i = 1; i <= m - 1; i++) {
      double const_ds_ds_i = const_ds_ds * i;
      double const_dv_dv_i = const_dv_dv * i;
      double const_dv_ds_i = const_dv_ds * i;

      for (int j = 1; j <= n - 1; j++) {
        int ij = get_Q_index(i, j);
        // The second term in LHS (operator L: convective + diffusive terms)
        int ipj = get_Q_index(i + 1, j);
        int ijp = get_Q_index(i, j + 1);
        int imj = get_Q_index(i - 1, j);
        int ijm = get_Q_index(i, j - 1);
        int ipjp = get_Q_index(i + 1, j + 1);
        int imjm = get_Q_index(i - 1, j - 1);
        int ipjm = get_Q_index(i + 1, j - 1);
        int imjp = get_Q_index(i - 1, j + 1);

        rhs[ij] = m_oe->r * current_state[ij];

        double temp_const = const_ds * j;
        rhs[ij] -= temp_const * (current_state[ijp] - current_state[ijm]);

        temp_const = const_ds_ds_i * std::pow(j, 2 * m_oe->beta);
        rhs[ij] -= temp_const * (current_state[ijp] - 2 * current_state[ij] + current_state[ijm]);

        temp_const = m_oe->tau * (m_oe->kappa / m_oe->h_v_max - i) / 2;
        rhs[ij] -= temp_const * (current_state[ipj] - current_state[imj]);

        temp_const = const_dv_dv_i;
        rhs[ij] -= temp_const * (current_state[ipj] - 2 * current_state[ij] + current_state[imj]);

        temp_const = const_dv_ds_i * std::pow(j, m_oe->beta);
        rhs[ij] -= temp_const * (current_state[ipjp] - current_state[ipjm] - current_state[imjp] + current_state[imjm]);
      }
    }

    update_RHS_boundary(rhs, current_time);


    for (int j = 0; j <= n; j++) {
      int oj = get_Q_index(0, j);
      rhs[oj] = (rhs[oj] - current_state[oj]) / m_oe->tau;

      int mj = get_Q_index(m, j);
      rhs[mj] = (rhs[mj] - current_state[mj]) / m_oe->tau;
    }

    for (int i = 1; i <= m - 1; i++) {
      int in = get_Q_index(i, n);
      rhs[in] = (rhs[in] - current_state[in]) / m_oe->tau;

      int i0 = get_Q_index(i, 0);
      rhs[i0] = (rhs[i0] - current_state[i0]) / m_oe->tau;
    }
  }

private:
  /**
   * @brief Constructs the coefficient matrix (A).
   *
   * Populates the internal matrix based on the specified
   * parameters and the boundary conditions derived from the problem.
   */
  void construct_SLAE() {
    // Inner matrix construction

    // We define many constants which don't have any sense, but
    // they are used many times in different places. It is only a matter of
    // optimization.
    double const_ds = m_oe->tau * (m_oe->r - m_oe->d) / 2;
    double const_ds_ds = m_oe->tau * m_oe->sigma * m_oe->sigma / 2 * m_oe->h_v_max *
                    std::pow(m_oe->h_s_max, 2 * m_oe->beta - 2);
    double const_dv_dv =
        m_oe->tau * m_oe->epsilon * m_oe->epsilon / (2 * m_oe->h_v_max);
    double const_dv_ds = m_oe->tau * m_oe->sigma * m_oe->epsilon * m_oe->rho *
                    std::pow(m_oe->h_s_max, m_oe->beta - 1) / 4;

    for (int i = 1; i <= m - 1; i++) {

      double const_ds_ds_i = const_ds_ds * i;
      double const_dv_dv_i = const_dv_dv * i;
      double const_dv_ds_i = const_dv_ds * i;

      for (int j = 1; j <= n - 1; j++) {
        int ij = get_Q_index(i, j);
        // The first term in LHS
        A[ij][ij] += 1;

        // The second term in LHS (operator L: convective + diffusive terms)
        int ipj = get_Q_index(i + 1, j);
        int ijp = get_Q_index(i, j + 1);
        int imj = get_Q_index(i - 1, j);
        int ijm = get_Q_index(i, j - 1);
        // int ipjp = get_Q_index(i + 1, j + 1);
        // int imjm = get_Q_index(i - 1, j - 1);
        // int ipjm = get_Q_index(i + 1, j - 1);
        // int imjp = get_Q_index(i - 1, j + 1);

        double temp_const = const_ds * j;
        A[ij][ijp] -= temp_const;
        A[ij][ijm] += temp_const;

        temp_const = const_ds_ds_i * std::pow(j, 2 * m_oe->beta);
        A[ij][ijp] -= temp_const;
        A[ij][ij] += 2 * temp_const;
        A[ij][ijm] -= temp_const;

        temp_const = m_oe->tau * (m_oe->kappa / m_oe->h_v_max - i) / 2;
        A[ipj][ij] -= temp_const;
        A[imj][ij] += temp_const;

        temp_const = const_dv_dv_i;
        A[ipj][ij] -= temp_const;
        A[ij][ij] += 2 * temp_const;
        A[imj][ipj] -= temp_const;

        temp_const = const_dv_ds_i * std::pow(j, m_oe->beta);
        A[ipj][ipj] -= temp_const;
        A[ipj][ijm] += temp_const;
        A[imj][ijp] += temp_const;
        A[imj][imj] -= temp_const;

        A[ij][ij] += m_oe->tau * m_oe->r;
      }
    }

    // Boundary conditions
    // Left Boundary
    for (int i = 0; i <= m; i++) {
      int i0 = get_Q_index(i, 0);
      A[i0][i0] = 0 ;
    }

    // Bottom Boundary
    for (int j = 0; j <= n; j++) {
      int oj = get_Q_index(0, j);
      A[oj][oj] = 1;
    }

    // Top Boundary
    for (int j = 1; j <= n - 1; j++) {
      int mj = get_Q_index(m, j);
      A[mj][mj] = 1;
      // int mmj = get_Q_index(m - 1, j);
      //  A[mj][mmj] = -1;
    }

    // Right Boundary
    for (int i = 1; i <= m; i++) {
      int in = get_Q_index(i, n);
      A[in][in] = 1;
    }
  }

public:
  constexpr static int N =
      (OptionEnvironment::n + 1) * (OptionEnvironment::m + 1); ///< Total number of equations (N = (n + 1) * (m + 1))
  double **A;                                                  ///< Coefficient matrix. Size = NxN

private:
  OptionEnvironment *m_oe; ///< Pointer to the OptionEnvironment object
  int n;                   ///< Number of divisions for the asset price grid
  int m;                   ///< Number of divisions for the volatility grid

  /**
   * @brief Computes the index in the linear system for given grid coordinates.
   *
   * @param v The index in the volatility dimension.
   * @param s The index in the asset price dimension.
   * @return The corresponding index in the solution vector.
   */
  [[nodiscard]] int get_Q_index(int v, int s) const { return v * (n + 1) + s; }
};
} // namespace ADAAI::LAB01
