#include <cmath>

namespace ADAII
{
    class SLAE 
    {
    /*
    This class defines a system of linear algebraic equations (Ax = b) that will be solved using Gaussian elimination. 
    The system arises from the numerical solution of a SDE (stochastic differential equation) using an implicit method.
    */

    public:
        SLAE (int n, int m, double r, double d, double sigma, double epsilon, double kappa, double theta, double beta, double h_s, double h_v, double tau, double rho, int p, int K) : 
            n(n), m(m), d(d), sigma(sigma), epsilon(epsilon), kappa(kappa), theta(theta), beta(beta), h_s(h_s), h_v(h_v), tau(tau), rho(rho), p(p), K(K) 
        {
            N = (n + 1) * (m + 1);

            x = new double[N]();
            b = new double[N]();
            A = new double*[N];
            for (int i = 0; i < N; ++i) {
                A[i] = new double[N]();
            }
        }

        ~SLAE() {
            delete[] x;
            delete[] b;
            for (int i = 0; i < N; ++i) {
                delete[] A[i];
            }
            delete[] A;
        }

        int get_Q_index (int i, int j) const {
            return i * (n + 1) + j;
        }

        void construct_SLAE(int k) {
            // Inner matrix construction

            // We define many constants which don't have any sense, but
            // they are used many times in different places. It is only a matter of optimization.
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
                    int ipjp = get_Q_index(i + 1, j + 1);
                    int imjm = get_Q_index(i - 1, j - 1);
                    int ipjm = get_Q_index(i + 1, j - 1);
                    int imjp = get_Q_index(i - 1, j + 1);

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

                    // TODO
                    // how to get previous xs? Do wee need to pass them as an argument?
                    b[ij] += prev_x[ij]; // RHS
                }
            }

            // Boundary conditions
            // Left Boundary
            for (int i = 0; i <= m; i++) {
                int i0 = get_Q_index(i, 0);
                A[i0][i0] += 1;
                b[i0] = 0;
            }

            // Lower Boundary
            double const7 = (p - k + 1) * tau;
            double const8 = K * std::exp((d - r) * const7) / h_s;
            double const9 = std::exp(-d * const7) * h_s;
            double const10 = std::exp(-r * const7) * K;
            for (int j = 0; j <= n; j++) {
                int oj = get_Q_index(0, j); 
                A[oj][oj] = 1;
                if (j >= const8) {
                    b[oj] = const9 * j - const10;
                } else {
                    b[oj] = 0;
                }
            }

            // Right Boundary
            double const11 = const9 * n - const10;
            for (int i = 1; i <= m; i++) {
                int i_n = get_Q_index(i, n);
                A[i_n][i_n] = 1;
                b[i_n] = const11;
            }

            // Upper Boundary
            for (int j = 1; j <= n - 1; j++) {
                int mj = get_Q_index(m, j);
                int mmj = get_Q_index(m - 1, j);
                A[mj][mj] = 1;
                A[mj][mmj] = -1;
                b[mj] = 0;
            }
        }

    private:
        int n; // rows
        int m; // columns
        int N; // Computed as (n + 1) * (m + 1)
        
        double r; // interest rate
        double d; // dividend rate
        double sigma; // volatility
        double epsilon;
        double kappa;
        double theta;
        double beta;
        double h_s; // S_max / n
        double h_v; // V_max / m
        double tau;
        double rho;
        int p; // MAX value for k (k = p...1)
        int K; // strike price

        double* x; // Solution vector
        double* b; // RHS vector
        double** A; // Coefficient matrix
    };
}
