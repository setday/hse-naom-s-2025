#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

class MomentumStrategy
{
public:
  // Parameters
  double a_;
  double beta_;
  double alpha1_;
  double alpha2_;
  double alpha3_;
  double alpha4_;

private:
  std::vector<int>    sell_volume_;
  std::vector<int>    buy_volume_;
  std::vector<double> mid_px_;
  std::vector<int>    volume_;

  std::vector<int> volume_prefix_sum_;
  std::vector<int> buy_volume_prefix_sum_;

  // Hyperparameters
  int    K_;
  double risk_control_;
  int    window_size_;
  int    tau_;
  int    start_tau_;
  double trading_fee_ = 0.0005; // 0.05%
  double delta_ = 0.5;

  std::vector<double> pos_history_;

public:
  MomentumStrategy( std::vector<int>& sell_volume, std::vector<int>& buy_volume, std::vector<int>& volume, std::vector<double>& mid_px, int K = 144, double risk_control = 1, int window_size = 300, int tau = 60, double a = 0, double beta = 0, double alpha1 = 0, double alpha2 = 0, double alpha3 = 0, double alpha4 = 0 )
  {
    // Initialize params and hyperparams
    a_      = a;
    beta_   = beta;
    alpha1_ = alpha1;
    alpha2_ = alpha2;
    alpha3_ = alpha3;
    alpha4_ = alpha4;

    K_            = K;
    risk_control_ = risk_control;
    window_size_  = window_size;
    tau_          = tau;
    start_tau_    = window_size_ * ( K_ + 2 );

    sell_volume_ = sell_volume;
    buy_volume_  = buy_volume;
    volume_      = volume;
    mid_px_      = mid_px;
    // Initialize prefix sums
    volume_prefix_sum_     = { 0 };
    buy_volume_prefix_sum_ = { 0 };
    for ( size_t i = 0; i < volume_.size(); ++i )
    {
      volume_prefix_sum_.push_back( volume_[i] + volume_prefix_sum_.back() );
      buy_volume_prefix_sum_.push_back( buy_volume_[i] + buy_volume_prefix_sum_.back() );
    }

    pos_history_.reserve(volume_.size());
  }

  double get_position( int t )
  {
    if ( t == start_tau_ )
    {
      return 0;
    }
    double x = get_x( t );
    if ( std::abs( x ) < a_ )
    {
      return 0;
    }
    if ( x >= a_ )
    {
      return 2 * risk_control_ / M_PI * std::atan( std::pow( ( x - a_ ), beta_ ) );
    }
    return -2 * risk_control_ / M_PI * std::atan( std::pow( ( -x - a_ ), beta_ ) );
  }

  double get_x( int t )
  {
    return alpha1_ * get_x1( t ) + alpha2_ * get_x2( t ) +
           alpha3_ * get_x3( t ) + alpha4_ * get_x4( t );
  }

  double get_x1( int t )
  {
    return get_mu( t ) / get_mu_std( t );
  }

  double get_x2( int t )
  {
    return ( 3 * get_x1( t ) - 4 * get_x1( t - window_size_ ) + get_x1( t - 2 * window_size_ ) ) /
           ( 2 * window_size_ );
  }

  double get_x3( int t )
  {
    return get_x3_tilda( t ) / get_x3_tilda_std( t );
  }

  double get_x4( int t )
  {
    return ( 3 * get_x3( t ) - 4 * get_x3( t - window_size_ ) + get_x3( t - 2 * window_size_ ) ) /
           ( 2 * window_size_ * get_x3_tilda_std( t ) );
  }

  double get_sigma_squared( int t )
  {
    double sigma_squared = 0;
    for ( int i = t; i > t - window_size_; --i )
    {
      sigma_squared += std::log( mid_px_[i] / mid_px_[i - 1] ) * std::log( mid_px_[i] / mid_px_[i - 1] );
    }
    return sigma_squared / window_size_;
  }

  double get_mu( int t )
  {
    double sigma_squared = get_sigma_squared( t );
    return sigma_squared / 2 + std::log( mid_px_[t] / mid_px_[t - window_size_] ) / window_size_;
  }

  double get_mu_std( int t )
  {
    return std::sqrt( get_sigma_squared( t ) / window_size_ );
  }

  double get_x3_tilda( int t )
  {
    int total_volume = volume_prefix_sum_[t] - volume_prefix_sum_[t - window_size_];
    if ( total_volume == 0 )
    {
      return 0;
    }
    int buy_volume = buy_volume_prefix_sum_[t] - buy_volume_prefix_sum_[t - window_size_];
    return 2.0 * buy_volume / total_volume - 1;
  }

  double get_x3_tilda_mean( int t )
  {
    double x3_tilda_mean = 0;
    for ( int k = 0; k < K_; ++k )
    {
      x3_tilda_mean += get_x3_tilda( t - k * window_size_ );
    }
    return x3_tilda_mean / K_;
  }

  double get_x3_tilda_std( int t )
  {
    double x3_tilda_mean = get_x3_tilda_mean( t );
    double x3_tilda_std  = 0;
    for ( int k = 0; k < K_; ++k )
    {
      x3_tilda_std += std::pow( get_x3_tilda( t - k * window_size_ ) - x3_tilda_mean, 2 );
    }
    return x3_tilda_std / ( K_ - 1 );
  }

  void compute_next_position( int t ) {
    assert(!pos_history_.empty());

    double pos = get_position( t );
    double prev_pos = pos_history_.back();

//    std::cout << std::abs(pos - prev_pos) << std::endl;

    if (std::abs(pos - prev_pos) <= delta_) {
      pos_history_.push_back(prev_pos);
      return;
    }

    pos_history_.push_back(pos);
  }

  double get_delta_PnL( int t1, int t2 )
  {
    int t1_pos = (t1 - start_tau_) / tau_;
    int t2_pos = (t2 - start_tau_) / tau_;

    assert( pos_history_.size() >= std::max( t1_pos , t2_pos ) );

    double pos_t1 = pos_history_[ t1_pos ];
    double pos_t2 = pos_history_[ t2_pos ];
    double delta_pos = pos_t2 - pos_t1;

    return pos_t1 * ( mid_px_[t2] - mid_px_[t1] ) - std::abs(delta_pos) * trading_fee_ * mid_px_[t2];
  }

  double get_PnL( int T )
  {
    double total_PnL = 0;
    for ( int i = start_tau_; i < T - tau_; i += tau_ )
    {
      total_PnL += get_delta_PnL( i, i + tau_ );
//      if ( total_PnL <= -risk_control_ ) {
//        std::cout << "Bankruptcy of a bot" << std::endl;
//        break;
//      }
    }
    return total_PnL;
  }

public:
  double get_sortino_ratio( int T )
  {
    double total_PnL    = 0;
    double max_drawdown = 0;

    pos_history_.clear();
    pos_history_.push_back(0);

    for ( int i = start_tau_; i < T - tau_; i += tau_ )
    {
      double delta_pnl;

      compute_next_position( i + tau_ );

      delta_pnl = get_delta_PnL( i, i + tau_ );
      total_PnL += delta_pnl;
      if ( total_PnL < 0 )
      {
        max_drawdown = std::max( max_drawdown, -total_PnL );
      }
//      if ( total_PnL <= -risk_control_ ) {
//        break;
//      }
    }
    if ( max_drawdown == 0 )
    {
      std::cout << "WARNING: Zero Drawdown!" << '\n';
    }
    return total_PnL / max_drawdown;
  }
};
