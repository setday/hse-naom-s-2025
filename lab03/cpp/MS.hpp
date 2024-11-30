#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

class MomentumStrategy
{
public:
  struct DataPoint
  {
    int    sell_volume;
    int    buy_volume;
    double mid_px;
  };

  std::vector<DataPoint>   data_;
  std::vector<int>         sell_volume_;
  std::vector<int>         buy_volume_;
  std::vector<int>         volume_;
  std::vector<int>         volume_prefix_sum_;
  std::vector<int>         buy_volume_prefix_sum_;
  std::vector<double>      mid_px_;

  // Parameters and hyperparameters
  std::unordered_map<std::string, double> params_;
  std::unordered_map<std::string, double> hyperparams_;

public:
  MomentumStrategy( int K = 144, double risk_control = 1, int window_size = 300, const std::string& file_path = "../train.csv", double a = 0, double beta = 0, double alpha1 = 0, double alpha2 = 0, double alpha3 = 0, double alpha4 = 0 )
  {
    params_      = { { "a", a }, { "beta", beta }, { "alpha1", alpha1 }, { "alpha2", alpha2 }, { "alpha3", alpha3 }, { "alpha4", alpha4 } };
    hyperparams_ = { { "K", K }, { "risk_control", risk_control }, { "window_size", window_size }, { "tau", 60 } };
    // Load data from CSV
    load_data( file_path );
    // Initialize prefix sums
    volume_prefix_sum_     = { 0 };
    buy_volume_prefix_sum_ = { 0 };
    for ( size_t i = 0; i < volume_.size(); ++i )
    {
      volume_prefix_sum_.push_back( volume_[i] + volume_prefix_sum_.back() );
      buy_volume_prefix_sum_.push_back( buy_volume_[i] + buy_volume_prefix_sum_.back() );
    }
  }

  // CSV Parser Function
  std::vector<std::string> getNextLineAndSplitIntoTokens( std::istream& str )
  {
    std::vector<std::string> result;
    std::string              line;
    std::getline( str, line );

    std::stringstream lineStream( line );
    std::string       cell;

    while ( std::getline( lineStream, cell, ',' ) )
    {
      result.push_back( cell );
    }
    // This checks for a trailing comma with no data after it.
    if ( !lineStream && cell.empty() )
    {
      result.push_back( "" ); // If there was a trailing comma then add an empty element.
    }
    return result;
  }

  // Function to load data from CSV file
  void load_data( const std::string& file_path )
  {
    std::ifstream file( file_path );
    std::string   line;
    bool          first_line = true;
    while ( std::getline( file, line ) )
    {
      std::stringstream        lineStream( line );
      std::vector<std::string> tokens = getNextLineAndSplitIntoTokens( lineStream );

      if ( first_line )
      {
        first_line = false; // Skip the header
        continue;
      }

      DataPoint dp;
      dp.sell_volume = std::stoi( tokens[0] );
      dp.buy_volume  = std::stoi( tokens[1] );
      dp.mid_px      = std::stod( tokens[2] );

      data_.push_back( dp );
      sell_volume_.push_back( dp.sell_volume );
      buy_volume_.push_back( dp.buy_volume );
      volume_.push_back( dp.sell_volume + dp.buy_volume );
      mid_px_.push_back( dp.mid_px );
    }
  }

  double get_position( int t )
  {
    if ( t == 0 )
    {
      return 0;
    }
    double x = get_x( t );
    if ( std::abs( x ) < params_["a"] )
    {
      return 0;
    }
    if ( x >= params_["a"] )
    {
      return 2 * hyperparams_["risk_control"] / M_PI * std::atan( std::pow( ( x - params_["a"] ), params_["beta"] ) );
    }
    return -2 * hyperparams_["risk_control"] / M_PI * std::atan( std::pow( ( -x - params_["a"] ), params_["beta"] ) );
  }

  double get_x( int t )
  {
    return params_["alpha1"] * get_x1( t ) + params_["alpha2"] * get_x2( t ) +
           params_["alpha3"] * get_x3( t ) + params_["alpha4"] * get_x4( t );
  }

  double get_x1( int t )
  {
    return get_mu( t ) / get_mu_std( t );
  }

  double get_x2( int t )
  {
    return ( 3 * get_x1( t ) - 4 * get_x1( t - hyperparams_["window_size"] ) + get_x1( t - 2 * hyperparams_["window_size"] ) ) /
           ( 2 * hyperparams_["window_size"] );
  }

  double get_x3( int t )
  {
    return get_x3_tilda( t ) / get_x3_tilda_std( t );
  }

  double get_x4( int t )
  {
    return ( 3 * get_x3( t ) - 4 * get_x3( t - hyperparams_["window_size"] ) + get_x3( t - 2 * hyperparams_["window_size"] ) ) /
           ( 2 * hyperparams_["window_size"] * get_x3_tilda_std( t ) );
  }

  double get_sigma_squared( int t )
  {
    double sigma_squared = 0;
    for ( int i = t; i > t - hyperparams_["window_size"]; --i )
    {
      sigma_squared += std::log( mid_px_[i] / mid_px_[i - 1] ) * std::log( mid_px_[i] / mid_px_[i - 1] );
    }
    return sigma_squared / hyperparams_["window_size"];
  }

  double get_mu( int t )
  {
    double sigma_squared = get_sigma_squared( t );
    return sigma_squared / 2 + std::log( mid_px_[t] / mid_px_[t - hyperparams_["window_size"]] ) / hyperparams_["window_size"];
  }

  double get_mu_std( int t )
  {
    return std::sqrt( get_sigma_squared( t ) / hyperparams_["window_size"] );
  }

  double get_x3_tilda( int t )
  {
    int total_volume = volume_prefix_sum_[t] - volume_prefix_sum_[t - hyperparams_["window_size"]];
    if ( total_volume == 0 )
    {
      return 0;
    }
    int buy_volume = buy_volume_prefix_sum_[t] - buy_volume_prefix_sum_[t - hyperparams_["window_size"]];
    return 2.0 * buy_volume / total_volume - 1;
  }

  double get_x3_tilda_mean( int t )
  {
    double x3_tilda_mean = 0;
    for ( int k = 0; k < hyperparams_["K"]; ++k )
    {
      x3_tilda_mean += get_x3_tilda( t - k * hyperparams_["window_size"] );
    }
    return x3_tilda_mean / hyperparams_["K"];
  }

  double get_x3_tilda_std( int t )
  {
    double x3_tilda_mean = get_x3_tilda_mean( t );
    double x3_tilda_std  = 0;
    for ( int k = 0; k < hyperparams_["K"]; ++k )
    {
      x3_tilda_std += std::pow( get_x3_tilda( t - k * hyperparams_["window_size"] ) - x3_tilda_mean, 2 );
    }
    return x3_tilda_std / ( hyperparams_["K"] - 1 );
  }

  double get_delta_PnL( int t1, int t2 )
  {
    double pos = get_position( t1 );
    return pos * ( mid_px_[t2] - mid_px_[t1] );
  }

  double get_PnL( int T )
  {
    double total_PnL = 0;
    for ( int i = hyperparams_["window_size"] * ( hyperparams_["K"] + 1 ); i < T - hyperparams_["tau"]; i += hyperparams_["tau"] )
    {
      total_PnL += get_delta_PnL( i, i + hyperparams_["tau"] );
    }
    return total_PnL;
  }

  double get_sortino_ratio( int T )
  {
    double total_PnL    = 0;
    double max_drawdown = 0;
    for ( int i = hyperparams_["window_size"] * ( hyperparams_["K"] + 1 ); i < T - hyperparams_["tau"]; i += hyperparams_["tau"] )
    {
      double delta_pnl = get_delta_PnL( i, i + hyperparams_["tau"] );
      total_PnL += delta_pnl;
      if ( total_PnL < 0 )
      {
        max_drawdown = std::max( max_drawdown, -total_PnL );
      }
    }
    if (max_drawdown == 0) {
      std::cout << "WARNING: Zero Drawdown!" << '\n';
    }
    return total_PnL / max_drawdown;
  }
};
