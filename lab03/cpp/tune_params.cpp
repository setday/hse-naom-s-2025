#include <iostream>
#include <unordered_map>
#include <vector>

#include "MS.hpp"
#include "parse_data.hpp"

// Function to perform parameter tuning
void tune_params( const std::string& file_path )
{
  // Define parameter grids
  std::vector<double> a_values      = { 2, 2.5, 3 };
  std::vector<double> beta_values   = { 1, 2, 3, 4, 5 };
  std::vector<double> alpha1_values = { 1, 3, 5 };
  std::vector<double> alpha2_values = { 1, 3, 5 };
  std::vector<double> alpha3_values = { 1, 3, 5 };
  std::vector<double> alpha4_values = { 1, 3, 5 };


  // Iterate through all parameter combinations
  double                                  best_sortino = -std::numeric_limits<double>::infinity();
  std::unordered_map<std::string, double> best_params;

  Records data;
  parseCSV( file_path, data );
  int              K           = 144;
  float            M           = 1;
  int              window_size = 300;
  int              tau         = 120;
  MomentumStrategy strategy( data.sell_volume, data.buy_volume, data.volume, data.mid_px, K, M, window_size, tau, 0, 0, 0, 0, 0, 0 );
  int              T = data.mid_px.size();

  for ( size_t j = 0; j < beta_values.size(); ++j )
  {
    for ( size_t i = 0; i < a_values.size(); ++i )
    {
      for ( size_t k = 0; k < alpha1_values.size(); ++k )
      {
        for ( size_t l = 0; l < alpha2_values.size(); ++l )
        {
          for ( size_t m = 0; m < alpha3_values.size(); ++m )
          {
            for ( size_t n = 0; n < alpha4_values.size(); ++n )
            {
              // Create strategy with current parameter set
              strategy.a_      = a_values[i];
              strategy.beta_   = beta_values[j];
              strategy.alpha1_ = alpha1_values[k];
              strategy.alpha2_ = alpha2_values[l];
              strategy.alpha3_ = alpha3_values[m];
              strategy.alpha4_ = alpha4_values[n];

              // Get Sortino ratio
              double sortino = strategy.get_sortino_ratio( T );
              std::cout << "Sortino Ratio: " << sortino << '\n';

              // Update best parameters if needed
              if ( sortino > best_sortino )
              {
                best_sortino = sortino;
                best_params  = { { "a", a_values[i] }, { "beta", beta_values[j] }, { "alpha1", alpha1_values[k] }, { "alpha2", alpha2_values[l] }, { "alpha3", alpha3_values[m] }, { "alpha4", alpha4_values[n] } };
              }
            }
          }
        }
      }
    }
  }

  // Output the best parameters
  std::cout << "Best Parameters: ";
  for ( const auto& param : best_params )
  {
    std::cout << param.first << "=" << param.second << " ";
  }
  std::cout << "\nBest Sortino Ratio: " << best_sortino << std::endl;
}

int main()
{
  tune_params( "../data/train.csv" );
  return 0;
}
