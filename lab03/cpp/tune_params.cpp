#include <iostream>
#include <unordered_map>
#include <vector>

#include "MS.hpp"
#include "parse_data.hpp"

// Function to perform parameter tuning
void tune_params( const std::string& file_path )
{
  // Define parameter grids
  std::vector<double> a_values      = { 4.5, 4, 5 };
  std::vector<double> beta_values   = { 3 };
  std::vector<double> alpha1_values = { 4.5, 4, 5 };
  std::vector<double> alpha2_values = { -5, -4, -6 };
  std::vector<double> alpha3_values = { -1, 0, -2 };
  std::vector<double> alpha4_values =  { -10, -9, -8 };

//  Best Parameters: alpha4=3 alpha3=-1 alpha2=-1 alpha1=5 beta=5 a=3

  // Iterate through all parameter combinations
  double                                  best_sortino = -std::numeric_limits<double>::infinity();
  double                                  best_PnL = -std::numeric_limits<double>::infinity();
  std::unordered_map<std::string, double> best_params;

  Records data;
  parseCSV( file_path, data );
  int              K           = 150;
  float            M           = 1;
  int              window_size = 500;
  int              tau         = 60;
  MomentumStrategy strategy( data.sell_volume, data.buy_volume, data.volume, data.mid_px, K, M, window_size, tau, 0, 0, 0, 0, 0, 0 );
  int              T = data.mid_px.size();

  for (double & beta_value : beta_values)
  {
    for (double & a_value : a_values)
    {
      for (double & alpha1_value : alpha1_values)
      {
        for (double & alpha2_value : alpha2_values)
        {
          for (double & alpha3_value : alpha3_values)
          {
            for (double & alpha4_value : alpha4_values)
            {
              // Create strategy with current parameter set
              strategy.a_      = a_value;
              strategy.beta_   = beta_value;
              strategy.alpha1_ = alpha1_value;
              strategy.alpha2_ = alpha2_value;
              strategy.alpha3_ = alpha3_value;
              strategy.alpha4_ = alpha4_value;

              // Get Sortino ratio
              double sortino = strategy.get_sortino_ratio( T );

              double PnL = strategy.get_PnL( T );
              std::cout << "Sortino Ratio: " << sortino << '\n';
              std::cout << "PnL: " << PnL << " RUB\n";

              // Update the best parameters if needed
              if ( sortino > best_sortino )
              {
                best_sortino = sortino;
                best_params  = { { "a", a_value }, { "beta", beta_value }, { "alpha1", alpha1_value }, { "alpha2", alpha2_value }, { "alpha3", alpha3_value }, { "alpha4", alpha4_value } };
              }

              if ( PnL > best_PnL )
              {
                best_PnL = PnL;
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
  tune_params( "lab03/data/train.csv" );
  return 0;
}
