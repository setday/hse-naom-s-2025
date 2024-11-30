#include <iostream>
#include <unordered_map>
#include <vector>

#include "MS.hpp"

// Function to perform parameter tuning
void tune_params( const std::string& file_path )
{
  // Define parameter grids
  std::vector<double> a_values      = { 2, 2.5, 3 };
  std::vector<double> alpha1_values = { 1, 3, 5 };
  std::vector<double> alpha2_values = { 1, 3, 5 };
  std::vector<double> alpha3_values = { 1, 3, 5 };
  std::vector<double> alpha4_values = { 1, 3, 5 };
  std::vector<double> beta_values   = { 1, 3, 5 };

  // Iterate through all parameter combinations
  double                                  best_sortino = -std::numeric_limits<double>::infinity();
  std::unordered_map<std::string, double> best_params;

  for ( double a : a_values )
  {
    for ( double beta : beta_values )
    {
      for ( double alpha1 : alpha1_values )
      {
        for ( double alpha2 : alpha2_values )
        {
          for ( double alpha3 : alpha3_values )
          {
            for ( double alpha4 : alpha4_values )
            {
              // Create strategy with current parameter set
              MomentumStrategy strategy( 144, 1, 150, file_path, a, beta, alpha1, alpha2, alpha3, alpha4 );
              int              T = strategy.data_.size();

              // Get Sortino ratio
              double sortino = strategy.get_sortino_ratio( T );
              std::cout << "Sortino Ratio: " << sortino << '\n';

              // Update best parameters if needed
              if ( sortino > best_sortino )
              {
                best_sortino = sortino;
                best_params  = { { "a", a }, { "beta", beta }, { "alpha1", alpha1 }, { "alpha2", alpha2 }, { "alpha3", alpha3 }, { "alpha4", alpha4 } };
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