#include "MS.hpp"
#include "parse_data.hpp"
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

// Function to perform parameter tuning
void tune_params( const std::string& file_path )
{
  // Define parameter grids
  std::vector<double> a_values      = { 2, 3, 4 };
  std::vector<double> beta_values   = { 2, 3, 4, 5, 7, 10 };
  std::vector<double> alpha1_values = { -1, 2, 3, 4, 5, 6, 10 };
  std::vector<double> alpha2_values = { -1, 2, 3, 4, 5, 6, 10 };
  std::vector<double> alpha3_values = { -1, 2, 3, 4, 5, 6, 10 };
  std::vector<double> alpha4_values = { -1, 2, 3, 4, 5, 6, 10 };
  std::vector<double> delta_values  = { 0.1, 0.2, 0.3, 0.4, 0.5 };

  Records data;
  parseCSV( file_path, data );
  int              K           = 144;
  float            M           = 100;
  int              window_size = 300;
  int              tau         = 120;
  MomentumStrategy strategy( data.sell_volume, data.buy_volume, data.volume, data.mid_px, K, M, window_size, tau, 0, 0, 0, 0, 0, 0 );
  int              T = data.mid_px.size();

  // Open the CSV file for logging
  std::ofstream log_file( "logs.csv" );
  if ( !log_file.is_open() )
  {
    std::cerr << "Failed to open log file!" << std::endl;
    return;
  }

  // Write the header to the CSV file
  log_file << "a,beta,alpha1,alpha2,alpha3,alpha4,SortinoRatio\n";

  for ( auto delta : delta_values )
  {
    for ( auto beta : beta_values )
    {
      for ( auto a : a_values )
      {
        for ( auto alpha1 : alpha1_values )
        {
          for ( auto alpha2 : alpha2_values )
          {
            for ( auto alpha3 : alpha3_values )
            {
              for ( auto alpha4 : alpha4_values )
              {
                // Set current parameters
                strategy.a_      = a;
                strategy.beta_   = beta;
                strategy.alpha1_ = alpha1;
                strategy.alpha2_ = alpha2;
                strategy.alpha3_ = alpha3;
                strategy.alpha4_ = alpha4;
                strategy.delta_  = delta;

                // Get Sortino ratio
                double sortino = strategy.get_sortino_ratio( T );

                if ( sortino > 0 )
                {
                  // Log the parameters and Sortino ratio to the CSV file
                  log_file << strategy.a_ << ","
                           << strategy.beta_ << ","
                           << strategy.alpha1_ << ","
                           << strategy.alpha2_ << ","
                           << strategy.alpha3_ << ","
                           << strategy.alpha4_ << ","
                           << sortino << "\n";
                }
              }
            }
          }
        }
      }
    }
  }

  // Close the log file
  log_file.close();
}

int main()
{
  tune_params( "../data/train.csv" );
  return 0;
}
