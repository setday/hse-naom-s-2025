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
  std::vector<double> beta_values   = { 2, 3, 10 };
  std::vector<double> alpha1_values = { -1, 3, 5 };
  std::vector<double> alpha2_values = { -1, 3, 5 };
  std::vector<double> alpha3_values = { -1, 3, 5 };
  std::vector<double> alpha4_values = { -1, 3, 5 };

  std::vector<double> alpha5_values = { -1, 3, 5 };
  std::vector<double> alpha6_values = { -1, 3, 5 };
  std::vector<double> alpha7_values = { -1, 3, 5 };
  std::vector<double> alpha8_values = { -1, 3, 5 };

  std::vector<double> delta_values = { 0.3, 0.4, 0.5 };

  Records data;
  parseCSV( file_path, data );
  int              K           = 144;
  float            M           = 100;
  int              window_size = 300;
  int              tau         = 120;
  MomentumStrategy strategy( data.sell_volume, data.buy_volume, data.volume, data.mid_px, K, M, window_size, tau, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );
  int              T = data.mid_px.size();

  // Open the CSV file for logging
  std::ofstream log_file( "logs.csv" );
  if ( !log_file.is_open() )
  {
    std::cerr << "Failed to open log file!" << std::endl;
    return;
  }

  // Write the header to the CSV file
  log_file << "a,beta,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8,SortinoRatio\n";

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
                for ( auto alpha5 : alpha5_values )
                {
                  for ( auto alpha6 : alpha6_values )
                  {
                    for ( auto alpha7 : alpha7_values )
                    {
                      for ( auto alpha8 : alpha8_values )
                      {
                        // Set current parameters
                        strategy.a_      = a;
                        strategy.beta_   = beta;
                        strategy.alpha1_ = alpha1;
                        strategy.alpha2_ = alpha2;
                        strategy.alpha3_ = alpha3;
                        strategy.alpha4_ = alpha4;

                        strategy.alpha5_ = alpha5;
                        strategy.alpha6_ = alpha6;
                        strategy.alpha7_ = alpha7;
                        strategy.alpha8_ = alpha8;

                        strategy.delta_ = delta;

                        // Get Sortino ratio
                        double sortino = strategy.get_sortino_ratio( T );
                        std::cout << sortino << '\n';
                        if ( sortino > 0 )
                        {
                          // Log the parameters and Sortino ratio to the CSV file
                          log_file << strategy.a_ << ","
                                   << strategy.beta_ << ","
                                   << strategy.alpha1_ << ","
                                   << strategy.alpha2_ << ","
                                   << strategy.alpha3_ << ","
                                   << strategy.alpha4_ << ","

                                   << strategy.alpha5_ << ","
                                   << strategy.alpha6_ << ","
                                   << strategy.alpha7_ << ","
                                   << strategy.alpha8_ << ","

                                   << sortino << "\n";
                        }
                      }
                    }
                  }
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
