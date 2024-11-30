#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Records
{
  std::vector<int>    sell_volume;
  std::vector<int>    buy_volume;
  std::vector<int>    volume;
  std::vector<double> mid_px;
};


void parseCSV( const std::string& filename, Records& records )
{
  std::ifstream file( filename );
  std::string   line;

  // Skip the header
  std::getline( file, line );

  while ( std::getline( file, line ) )
  {
    std::stringstream ss( line );
    std::string       temp;

    // Parse first column (sell_volume)
    std::getline( ss, temp, ',' );
    records.sell_volume.push_back( std::stoi( temp ) );

    // Parse second column (buy_volume)
    std::getline( ss, temp, ',' );
    records.buy_volume.push_back( std::stoi( temp ) );

    // Add total volume
    records.volume.push_back( records.sell_volume.back() + records.buy_volume.back() );

    // Parse third column (mid_px)
    std::getline( ss, temp, ',' );
    records.mid_px.push_back( std::stof( temp ) );
  }
  file.close();
}
