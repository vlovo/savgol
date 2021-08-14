
 
#include "savgol.hpp"

#include <vector>


int main(int argc, char **argv[])
{
 
using namespace  filter;

  std::vector<float> data{1.0, 2, 1, 2, 1.0, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2};
  
  std::vector<float> r;
   

  savgol(data.begin(), data.end(), std::back_inserter(r), SmoothQuadCubic(5));


  std::vector<double> data2{1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2};
  
  
  std::vector<double> r2(data2);
  savgol(data.begin(), data.end(), r2.begin(), DeriveQuadFirst(7));

  
  return 0;
}
