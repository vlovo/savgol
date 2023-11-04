

#include "savgol.hpp"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>

int main( )
{
  // comment PR01
  using namespace filter;

  std::ofstream file;

  std::vector<float> data{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  std::vector<float> gauss;
  std::vector<float> savgolOrder2;
  std::vector<float> average;

  const int windowSize = 5;

  savgol(data.begin(), data.end(), std::back_inserter(gauss), SmoothGaussian(windowSize));

  savgol(data.begin(), data.end(), std::back_inserter(savgolOrder2), SmoothQuadCubic(windowSize));

  savgol(data.begin(), data.end(), std::back_inserter(average), SmoothAverage(windowSize));

  file.open("savgol_demo.dat");

  int i = 0;
  for (auto e : data)
  {
    file << e << " " << average[i] << " " << gauss[i] << " " << savgolOrder2[i] << "\n";

    ++i;
  }

  file.close();

  return 0;
}
