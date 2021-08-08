# savgol

C++ implementation of Savitzky Golay  Filter based on the original paper : 
"Smoothing and Differentiation of Data by Simplified Least Squares Procedures" by  A.Savitzky and M. J.E Golay". [1]

* The implementation is close to the original publication. The convolution weights for each given window size are stored
in a constexpr map.   

* there is no padding implemented so far. The data points at the bounds left unchanged.
* 
## usage

```cpp
#include <vevtor>
#include "savgol.hpp"

std::vector<double> data = {-0.2,0.4,0.7,1.2,-0.1,-0.2,-0.2,-0.2,-0.2,-0.2};

std::array<double> filtered_data;

filter::savgol(data.begin(), data.end(), sd::back_inserter(filtered_data),filter::SmoothQuadCubic(7) );_

```

## implementation

###  available Convolution Operators

> SmoothQuadCubic : 

> SmoothQuarticQuintic :  

> DeriveQuadFirst :
###  available window sizes

* as presentented in the paper the window sizes are odd numbers, up to a size  of 25 


[1] Abraham. Savitzky, M. J. E. Golay: Smoothing and Differentiation of Data by Simplified Least Squares Procedures. In: Analytical Chemistry. Band 36, Nr. 8, 1. Juni 1964, S. 1627–1639, [doi:10.1021/ac60214a047](doi:10.1021/ac60214a047).