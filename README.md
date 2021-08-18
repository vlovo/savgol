# savgol

C++ implementation of Savitzky Golay  Filter based on the original paper : 
"Smoothing and Differentiation of Data by Simplified Least Squares Procedures" by  A.Savitzky and M. J.E Golay". [1]

* The implementation is close to the original publication with corrections taken from [2] . The convolution weights for each given window size are stored
in a constexpr map.   

* there is no padding implemented so far. The data points at the bounds left unchanged.
* 
## usage

```cpp
#include <vector>
#include "savgol.hpp"

std::vector<double> data = {-0.2,0.4,0.7,1.2,-0.1,-0.2,-0.2,-0.2,-0.2,-0.2};

std::vector<double> filtered_data;

filter::savgol(data.begin(), data.end(), sd::back_inserter(filtered_data),filter::SmoothQuadCubic(7) );

```

## implementation

###  available Convolution Operators


> SmoothQuadCubic : 

> SmoothQuarticQuintic :  

> DeriveQuadFirst :

--- 

> SmoothGaussian: to make comparisions to SG filter, the gaussian kernel is approximated by binominal coefficients, taken from Pascals triangle.

> SmoothAverage: to make comparisions to SG filter,  convolution weights are constant  =1.
> 
> 
###  available window sizes

* as presentented in the original paper the window sizes are odd numbers, up to a size  of 25 



![](filter_response.svg) 

[1] A. Savitzky, M. J. E. Golay: Smoothing and Differentiation of Data by Simplified Least Squares Procedures. In: Analytical Chemistry. Volume 36, No. 8, June 1964,  p. 1627–1639, [doi:10.1021/ac60214a047](doi:10.1021/ac60214a047).

[2] J. Steiner, Y. Termonia, J. Deltour: Comments on Smoothing And Differentiation Of Data By Simplified Least Square Procedure In: Analytical Chemistry. Volume 44, No. 11, 1. Sept 1972, p. 1906–1909