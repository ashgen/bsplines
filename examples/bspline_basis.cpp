// Creating a set of b-spline basis and calling them. 

#include "../bsplines.hpp"
#include <iostream>
using namespace std;
#include "../../armadillo/include/armadillo"

int main(){

	// Construct a vector of breakpoints.  See bsplines.hpp file for alternative definition. 
	vector<double> breakpts {-0.06,-0.03,-0.02,-0.01,0.0,0.01,0.02,0.03,0.06};


	// Define order of b-spline
	int k=3;

	// Define bspline basis:
    bspline_basis mybasis(breakpts,k);
  typedef std::chrono::high_resolution_clock Clock;
  auto x=arma::randu(1e9);
  for(int i=0;i<1e9;i++){
    auto begin = Clock::now();
    auto t=mybasis.get_Bix(x[i]/100);
    auto end = Clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << endl;

  }
}
