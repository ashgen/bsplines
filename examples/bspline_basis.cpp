// Creating a set of b-spline basis and calling them. 

#include "../bsplines.hpp"
#include <iostream>
using namespace std;


int main(){

	// Construct a vector of breakpoints.  See bsplines.hpp file for alternative definition. 
	vector<double> breakpts {-0.06,-0.03,-0.02,-0.01,0.0,0.01,0.02,0.03,0.06};


	// Define order of b-spline
	int k=3;

	// Define bspline basis: 
	bspline_basis mybasis(breakpts,k); 
    for(int i=0;i < mybasis.get_nbasis();i++){
      auto t=mybasis.get_Bix(i,0.07);
      std::cout<<t<<"\t";

    }




}
