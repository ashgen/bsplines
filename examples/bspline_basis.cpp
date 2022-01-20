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
    int N=1000000;
	arma::vec test(N,arma::fill::randn);
    bspline_basis mybasis(breakpts,k);
    arma::vec coefs(10,arma::fill::randn);
    for(int i=0;i<N;i++){

        auto start = std::chrono::high_resolution_clock::now();
        auto t=arma::dot(mybasis.basis_vector(test[i]),coefs);
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed_seconds=std::chrono::duration_cast<std::chrono::nanoseconds>(
                end - start);
        std::cout<<elapsed_seconds.count()<<std::endl;
    }
}
