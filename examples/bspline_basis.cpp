// Creating a set of b-spline basis and calling them. 

#include "../bsplines.hpp"
#include "../splines2Armadillo.h"
#include <iostream>
using namespace std;

int main(){

  // Construct a vector of breakpoints.  See bsplines.hpp file for alternative definition.
  arma::vec breakpts {-0.03,-0.02,-0.01,0.0,0.01,0.02,0.03};
  arma::vec breakptsB {-0.06,0.06};
  auto x=0.003;
  splines2::NaturalSpline nat_spline;
  nat_spline.set_internal_knots(breakpts);
  nat_spline.set_boundary_knots(breakptsB);
  vector<double> breakptsf {-0.06,-0.03,-0.02,-0.01,0.0,0.01,0.02,0.03,0.06};


  // Define order of b-spline
  int k=4;
  bspline_basis mybasis(breakptsf,k);
  std::cout<<mybasis.basis_vector(x).t()<<std::endl;
  nat_spline.set_x(x);
  std::cout<<nat_spline.basis(true)<<std::endl;
}
