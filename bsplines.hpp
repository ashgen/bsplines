/*

    Copyright (c) F.I.Diakogiannis  2015


    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with popmcmc++. If not, see <http://www.gnu.org/licenses/>.

*/


/* File Description:
	
This file contains definition of the bspline basis class and bspline functions (pp polynomials). 
Sofisticated algorithms are based mainly in the book "The NURBS Book" - Piegl, Tiller, 2nd edition. 


	References: 
		[1] The NURBS Book, Piegl - Tiller, 1996, 2nd edition 
		[2] Curves and surfaces for CAGD, G. Farin, 2002. 
*/
 

#ifndef _bspline_basis_
#define _bspline_basis_

#include <iostream> 
#include <utility>
#include <algorithm>
#include <vector>
#include <cmath>
#include "macros.h"
#include <functional>
#include "armadillo"
using namespace std; 




class bspline_basis{

	private:
	//public:
		int k; 						/*! order of Bspline basis */
		int nbreak; 					/*! Dimension of breakpoints vector. */
		int nknots; 					/*! Dimension of knots vector. */
		int nbasis; 					/*! Number of basis functions */
		int nderivs; 					/*! Total number of derivatives evaluated. Allow (nderivs < k) ?  */

/*add*/		int nintegrals; 				/*! Total number of evaluated integrals. Default: 1 */		

		/* Single knot <--> multiplicity =1 */
		vector<pair<double,int> > knots_wth_mult; 	/*! Vector of knots, each entry is the pair (knot value, multiplicity)*/
		vector<double> breakpts;			/*! Represents strictly increasing values of knots, excluding all multiplicities */
		vector<double> knots; 				/*! Raw knot vector of BSpline basis, it may contain multiple entries due to multiplicity */
		
		bool greville_evaluated; 				/*! Logical operator, informs if Greville abscissa have been evaluated */
		vector<double> greville; 			/* Contains values of Greville abscissa according to definitions [2] */

		arma::vec Bix_nonzero; 			/*! size: (k,1). Stores  nonzero components of bspline basis */
		arma::vec Bix; 				/*! size: nbasis Stores  all  components of bspline basis.  Not necessary - remove? */
		
		arma::vec DjBix_nonzero; 			/*! size: (nderivs+1,k). Stores  the jth derivative of the ith bspline basis function at x, nonzero terms only */
		arma::vec DjBix; 				/*! Stores  the jth derivative of the ith bspline basis function at x */

        arma::vec Bix_lower; 				/*! Stores  the jth derivative of the ith bspline basis function at x */
        arma::vec Bix_upper; 				/*! Stores  the jth derivative of the ith bspline basis function at x */
        arma::vec DBix_lower; 				/*! Stores  the jth derivative of the ith bspline basis function at x */
        arma::vec DBix_upper; 				/*! Stores  the jth derivative of the ith bspline basis function at x */


		// Need custom copy and assignment constructor for this, I think. 
		int i_saved; // initialize to negative; 					/*! Saved index i of Bix evaluation, avoids repetetion of find_nonzero_interval */ 
		double x_saved; // initialize to NAN; 						/*! Temporary stored x value, used for avoiding unnecessary function calls */
		bool eval_Bix_flag; 	// Informs program if previous calls to get_Bix took place 
		bool eval_DjBix_flag;   // Informs program if previous calls to get_DjBix took place

		/* Experimental */
		bool eval_DjBix_lim_flag; 	// Informs program if previous calls to get_DjBix_lim took place


	public:
		/* Experimental */
		int find_knot_span_of_x_lim(const double &x, const char * lim); 	/*! Returns integer i: t_i < x <= t_{i+k}. Upon call it stores i in i_saved */		
		const char * low = "low";
		const char * high = "high";



		int find_knot_span_of_x(const double &x); 					/*! Returns integer i: t_i <= x < t_{i+k}. Upon call it stores i in i_saved */		
		int find_knot_index (const double &x );						/*! Returns integer i that corresponds to the first instance of x in knot vector, t_i==x */
		pair<int,int> find_nonzero_basis_at_x(const double &x); 			/*! Returns first, last index of nonzero basis B_i(x) at particular x. */
		pair<int,int> find_base_nonzero_interval(const double &x); 			/*! Returns first (i) , last (i+k) index of knots t_i at particular x. */


	private:
		/* !ESSENTIAL ROUTINES FOR EVALUATION! Add as optional argument another knot vector for use in evaluation of integrals */
		void eval_nonzero_basis(const int &i, const double &x); 	/*! Evaluates non zero basis functions at x */
		void eval_nonzero_Djbasis(const int &i, const double &x);	/*! Evaluates non zero basis functions and their derivatives at x */

		void eval_Bix(const int &i, const double &x);			/*! Evaluates all basis functions at x */

		int idx_nonzero(const int &i, const int &j);			/*! collective index for DjBix_nonzero 2D --> 1D */
		int idx_all(const int &i, const int &j); 			/*! collective index for DjBix 2D --> 1D */
		void eval_DjBix(const int &i, const double &x);				/*! Evaluates all basis functions and their derivatives at x */
        void eval_DBix(const int &i, const double &x);				/*! Evaluates all basis functions and their derivatives at x */


	public:
		bspline_basis(){};
		virtual ~bspline_basis(){};
		/*! Constructor from custom-multiplicity knot vector */
		bspline_basis(const vector<pair<double,int> > &_knots_wth_mult, const int &_k);
		/*! Default clamped  knot vector constructor */
		bspline_basis(const vector<double> &_breakpts, int _k);

        void initialize(const vector<double> &_breakpts, int _k);
		//bspline_basis(const vector<double> &_knots, const int &_k);	/*! Need to add this constructor */

		void set_nderivs (const int & new_nderivs);			/*! Sets total number of derivatives to be evaluated. Default nderivs = min(3,k-1) */
		void set_nintegrals(const int & new_nintegrals); 		/*! Sets total number of Bix integrals. Default = 1 */


		int get_nbasis();						/*! Returns number of Bix basis functions */
		int get_nbreak();						/*! Returns dimension of breakpoints vector. */
		int get_order();						/*! Returns number of Bix basis functions */
		vector<pair<double,int> >  get_knots_wth_mult();		/*! Returns knot vector with multiplicities */
		vector<double> 	get_knots();					/*! Returns knot vector in vector<double> format, with possible multiple knots */
		vector<double> 	get_breakpts();					/*! Returns breakpoints vector */
		vector<double> get_abscissa();					/*! If necessary, evaluates and returns vector of Greville abscissa */


		void knots_TO_knots_wth_mult();					/*! Pass information from vector<double> to vector<pair<double,int> > */
		void knots_wth_mult_TO_knots();					/*! Pass information from vector<pair<double,int> > to vector<double> */
		void set_knots(vector<double> & knots_new);			/*! Update knot vector in vector<double> form  */
		void set_knots(vector<pair<double,int> > & knots_wth_mult_new);	/*! Update knot vector in vector<pair<double,int>> form  */



		/* Experimental */
		double get_DjBix_lim(int j, int i, const double &x, const char * lim); 	/*! Value d^j B_i(x)/ dx^j for upper or lower limit.  */

		/* Evaluation functions */
		double get_Bix(const int &i, const double &x); 			/*! Value B_i(x) */
		double get_DjBix(int j,int i, const double &x); 	/*! Value d^j B_i(x)/ dx^j */
		void eval_greville();						/*! Evaluates all Greville abscissa according to [2] */


		// TO BE ADDED, requires vectors of fixed coefficients to be stored for faster evaluation.
		double get_IntjBix(const int &j, const int &i, const double &x); /*! Value of j nested integrals \int_0^x B_i(u) du */

        const arma::vec &get_Bix(const double &x);
        const vector<double> &basis(const double &x);
        const arma::vec &eval(const int &ii, const double &x);
        const arma::vec &get_DjBix(int j, const double &x);
        const arma::vec &get_DBix(const double &x);
        arma::vec basis_vector(const double &x);
        arma::mat basis_matrix(const arma::vec &x);
};



#endif
