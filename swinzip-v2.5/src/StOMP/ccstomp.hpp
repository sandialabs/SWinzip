#include <iostream>
#include <chrono>
#include <boost/numeric/bindings/lapack/lapack.h>
#include <boost/numeric/bindings/blas/blas.h>
#include <algorithm>

  int cc_find_new_j(double val [],
          int vct,
          double as [],
          int nas,
          double tol,
          double res []);
  
  
  void cc_solve_stomp(
          double * a,
          double * y,
          int * asize,
          int * ressize,
	  double OptTol,
	  int maxIters,
          int *niter,
          int *naSet,
          double *asol,
	  double *activeSet,
	  double *aresult);
   
  double cmeana_ms(double * vals, int size,double * res);
  