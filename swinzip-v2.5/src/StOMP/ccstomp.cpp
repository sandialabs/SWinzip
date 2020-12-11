#undef VERBOSE
#include <set>
#include <math.h>
#include <iostream>
#include "ccstomp.hpp"

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
	  double *aresult)
  { 
// Initialize threshold parameters

    int active_size = 0;
    int new_active_size = 0;
    double * Ared = NULL;
    double * Ared2 = NULL;
    double normy, my, thr;
    int iter = 1;
    double * res = new double[ressize[0]*ressize[1]];
    double * newJ = new double[asize[1]*ressize[1]];
    double normres;
    
    int ione;
    int M, N, K, LWORK, LDB, MM, NN, inc;
    double alpha, beta, malpha, alphar;
    M = asize[0];
    K = ressize[1];
    double * WORK;

    NN = asize[1];
    MM = asize[0];
    alpha = 1.0e0;
    malpha = -1.0e0;
    beta = 0.0e0;
    inc =1;
    
    normy = BLAS_DNRM2(&M, y, &inc);
    my = BLAS_DASUM(&M, y, &inc)/double(M);
    thr=3.0e0*normy/my/sqrt(double(asize[1]))+0.03e0;
    thr/=sqrt(double(asize[0]));  
    
    for(int ict = 0; ict < ressize[0]*ressize[1];ict++)res[ict] = y[ict];
// % Initialize

    int done = 0;
    double * corr = new double[asize[1]*ressize[1]];
    
    normres=BLAS_DNRM2(&MM, res, &inc);
    
    while (!done){
            
      for(int ict  = 0; ict < asize[1]*ressize[1]; ict ++) corr[ict] =0.0e0;
      
      alphar = 1.0e0/normres;
      BLAS_DGEMV("N", &NN, &MM, &alphar, a, &NN, res, &inc, &beta, corr, &inc);
      
      new_active_size = cc_find_new_j(corr,
              asize[1]*ressize[1],
              activeSet,
              active_size,
              thr,
              newJ);
      
      if(active_size == new_active_size){
        done = 1;
      }
      
      else{ /*not done*/
        active_size = new_active_size;
        if(activeSet)   delete[]activeSet;
        activeSet = new double[new_active_size];
	BLAS_DCOPY( &new_active_size, newJ, &inc, activeSet, &inc);
        N = new_active_size;
	
//         % Compute current estimate and residual
        if(Ared)delete [] Ared;        
        Ared = new double[asize[0]*new_active_size];

	for(int jct = 0; jct < new_active_size; jct ++)
	  BLAS_DCOPY( &M, a+((int)activeSet[jct]-1), &NN, Ared+jct*asize[0], &inc);

        if(aresult)delete [] aresult;        
	LWORK=2*std::min(M,N);	
	WORK = new double[LWORK];
	LDB=std::max(M,N);
        aresult  = new double [LDB];            
	
	BLAS_DCOPY( &M, y, &inc, aresult, &inc);
	
	if(Ared2)delete [] Ared2;
	Ared2 = new double[M*N];
	N=N*M;
        BLAS_DCOPY( &N, Ared, &inc, Ared2, &inc);
	N=N/M;
        LAPACK_DGELS("N", &M, &N, &K, Ared2, &M, aresult, &LDB, WORK, &LWORK, &ione);
	delete[]WORK;
	
        
	BLAS_DCOPY( &M, y, &inc, res, &inc);
	BLAS_DGEMV("N", &M, &N, &malpha, Ared, &M, aresult, &inc, &alpha, res, &inc);
              
      }
      
      normres=BLAS_DNRM2(&MM, res, &inc);
      
      iter = iter+1;
      std::cout << iter << " " << new_active_size << " " << normres << std::endl;
      
      // Check stopping criteria
      if ((normres <= OptTol*normy)||(iter > maxIters)) 
        done = 1;      
      
    }
    for(int ict = 0; ict < new_active_size;ict ++)
          asol[(int)activeSet[ict]-1] = aresult[ict];
    
    *niter=iter;
    *naSet=new_active_size;
	
    delete [] Ared;
    delete [] aresult;
    delete [] newJ;
    delete [] res;
    delete [] corr;
    
  }
  
int cc_find_new_j(double val [],
          int vct,
          double as [],
          int nas,
          double tol,
          double res []) {
    
    int ict;
    std::set<int>Aset;
    for(ict = 0; ict < nas;ict ++)Aset.insert(as[ict]);

    for(ict = 0; ict < vct; ict ++)
      if(fabs(val[ict])>tol)
	Aset.insert(ict+1);
    
    std::set<int>::iterator sit;
    ict = 0;
    for(sit = Aset.begin() ; sit != Aset.end(); sit++){
      res[ict] = *sit;
      ict ++;
    }
    return ict;
  }
  