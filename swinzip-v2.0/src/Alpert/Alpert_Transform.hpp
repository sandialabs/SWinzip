#include <cstdio> 
#include <cstdlib> 
#include <iostream>
#include <vector>
#include <list>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/lapack/geqrf.hpp> 
#include <boost/numeric/bindings/lapack/orgqr.hpp> 
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector_of_vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>

#ifdef OPENMP
#include <omp.h>
#endif

void Perform_Alpert_transform(boost::numeric::ublas::vector<double> f, std::vector<int> ki, int J,int dir, boost::numeric::ublas::vector<double>& w, std::list<std::list<boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> > > Uj, std::list<std::vector<int> > part, std::vector<int> G);
void Compute_Uj(std::vector<std::vector<double> >p, std::vector<int> ki, int &J, std::list<std::list<boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> > > &Uj, std::list<std::vector<int> > &part, std::vector<int> &G);
void Optimal_Alpert_transform(boost::numeric::ublas::vector<double> f, std::vector<std::vector<double> >p, std::vector<int> ki, int J, boost::numeric::ublas::vector<double>& w, std::list<std::list<boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> > > Uj, std::list<std::vector<int> > part, std::vector<int> G, double NRMSE, double &R, double mx, double mn, double mxg, double mng, double alpha);