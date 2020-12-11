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

void boostbuild_Alpert_matrix(std::vector<std::vector<double> >p, std::vector<int> ki, boost::numeric::ublas::compressed_matrix <double> &  bU , int J,int N);
void build_Alpert_Sampling_product(boost::numeric::ublas::compressed_matrix <double> U, boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> P, boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> &T);
