#include <iostream>
#include <vector>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>

void ConstructBernoulliMatrix(int M, int N, unsigned seed, boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> &P);

void CompressBernoulli(std::vector<double> x, double * y, int M, unsigned seed);
