#include <cstdio> 
#include <cstdlib> 
#include <iostream>
#include <math.h>
#include <random>
#include <chrono>

#include "Sampling_Matrix.hpp"

using namespace std;

void ConstructBernoulliMatrix(int M, int N, unsigned seed, boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> &P) {

    boost::numeric::ublas::vector<double> rowelems(N);
    for(int i = 0 ; i < N ; ++i) rowelems(i)=1.0e0/sqrt(double(N));
      
    for (int i=floor(double(N)/2.0);i<N;i++)
        rowelems(i)=-rowelems(i);
    
    std::default_random_engine generator (seed);
    uniform_int_distribution<int> distribution(0,N-1);
    
    for(int j = 0 ; j < M ; ++j)
        for(int i = 0 ; i < N ; ++i)
            P(i,j)=rowelems(distribution(generator));
    
}


void CompressBernoulli(std::vector<double> x, double * y, int M, unsigned seed) {
    
    int N=x.size();
    vector<double> rowelems(N,1.0/sqrt(double(N)));
    
    for (int i=floor(double(N)/2.0);i<N;i++)
        rowelems[i]=-rowelems[i];
    
    std::default_random_engine generator (seed);
    uniform_int_distribution<int> distribution(0,N-1);
    
    cerr << "Sampling " << M << "x" << N << " total random bits\n";

    for(int j = 0 ; j < M ; ++j) {
        y[j]=0.0e0;
        for(int i = 0 ; i < N ; ++i)
            y[j]+=x[i]*rowelems[distribution(generator)];
    }
}
