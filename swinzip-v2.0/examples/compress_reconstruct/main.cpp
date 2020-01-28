#include <cstdio> 
#include <cstdlib> 
#include <iostream>
#include <fstream> 
#include <vector> 
#include <chrono>

#include <getopt.h>

#include "Sampling_Matrix.hpp"
#include "Alpert_Matrix.hpp"
#include "Utils.hpp"
#include "ccstomp.hpp"

#include <boost/numeric/bindings/blas/blas.h>

using namespace std;

int main(int argc, char** argv) {
    
  if(argc == 1 ){
    cout << argv[0] << " usage is:" << endl;
    cout << argv[0] << " --seed 'seed_value' --reduction 'reduction_value' [mesh input_file] [data input_file]" << endl;
    cout << "If the data input file is provided it should have the same number of points as the mesh, if not a default data field will be assigned" << endl;
    cout << "Please refer to README.txt in this folder for more info" << endl;
    exit (0);
  }


  /*Seed is overwritten if set as a command line option*/
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    chrono::duration <double> duration;

    int d, N;
    double R=10.0e0;
    int c;
    double alpha=1.0e0, beta=0.0e0;
    int inc=1;
    
    {
          static struct option long_options[] =
        {
          /* These options set a flag. */
          /* These options don’t set a flag.
             We distinguish them by their indices. */
          {"seed",  required_argument, 0, 's'},
          {"reduction",  required_argument, 0, 'r'},
          {0, 0, 0, 0}
        };
      /* getopt_long stores the option index here. */
      int option_index = 0;

      while (c = getopt_long (argc, argv, "s",
			      long_options, &option_index) != -1){

	if(option_index == 0 && optarg){
	  seed  = atoi(optarg); 
	  cout << "Seed set to " << seed << endl;
	}
	if(option_index == 1 && optarg){
	  R  = atof(optarg); 
	  cout << "Reduction set to " << R << endl;
	}
      }

    }


    if(optind >= argc ) {
      cout << "Expected an input file name " << endl;
      exit(0);
    }       
    
    cout << "Reading mesh input file" << argv[optind] << endl;
    ifstream f(argv[optind]);
    f >> N >> d;

    cout << "Number of Points " << N << endl;
    cout << "Number of dimensions " << d << endl;
    cout << "R " << R << endl;
        
    int M=ceil(double(N)/R);
    int NR=N-M;
    
    vector<vector<double> >p(N, vector<double>(d));
    
    // Read mesh into array
    for (int i = 0; i < N; i++)
        for (int j = 0; j < d; j++)
            f >> p[i][j];
    
    f.close();
    
    // Setting Wavelet orders
    vector<int> ki(d);
    for (int j = 0; j < d; j++)
      ki[j] = 5;
    
    cout << "Constructing Alpert Matrix" << endl;
    boost::numeric::ublas::compressed_matrix <double>  U(N, N);
    chrono::system_clock::time_point t1 = chrono::system_clock::now();
    boostbuild_Alpert_matrix(p, ki,  U, -1, N);
    chrono::system_clock::time_point t2 = chrono::system_clock::now();
    duration=t2-t1;
    cout << "Total wavelet matrix construction time " << duration.count() << endl;
    
    boost::numeric::ublas::vector<double> x(N);
    boost::numeric::ublas::vector<double> w(N);
    boost::numeric::ublas::vector<double> xw(N);
    boost::numeric::ublas::vector<double> ys(N);
    
    // Setting data, this step can be replaced by some commands to read an external data file
    double pi=4.0e0*atan(1.0e0);
    if (argc>6) {
      cout << "Reading data input file" << argv[optind+1] << endl;
      f.open(argv[optind+1]);
      for (int i = 0; i < N; i++)
	f >> x(i);
    }
    else
      for (int i = 0; i < N; i++) x(i)=(4.0e0*sin(8.0e0*pi*p[i][0]))*(4.0e0*sin(7.0e0*pi*p[i][1]) )*3.0e0*sin(6.0e0*pi*p[i][0]);
//     for (int i = 0; i < N; i++) x[i]=(4.0*sin(2.0*pi*p[i][0]) - 4.0*sin(2.0*pi*p[i][1]))*3.0*sin(2.0*pi*p[i][0]);
    
    // Forward wavelet transform
    t1 = chrono::system_clock::now();
    boost::numeric::ublas::axpy_prod(U, x, w, true);    
    
    // Sorting, wavelet compression and inverse wavelet transform
    vector<double> yv(N);
    vector<size_t> iyv(N);
    for (int i = 0; i < N; i++) {yv[i]=fabs(w(i));iyv[i]=i;}
    sortandindices(yv, iyv);
    for (int i = NR; i < N; i++) ys(iyv[i])=w(iyv[i]);
    t2 = chrono::system_clock::now();
    duration=t2-t1;
    cout << "Total wavelet compression time " << duration.count() << endl;
    
    // xw is the reconstruction by wavelets
    boost::numeric::ublas::axpy_prod(ys, U, xw, true);
    
    cout << "Constructing Bernoulli Sampling Matrix" << endl;
    boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> P(N,M,0.0e0);
    t1 = chrono::system_clock::now();
    ConstructBernoulliMatrix(M, N, seed, P);
    t2 = chrono::system_clock::now();
    duration=t2-t1;
    cout << "Total compression matrix construction time " << duration.count() << endl;
    
    cout << "Performing CS compression" << endl;
    boost::numeric::ublas::vector<double> y(M);
    t1 = chrono::system_clock::now();
    BLAS_DGEMV("T", &N, &M, &alpha, &P(0,0), &N, &x(0), &inc, &beta, &y(0), &inc);
    t2 = chrono::system_clock::now();
    duration=t2-t1;     
    cout << "Total CS compression time " << duration.count() << endl;

    cout << "Constructing Alpert-Bernoulli Matrix product" << endl;
    t1 = chrono::system_clock::now();
    build_Alpert_Sampling_product(U,P,P);
    t2 = chrono::system_clock::now();
    duration=t2-t1;     
    cout << "Alpert-Compression matrix product time " << duration.count() << endl;

    // Computing incoherence between sampling and wavelet matrices
    double maxT=0.0;
    for (int j = 0; j < M; j++) 
      for (int i = 0; i < N; i++) 
	if (fabs(P(i,j))>maxT)
	  maxT=fabs(P(i,j));
	
    
    cout << "Performing StOMP" << endl;
    t1 = chrono::system_clock::now();
    boost::numeric::ublas::vector<double> ywp(N);
    int * asize= new int[2]; asize[0]=M; asize[1]=N;
    int * rsize= new int[2]; rsize[0]=M; rsize[1]=1;
    int niter=0, nywp=0;
    double *activeSet=NULL;
    double *aresult=NULL;
    // activeSet and aresult are the initial guess of the solution, here we assume the coefficient
    // vector is empty
    cc_solve_stomp(&P(0,0), &y(0), asize, rsize, 1.0e-8, 50, &niter, &nywp, &ywp(0), activeSet, aresult);
    
    boost::numeric::ublas::vector<double> xs(N);
    boost::numeric::ublas::axpy_prod(ywp, U, xs, true);
    t2 = chrono::system_clock::now();
    duration=t2-t1;     
    cout << "StOMP time " << duration.count() << endl;
    
    // Computing NRMSE
    double es=0.0e0, ew=0.0e0;
    for (int i = 0; i < N; i++) {
      ew+=pow(x(i)-xw(i),2.0e0);
      es+=pow(x(i)-xs(i),2.0e0); 
    }
    
    // Computing min and max values in the data x
    double mx=x(0);
    double mn=x(0);
    for (int i = 0; i < N; i++) {
      if (x(i)>mx)
	mx=x(i);
      if (x(i)<mn)
	mn=x(i);
    }
    
    es=sqrt(es/double(N))/(mx-mn);
    ew=sqrt(ew/double(N))/(mx-mn);
    
    cout << "Alpert Wavelets Compression NRMSE = " << ew << endl;
    cout << "CS Compression NRMSE = " << es << endl;
}
