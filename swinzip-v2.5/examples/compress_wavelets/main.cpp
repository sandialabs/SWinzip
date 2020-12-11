#include <cstdio> 
#include <cstdlib> 
#include <iostream>
#include <fstream> 
#include <vector> 
#include <chrono>

#include <getopt.h>

#include "Alpert_Transform.hpp"
#include "Utils.hpp"

#include <boost/numeric/bindings/blas/blas.h>

using namespace std;

int main(int argc, char** argv) {
    
  if(argc == 1 ){
    cout << argv[0] << " usage is:" << endl;
    cout << argv[0] << " --reduction 'reduction_value' [mesh input_file] [data input_file]" << endl;
    cout << "If the data input file is provided it should have the same number of points as the mesh, if not a default data field will be assigned" << endl;
    cout << "Please refer to README.txt in this folder for more info" << endl;
    exit (0);
  }

    chrono::duration <double> duration;
    chrono::system_clock::time_point t1, t2;
    
    int d, N;
    double R=10.0e0;
    int c;
    
    {
          static struct option long_options[] =
        {
          /* These options set a flag. */
          /* These options don’t set a flag.
             We distinguish them by their indices. */
          {"reduction",  required_argument, 0, 'r'},
          {0, 0, 0}
        };
      /* getopt_long stores the option index here. */
      int option_index = 0;

      while (c = getopt_long (argc, argv, "r",
			      long_options, &option_index) != -1){

	if(option_index == 0 && optarg){
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
   // meshpoints 
    vector<vector<double> >meshPoints(N, vector<double>(d));
    
    // Read mesh into array
    for (int i = 0; i < N; i++)
        for (int j = 0; j < d; j++)
            f >> meshPoints[i][j];
    
    f.close();
    
    // Setting Wavelet orders
    vector<int> ki(d);
    for (int j = 0; j < d; j++)
      ki[j] = 5;
    
    boost::numeric::ublas::vector<double> x(N);
    boost::numeric::ublas::vector<double> w(N);
    boost::numeric::ublas::vector<double> xw(N);
    boost::numeric::ublas::vector<double> ys(N);
    
    // Setting data, this step can be replaced by some commands to read an external data file
    double pi=4.0e0*atan(1.0e0);
    if (argc>4) {
      cout << "Reading data input file" << argv[optind+1] << endl;
      f.open(argv[optind+1]);
      for (int i = 0; i < N; i++)
	f >> x(i);
    }
    else
      for (int i = 0; i < N; i++) x(i)=(4.0e0*sin(8.0e0*pi*meshPoints[i][0]))*(4.0e0*sin(7.0e0*pi*meshPoints[i][1]) )*3.0e0*sin(6.0e0*pi*meshPoints[i][0]);
//     for (int i = 0; i < N; i++) x[i]=(4.0*sin(2.0*pi*p[i][0]) - 4.0*sin(2.0*pi*p[i][1]))*3.0*sin(2.0*pi*p[i][0]);
    
    int J=-1;
    std::list<std::list<boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> > > Uj;
    std::list<std::vector<int> > part;
    std::vector<int> G;
    //void Compute_Uj(std::vector<std::vector<double> >meshPoints, std::vector<int> ki, int &J, std::list<std::list<boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> > > &Uj, std::list<std::vector<int> > &part, std::vector<int> &G) 
    
    // Forward wavelet transform
    t1 = chrono::system_clock::now();
    Compute_Uj(meshPoints,ki,J,Uj, part, G);
    Perform_Alpert_transform(x, ki, J,1, w, Uj, part, G);
    
    // Sorting, wavelet compression and inverse wavelet transform
    vector<double> yv(N);
    vector<size_t> iyv(N);
    for (int i = 0; i < N; i++) {yv[i]=fabs(w(i));iyv[i]=i;}
    sortandindices(yv, iyv);
    t2 = chrono::system_clock::now();
    duration=t2-t1;
    cout << "Total wavelet compression time " << duration.count() << endl;
    
//     for (int i = 0; i < N; i++)
//       std::cout << w(i) << std::endl;
    
    t1 = chrono::system_clock::now();
    for (int i = NR; i < N; i++) ys(iyv[i])=w(iyv[i]);
    
    // xw is the reconstruction by wavelets
    Perform_Alpert_transform(ys,  ki, J,-1, xw, Uj, part, G);
    
    // Computing NRMSE
    double ew=0.0e0;
    for (int i = 0; i < N; i++) 
      ew+=pow(x(i)-xw(i),2.0e0);
    
    t2 = chrono::system_clock::now();
    duration=t2-t1;
    cout << "Total wavelet decompression time " << duration.count() << endl;
    
    // Computing min and max values in the data x
    double mx=x(0);
    double mn=x(0);
    for (int i = 0; i < N; i++) {
      if (x(i)>mx)
	mx=x(i);
      if (x(i)<mn)
	mn=x(i);
    }
    
    ew=sqrt(ew/double(N))/(mx-mn);
    
    cout << "Alpert Wavelets Compression NRMSE = " << ew << endl;
}
