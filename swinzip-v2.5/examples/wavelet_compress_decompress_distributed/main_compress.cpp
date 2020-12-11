#include <cstdio> 
#include <cstdlib> 
#include <iostream>
#include <fstream> 
#include <vector> 
#include <cmath>

#include <getopt.h>

#include "Alpert_Transform.hpp"
#include "Utils.hpp"

#include "mpi.h"

using namespace std;

int main(int argc, char** argv) {
    
  if(argc == 1 ){
    cout << argv[0] << " usage is:" << endl;
    cout << argv[0] << "mpirun -np [number of processes] ./wcompress [threshold value] [path to data files+filename root] [number of physical dimensions] [wavelet order]" << endl;
    cout << "In this example, the mesh data is included in the data files as the first columns" << endl;
    cout << "The number of physical dimensions and wavelet order should be provided in the last two arguments, respectively" << endl;
    cout << "Please refer to README.txt in this folder for more info" << endl;
    exit (0);
  }
  
    MPI_Init(&argc, &argv);
    int rank, np; // Number of MPI processes, My process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    // Nv is the number of variables
    // N is the number of mesh points
    int dd, N,Nv;
    
    // d is the number of dimensions
    int d=atoi(argv[3]);
    double R;
    
    std::string fa = argv[2];
    std::string fb;
    char buffer [5];
    std::string txt=".txt";
    std::string fc;
    if (rank+1<10)
      sprintf(buffer,"000%d",rank+1);
    else
      sprintf(buffer,"00%d",rank+1);
    fc=fa+buffer+txt;
    
    // Reading the number of points and dd = d + Nv
    ifstream f(fc);
    f >> N >> dd;
    
    Nv=dd-d;
       
    vector<vector<double> > p(N, vector<double>(d));
    vector<vector<double> > data(N, vector<double>(Nv));
        
    // Read mesh and data into arrays
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < d; j++)
            f >> p[i][j];
	for (int j = 0; j < Nv; j++)
            f >> data[i][j];
    }
    
    f.close();
    
    // Setting Wavelet orders
    vector<int> ki(d);
    int kk=atoi(argv[4]);
    if ((kk<1) || (kk>6)) {
      cout << "Invalid wavelet order" << endl;
      exit (0); 
    }      
    
    for (int j = 0; j < d; j++)
      ki[j] = kk;
    
    boost::numeric::ublas::vector<double> x(N); // original data vector
    boost::numeric::ublas::vector<double> w(N); // transformed vector
    boost::numeric::ublas::vector<double> xw(N);// reconstructed data vector
        
    int J=-1;
    std::list<std::list<boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> > > Uj;
    std::list<std::vector<int> > part;
    std::vector<int> G;
    
    // reading the threshold
    double th0=atof(argv[1]);
    
    double mn,mx,mng,mxg;
    int M;
    
    FILE *mcw;
    fa = "mesh_";
    fb = ".dat";
    fc=fa+buffer+fb;
    mcw=fopen(fc.c_str(), "w");
    
    // First dump the mesh in the output file
    // because it is needed for decompression
    // the size overhead caused by the mesh amortizes
    // with more variables and time steps
    fwrite (&d , sizeof(int), 1, mcw);
    fwrite (&N , sizeof(int), 1, mcw);
    for (int i = 0; i < d; i++)
      fwrite (&ki[i] , sizeof(int), 1, mcw);
    for (int i = 0; i < N; i++)
      fwrite (&p[i][0], sizeof(double), d, mcw);
    fclose(mcw);
    
    // this files contains the compressed data
    FILE *fcw;
    fa = "fw_";
    fb = ".dat";
    fc=fa+buffer+fb;
    fcw=fopen(fc.c_str(), "w");
    fwrite (&Nv , sizeof(int), 1, fcw);
    
    // pre-computing the moment matrices once for all since the mesh does not change with time
    Compute_Uj(p,ki,J,Uj, part, G);

    // put in this file the resulting compression ratio and error
    FILE *fo;
    fa = "R_NRMSE_";
    fb = ".txt";
    fc=fa+buffer+fb;
    fo=fopen(fc.c_str(), "w");
        
    if (rank==0) {
      cout << "Number of dimensions = " << d << endl;
      cout << "Number of variables = " << Nv << endl;
      cout << "Wavelet order = " << kk << endl;
    }
    cout << "Processor " << rank << " , Number of Points = " << N << endl;
    
    double rr,Rg,ewg;
    for (int v=0;v<Nv;v++) {      
      
        for (int i = 0; i < N; i++)
	  x(i)=data[i][v];
	
	// Computing min and max values in the data x
	// We use percentiles to get rid of outliers
	prctile(&x[0], N, 97.0e0, &mx);
	prctile(&x[0], N, 3.0e0, &mn);
	
	// Computing the global maximum and minimum in the data
	MPI_Allreduce(&mn, &mng, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&mx, &mxg, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	
	// Forward wavelet transform
	Optimal_Alpert_transform(x, p, ki, J, w, Uj, part, G, 0.0e0, R,mx,mn,mxg,mng, th0);
	
	rr=1/R;
	MPI_Allreduce(&rr, &Rg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	// Collecting large wavelet coefficients
	vector<double> wc;
	vector<int> ic;
	for (int i = 0; i < N; i++) {
	  if (w(i)!=0.0e0) {
	    wc.push_back(w(i));
	    ic.push_back(i);
	  }
	}
	
	M=wc.size();
	
	fwrite (&M , sizeof(int), 1, fcw);       // Dump the number of large coefficients
	fwrite (&wc[0] , sizeof(double), M, fcw);// Dump the large coefficients (doubles)
	fwrite (&ic[0] , sizeof(int), M, fcw);   // Dump the large coefficient indices (integers)
	
	// Inverse wavelet transform
	Perform_Alpert_transform(w, ki, J,-1, xw, Uj, part, G);
	
	// Computing NRMSE
	double ew=0.0e0;
	for (int i = 0; i < N; i++) 
	  ew+=pow(x(i)-xw(i),2.0e0);
	
	if (mx==mn)
	  ew=0.0e0;
	else
	  ew=sqrt(ew/double(N))/(mxg-mng);
	
	fwrite (&ew , sizeof(double), 1, fcw);  // // Dump the resulting error
	MPI_Allreduce(&ew, &ewg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	// saving error and compression ratio to file
	fprintf(fo, "%d %e %e\n",v,R,ew);
	
	Rg=double(np)/Rg;
        ewg/=double(np);
        if (rank==0) cout << v << " " << Rg << " " << ewg << endl;
    }
        
    fclose(fcw);
    
    MPI_Finalize();
}
