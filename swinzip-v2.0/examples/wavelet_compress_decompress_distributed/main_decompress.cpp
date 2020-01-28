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
    cout << argv[0] << "mpirun -np [number of processes] ./wdecompress [mesh filename root] [compressed data filename root] [reconstructed data filename root]" << endl;
    cout << "In this example..." << endl;
    cout << "the mesh and compressed data files are obrained from wcompress" << endl;
    cout << "Please refer to README.txt in this folder for more info" << endl;
    exit (0);
  }
  
    MPI_Init(&argc, &argv);
    int rank, np; // Number of MPI processes, My process ID
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    // Nv is the number of variables
    // N is the number of mesh points
    // d is the number of dimensions
    int d, N,Nv,kk;
    double R;
    
    std::string fa = argv[1];
    std::string fb = ".dat";
    char buffer [5];
    std::string fc;
        
    if (rank+1<10)
      sprintf(buffer,"000%d",rank+1);
    else
      sprintf(buffer,"00%d",rank+1);
    
    FILE *mcw;       
    fc=fa+buffer+fb;
    mcw=fopen(fc.c_str(), "r");
        
    // Read the N and d from the mesh file
    fread (&d , sizeof(int), 1, mcw);
    
    // Setting Wavelet orders
    vector<int> ki(d);
    
    fread (&N , sizeof(int), 1, mcw);
    for (int i = 0; i < d; i++)
      fread (&ki[i] , sizeof(int), 1, mcw);
    kk=ki[0];
    
    cout << "Processor " << rank << " , Number of Points = " << N << endl;
    
    // read the mesh
    vector<vector<double> > p(N, vector<double>(d));    
    for (int i = 0; i < N; i++)
      fread(&p[i][0], sizeof(double), d, mcw);
    fclose(mcw);
    
    // read the number of variables from compressed data
    FILE *fcw;
    fa = argv[2];
    fb = ".dat";
    fc=fa+buffer+fb;
    fcw=fopen(fc.c_str(), "r");    
    fread (&Nv , sizeof(int), 1, fcw);
    
    vector<vector<double> > data(N, vector<double>(Nv));
    
    boost::numeric::ublas::vector<double> w(N); // transformed vector
    boost::numeric::ublas::vector<double> xw(N);// reconstructed data vector
        
    int J=-1;
    std::list<std::list<boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> > > Uj;
    std::list<std::vector<int> > part;
    std::vector<int> G;
        
    double ewg;
    int Ng,NR,NRg,M;
    
    if (rank==0) {
      cout << "Number of dimensions = " << d << endl;
      cout << "Number of variables = " << Nv << endl;
      cout << "Wavelet order = " << kk << endl;
    }
        
    // file write the reconstructed data
    fa = argv[3];
    fb = ".txt";
    fc=fa+buffer+fb;
    ofstream f(fc);
    f << N << " " << Nv << endl;
    
    // pre-computing the moment matrices once for all since the mesh does not change with time
    Compute_Uj(p,ki,J,Uj, part, G);

    for (int v=0;v<Nv;v++) {
        
        for(int i = 0; i < N; i++)
	  w(i)=0.0e0;
	
	fread (&M , sizeof(int), 1, fcw);
	
	// reading wavelet coefficients and their indices
	vector<double> wc(M);
	vector<int> ic(M);
	fread (&wc[0] , sizeof(double), M, fcw);
	fread (&ic[0] , sizeof(int), M, fcw);
	
	for(int i = 0; i < M; i++)
	  w(ic[i])=wc[i];
		
	// Inverse wavelet transform
	Perform_Alpert_transform(w, ki, J,-1, xw, Uj, part, G);
	for(int i = 0; i < N; i++)
	  data[i][v]=xw(i);
	  
	double ew;
	fread (&ew , sizeof(double), 1, fcw);
	R=double(N)/double(M);
	
	// computing global error and compression ratio
	MPI_Allreduce(&N, &Ng, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	NR=floor(double(N)/R);
	MPI_Allreduce(&NR, &NRg, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	ew=ew*ew*double(N);
	MPI_Allreduce(&ew, &ewg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	if (rank==0) cout << v << " " << double(Ng)/double(NRg) << " " << sqrt(ewg/double(Ng)) << endl;  
    }
    
    // Write reconstructed data to file
    for (int i = 0; i < N; i++) {
	for (int j = 0; j < Nv; j++)
            f << data[i][j] << " ";
	f << endl;
    }
    
    fclose(fcw);
    f.close();
    
    MPI_Finalize();
}
