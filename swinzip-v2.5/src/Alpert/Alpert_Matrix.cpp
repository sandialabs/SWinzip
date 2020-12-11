#include <cstdio> 
#include <cstdlib> 
#include <iostream> 
#include <algorithm> 
#include <vector> 
#include <math.h> 
#include <list>
#include <iterator>
#include <numeric>
#include <functional> 

#include "Alpert_Matrix.hpp"
#include "Utils.hpp"

using namespace std;

void boostbuild_Alpert_matrix(vector<vector<double> >p, vector<int> ki,
			  boost::numeric::ublas::compressed_matrix <double> &  CWlast, 
			  int J,int Numnodes) {
  int d=(p[0]).size();
  if (d!=ki.size()) {
    cout << "Error: The size of vanishing moments array should be identical to the number of dimensions" << endl;
    abort();
  }
  
  int N = p.size();
  list<vector<int> > part;
  vector<int> G, cG, seli;
  int k2=1;
  for (int i = 0; i < ki.size(); i++)
        k2*=ki[i];


  dichotomic_grouping(p, k2, part, G);
    
  int P=G.size();
  
  // The required wavelet level is now supplied as an input
  // If we supply a negative level, the code computes and uses the maximum number of levels
  if (J<=0)
    J=int(log2(double(P))+1.0);

  cG=G;
  partial_sum(G.begin(), G.end(), cG.begin());
  vector<int> si(G.size()+1);
  si[0]=1;
  for (int i = 0; i < cG.size(); i++)
    si[i+1]=cG[i]+1;

  list<vector<int> >::iterator itp;
  itp=part.begin();
  
  list<boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> > Mi,Ui,Uti;
  list<boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> >::iterator itM, itU;

  // building a table for the vanishing moments orders in each direction
  std::vector< std::vector < int > > table;
  {
    std::vector < int > id;
    for (int i = 0; i < d; i++)
      id.push_back(0);//id[i]=0;
    
    int cnt=0;
    for (int i = 0; i < d; ++i){
      table.push_back(std::vector<int>(k2));
    }
    
    while (cnt<k2) {
      for (int j = 0; j < d; ++j){
	table[j][cnt] = id[j];
      }
      id[0]+=1;
      int i = 0;
      while ((i<(d-1))&&(id[i]==ki[i])) {
	id[i]=0;
	i+=1;
	id[i]+=1;
      }
      cnt+=1;
    }
  }

  // initialization of moments matrices  
  for (int i = 0; i < P; i++)  {
      
    seli = *itp;
    int kk = seli.size();                   // nbr of point in this bin, ie. 2k^2 in regular case

     // build moment matrix 
    boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> M (kk, 2*k2, 0.0e0);
    for (int ii = 0; ii < kk; ii++)  {
        // first part    : polynoms X^j1*Y^j2 for 0<=j1,j2<k
        // second part   : polynoms X^j1*Y^j2 for k<=j1<2k and 0<=j2<k
        for (int j = 0; j < k2; j++)  {
            M(ii,j)=1.0e0;
            for (int k = 0; k < d; k++)
                // compute polynomials at the point coordinates
                M(ii,j) *= pow(p[seli[ii]][k],table[d-k-1][j]);
        }
    }
        
    Mi.push_back(M);
    M=boost::numeric::ublas::project(M,boost::numeric::ublas::range(0,kk),boost::numeric::ublas::range(0,kk)); 
    // orthogonalize
    boost::numeric::ublas::vector<double> tau(kk);
    boost::numeric::bindings::lapack::geqrf(M,tau);
    boost::numeric::bindings::lapack::orgqr(M,tau);
    Uti.push_back(M);
    Ui.push_back(boost::numeric::ublas::trans(M));
    ++itp;
  }
  
//  ////////////////////////////////
//// initialization of first matrix
  boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> Uii, Mii, Vii;
  boost::numeric::ublas::generalized_vector_of_vector< double, boost::numeric::ublas::row_major, boost::numeric::ublas::vector<boost::numeric::ublas::compressed_vector<double>  > > bWlast(Numnodes, Numnodes,0);
  
  itp=part.begin();
  itU=Uti.begin();
  for (int i = 0; i < P; i++)  {
      seli = *itp;                   // selected column
      // to keep : upper part is of size n-P*k^2
      int offs = si[i]-i*k2;   // offset on row
      int lo = G[i]-k2;          // length on row
      Uii=*itU; 

      int kk = seli.size();

      for (int j = 0; j < lo; j++) {
	  for( int jj = 0; jj < kk; jj++){
	    bWlast(offs+j-1,seli[jj]) = Uii(jj,k2+j);
	  }
      }

      // to retransform : lower part is of size P*k2
      offs = N - P*k2+1 + i*k2;
      if (offs>=0) {
	for (int j = 0; j < k2; j++){
	  for( int jj = 0; jj < kk; jj++){
	    bWlast(offs+j-1,seli[jj]) = Uii(jj,j);
	  }
	}
      }
      else
	  cout << "Pbm empty bin." << endl;
      
	  // bug ...
      ++itU;
      ++itp;
  }
    
  CWlast.assign(bWlast);
	
  if ((J>1)||(J<=0)) {
    ////   /////////////////////////////////////////////////
      for (int j = 2; j <= J; j++) {   // for each scale
	
	// at this scale, we have nj = P/2^(j-1) groups
	int nnj=1;  
	for (int k = 1; k <= (j-1); k++)  nnj*=2; // n/(2^j*k);
	int nj=floor(double(P)/double(nnj));
	int mj = floor(double(P)/double(nnj)*double(2*k2));     // total length of the blocks
    //     nj=P/nj;
    //     int mj = nj*2*k2;     // total length of the blocks
	
	// update each sub matrix
	itU=Ui.begin();
	itM=Mi.begin();
	list<boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> > MMi,UUi,UUti;
	for (int i = 0; i < nj; i++) {
		  
	    boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> M (2*k2, 2*k2, 0.0e0);
	    Uii=*itU;
	    Mii=*itM;
	    int nc=Uii.size2();
	    Vii=boost::numeric::ublas::prod(boost::numeric::ublas::project(Uii,boost::numeric::ublas::range(0,k2),boost::numeric::ublas::range(0,nc)),Mii);
	    nc=Vii.size2();
	    
	    boost::numeric::ublas::project(M,boost::numeric::ublas::range(0,k2),boost::numeric::ublas::range(0,nc))=Vii;
	    ++itU;
	    ++itM;
	    Uii=*itU;
	    Mii=*itM;
	    nc=Uii.size2();
	    Vii=boost::numeric::ublas::prod(boost::numeric::ublas::project(Uii,boost::numeric::ublas::range(0,k2),boost::numeric::ublas::range(0,nc)),Mii);
	    nc=Vii.size2();
	    
	    boost::numeric::ublas::project(M,boost::numeric::ublas::range(k2,2*k2),boost::numeric::ublas::range(0,nc))=Vii;
	    MMi.push_back(M);
	    
	    // orthogonalize
	    boost::numeric::ublas::vector<double> tau(2*k2);
	    boost::numeric::bindings::lapack::geqrf(M,tau);
	    boost::numeric::bindings::lapack::orgqr(M,tau);
	    UUti.push_back(M);
	    
	    UUi.push_back(boost::numeric::ublas::trans(M));
		
	    ++itU;
	    ++itM;
	    }
	  
	Mi = MMi;
	Ui = UUi;
	Uti = UUti;
	
	// lower part of the multiplicative matrix
	boost::numeric::ublas::generalized_vector_of_vector< double, boost::numeric::ublas::row_major, boost::numeric::ublas::vector<boost::numeric::ublas::compressed_vector<double> > >  bUU(mj, mj);

	itU=Uti.begin();
	for (int i = 0; i < nj; i++) {
	  Uii=*itU;
	  std::vector<int> Mind(2*k2);
	  for (int k = 0; k < 2*k2; k++)
	      Mind[k]=2*k2*i+k;
	    
	  for (int k = 0; k < k2; k++){
	      for( int jj = 0; jj < 2*k2; jj++){
		bUU(k2*i+k,Mind[jj]) = Uii(jj,k+k2);
	      }
	  }


		
	  for (int k = 0; k < k2; k++){
	      for( int jj = 0; jj < 2*k2; jj++){
		bUU(mj/2+k2*i+k,Mind[jj]) = Uii(jj,k);
	      }
	  }
	    ++itU;
	}

	// multiplicative matrix

	boost::numeric::ublas::generalized_vector_of_vector< double, boost::numeric::ublas::row_major, boost::numeric::ublas::vector<boost::numeric::ublas::compressed_vector<double> > >  bUj(Numnodes, Numnodes);

	vector<double> id(1);id[0]=1.0e0;
	vector<int> idi(1);
	for (int i = 0; i < (N-mj); i++) { 
	  idi[0]=i; 
	  bUj(i,i) = 1.0e0;
	}

	{
	  std::vector < int > nind(mj);
	  std::vector < double > doo(mj);
	  int no;
	  
	  for (int i = (N-mj); i < N; i++) {      
	    for (int jj = 0; jj < no; jj++) nind[jj]+=(N-mj); 
	    for( int ll = 0; ll < mj;ll++){
	      if (bUU(i+mj-N,ll)!=0.0e0)
		bUj(i,ll+N-mj) = bUU(i+mj-N,ll);//doo[ll];
	    }
	  }
	}

// update vectors
        boost::numeric::ublas::compressed_matrix <double> CUj(Numnodes, Numnodes);
        CUj.assign(bUj);

	boost::numeric::ublas::compressed_matrix <double> CWcurr(Numnodes, Numnodes);
	
	cout << "j= " << j << endl;
	boost::numeric::ublas::sparse_prod(CUj,CWlast,CWcurr);

	CWlast = CWcurr;

      }
  }
 
}

void build_Alpert_Sampling_product(boost::numeric::ublas::compressed_matrix <double> U, boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> P, boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> &T) {
  
  int M=P.size2();

#ifdef OPENMP
#pragma omp parallel for
#endif
  for (int j = 0; j < M; j++) {
     boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> > mcP (P, j);
     boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> > mcT (T, j);
     boost::numeric::ublas::axpy_prod(U,mcP,mcT,true);
  }
  
}
