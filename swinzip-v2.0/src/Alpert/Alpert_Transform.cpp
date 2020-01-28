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

#include "Alpert_Transform.hpp"
#include "../Utils/Utils.hpp"

void Optimal_Alpert_transform(boost::numeric::ublas::vector<double> x, std::vector<std::vector<double> >p, std::vector<int> ki, int J, boost::numeric::ublas::vector<double>& w, std::list<std::list<boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> > > Uj, std::list<std::vector<int> > part, std::vector<int> G, double NRMSE, double &R, double mx, double mn, double mxg, double mng, double alpha) {
// x: is the input data vector represented on the mesh in p (input)
// p: is a matrix containing the mesh points coordinates (input)
// ki: is a vector containing the wavelet orders in all directions (input)
// J: required level of detail, if J=-1 then full detail is considered (input-output)
// w: transformed vector (output)
// Uj, part and G contain the moment matrices and mesh hierarchy (inputs)
// NRMSE: required accuracy (input)
// R: Achieved compression ratio (output)
// mx, mn: maximum and minimum values of the local data x (inputs)
// mxg, mng: maximum and minimum values of the global data x (inputs)
// alpha: threshold coefficient (input)

  int NR;
  int N=x.size();
  
  // the universal thresholding technique
  double thresh0 = alpha/sqrt((mxg-mng)/(mx-mn))/0.6745e0*sqrt(2.0e0*log(double(N)));
  double thresh, ew;
  boost::numeric::ublas::vector<double> xw(N);
    
  NRMSE*=(mxg-mng);
  Perform_Alpert_transform(x, ki, J,1, w, Uj, part, G);

  double s=0.0e0;
  double m;
  for (int i = 0; i < N; i++) {
    xw(i)=fabs(w(i))/(mxg-mng);
    s+=xw(i);
  }
  m=s/double(N);
  
  s=0.0e0;
  for (int i = 0; i < N; i++) 
    s+=fabs(xw(i)-m);
  s/=double(N);
  thresh0*=s/m;
  
//   std::cout << thresh0 << std::endl;
  
  // if provided accuracy is negative then do normal hard-thresholding
  if (NRMSE<=0.0e0) {    
    NR=N;
    for (int i = 0; i < N; i++) {
      if (xw(i)<thresh0) {
	// coefficients of magnitudes less than the threshold are set to zero
	w(i)=0.0e0;
	NR-=1;
      }
    }
    
    R = double(N)/double(NR);
    return;
  }
  // if provided accuracy is positive then choose a threshold that satisfies it
  //// this is the recommended operation mode for the time being ////
  else {
    ew=0.0e0;
    thresh=thresh0/10000.0e0;
    
    // start with all coefficients in w
    // iterate as long the error is lower than the required NRMSE
    while ((ew<=NRMSE)&&(thresh<=thresh0)) {
      // increase the threshold by 20%, this part still needs some work
      thresh*=1.2e0;
      
      NR=N;
      for (int i = 0; i < N; i++) {
	if (std::abs(w(i))<thresh) {
	  w(i)=0.0e0;
	  NR-=1;
	}
      }      
      R = double(N)/double(NR);
      
      Perform_Alpert_transform(w, ki, J,-1, xw, Uj, part, G);
      
      double ew=0.0e0;
      for (int i = 0; i < N; i++) 
	ew+=(x(i)-xw(i))*(x(i)-xw(i));
      ew=std::sqrt(ew/double(N));
    }
    return;
  }
  
}

// Performs forwards and inverse wavelet transform 
void Perform_Alpert_transform(boost::numeric::ublas::vector<double> f, std::vector<int> ki, int J,int dir, boost::numeric::ublas::vector<double>& w, std::list<std::list<boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> > > Uj, std::list<std::vector<int> > part, std::vector<int> G) {
// f: is the input data vector represented on the mesh in p (input)
// ki: is a vector containing the wavelet orders in all directions (input)
// w: transformed vector (output)
// Uj, part and G contain the moment matrices and mesh hierarchy (inputs)

  std::vector<int> cG;
  
  int k2=1;
  for (int i = 0; i < ki.size(); i++)
        k2*=ki[i];
  
  int P=G.size();
  int nj, mj, nnj, offs;
  int N = f.size();
  
  cG=G;
  partial_sum(G.begin(), G.end(), cG.begin());
  std::vector<int> si(G.size()+1);
  si[0]=1;
  for (int i = 0; i < cG.size(); i++)
    si[i+1]=cG[i]+1;
  
  std::vector<int> seli, selj;
  std::list<std::list<boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> > >::iterator itUj;
  
  std::vector<int> Jl(J);
  itUj=Uj.begin();
  
  if (dir>=1) {    
    for (int j = 0; j < J; j++)
      Jl[j]=j+1;
  }
  else {
    for (int j = 0; j < Uj.size()-1; j++)
      ++itUj;
    for (int j = J-1; j >=0; j--)
      Jl[J-j-1]=j+1;
  }
    
  w = f;
  boost::numeric::ublas::vector<double> r, rs;
  
  std::list<std::vector<int> >::iterator itp;
  std::list<boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> > Mi,Ui,Uti;
  std::list<boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> >::iterator itM, itU;
  boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> Uii, Mii, Vii;
  
  for (int j = 0; j < J; j++) {   // for each scale   
       
    Ui=*itUj;
        
    if (Jl[j]==1) {
     
      // Special treatment for the first level
      r = w;
      itp=part.begin();
      itU=Ui.begin();
      
      for (int i = 0; i < P; i++)  {
	selj = *itp;
	
	boost::numeric::ublas::vector<double> ws(std::max(seli.size(),selj.size()));
	
	offs = si[i]-i*k2;   // offset on row
        int lo = G[i]-k2;          // length on row
        
        seli.clear();
	for (int l = 0; l < lo; l++)
	  seli.push_back(l+offs);
	
        Uii=*itU;
	int nc=Uii.size2();
	int nc1=Uii.size1();
	
	if (dir>=1) {
	  for (int l = 0; l < selj.size(); l++)
            ws(l)=w(selj[l]);
	  rs=boost::numeric::ublas::prod(boost::numeric::ublas::project(Uii,boost::numeric::ublas::range(k2,nc1),boost::numeric::ublas::range(0,nc)),ws);
	  for (int l = 0; l < seli.size(); l++)
	    r(seli[l])=rs(l);
	}
	else {	  
	  for (int l = 0; l < seli.size(); l++)
	    ws(l)=w(seli[l]);
	  rs=boost::numeric::ublas::prod(boost::numeric::ublas::trans(boost::numeric::ublas::project(Uii,boost::numeric::ublas::range(k2,nc1),boost::numeric::ublas::range(0,nc))),ws);
	  for (int l = 0; l < selj.size(); l++)
	    r(selj[l])=rs(l);
	}
	 
	offs = N - P*k2 + 1 + i*k2;
        ws.clear();
	
	if (offs>=0) {
	  seli.clear();
	  for (int l = 0; l < k2; l++)
	    seli.push_back(l+offs-1);
	  
	  if (dir>=1) { 
	    for (int l = 0; l < selj.size(); l++)
	      ws(l)=w(selj[l]);
	    rs=boost::numeric::ublas::prod(boost::numeric::ublas::project(Uii,boost::numeric::ublas::range(0,k2),boost::numeric::ublas::range(0,nc)),ws);
	    for (int l = 0; l < seli.size(); l++)
	      r(seli[l])=rs(l);  
	  }
	  else {
	    for (int l = 0; l < seli.size(); l++)
	      ws(l)=w(seli[l]);
	    rs=boost::numeric::ublas::prod(boost::numeric::ublas::trans(boost::numeric::ublas::project(Uii,boost::numeric::ublas::range(0,k2),boost::numeric::ublas::range(0,nc))),ws);
	    for (int l = 0; l < selj.size(); l++)
	      r(selj[l])+=rs(l);
	  }
	}
	else
	  std::cout << "Pbm empty bin." << std::endl;
	
        ++itU;
	++itp;
      }
      
      w = r;
      
    }
    
    else {
      
      nnj=1;  
      for (int k = 1; k <= (Jl[j]-1); k++)  nnj*=2; // n/(2^j*k);
      nj=floor(double(P)/double(nnj));
      mj = floor(double(P)/double(nnj)*double(2*k2));     // total length of the blocks
      
      r = w;
      offs = N-mj;
      
      itU=Ui.begin();
      for (int i = 0; i <nj; i++) {
	seli.clear();
	selj.clear();
	
	for (int l = 1; l <= k2; l++)
	  seli.push_back(offs+k2*i+l-1);
	for (int l = 1; l <= 2*k2; l++)
	  selj.push_back(offs+2*k2*i+l-1);
	
	boost::numeric::ublas::vector<double> ws(selj.size());

	Uii=*itU;
	int nc=Uii.size2();
	int nc1=Uii.size1();	
	
	if (dir>=1) { 
	    for (int l = 0; l < selj.size(); l++)
	      ws(l)=w(selj[l]);
	    rs=boost::numeric::ublas::prod(boost::numeric::ublas::project(Uii,boost::numeric::ublas::range(nc1/2,nc1),boost::numeric::ublas::range(0,nc)),ws);
	    for (int l = 0; l < seli.size(); l++)
	      r(seli[l])=rs(l);  
	  }
	else {
	    for (int l = 0; l < seli.size(); l++)
	      ws(l)=w(seli[l]);
	    rs=boost::numeric::ublas::prod(boost::numeric::ublas::trans(boost::numeric::ublas::project(Uii,boost::numeric::ublas::range(nc1/2,nc1),boost::numeric::ublas::range(0,nc))),ws);
	    for (int l = 0; l < selj.size(); l++)
	      r(selj[l])=rs(l);
	}
	
	ws.clear();
	seli.clear();
	for (int l = 1; l <= k2; l++)
	  seli.push_back(offs+mj/2+k2*i+l-1);
	
	if (dir>=1) { 
	    for (int l = 0; l < selj.size(); l++)
	      ws(l)=w(selj[l]);
	    rs=boost::numeric::ublas::prod(boost::numeric::ublas::project(Uii,boost::numeric::ublas::range(0,nc1/2),boost::numeric::ublas::range(0,nc)),ws);
	    for (int l = 0; l < seli.size(); l++)
	      r(seli[l])=rs(l);  
	  }
	else {
	    for (int l = 0; l < seli.size(); l++)
	      ws(l)=w(seli[l]);
	    rs=boost::numeric::ublas::prod(boost::numeric::ublas::trans(boost::numeric::ublas::project(Uii,boost::numeric::ublas::range(0,nc1/2),boost::numeric::ublas::range(0,nc))),ws);
	    for (int l = 0; l < selj.size(); l++)
	      r(selj[l])+=rs(l);
	}
	
	++itU;  
      }
      w=r;
    }
    
    if (dir>=1)
      ++itUj;
    else
      --itUj;
      
  }
}

// Computes the moment matrices and mesh hirearchy
void Compute_Uj(std::vector<std::vector<double> >p, std::vector<int> ki, int &J, std::list<std::list<boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> > > &Uj, std::list<std::vector<int> > &part, std::vector<int> &G) {
// p: is a matrix containing the mesh points coordinates (input)
// ki: is a vector containing the wavelet orders in all directions (input)
// J: required level of detail, if J=-1 then full detail is considered (input-output)
// Uj, part and G contain the moment matrices and mesh hierarchy (outputs)

  int d=(p[0]).size();
  if (d!=ki.size()) {
    std::cout << "Error: The size of vanishing moments array should be identical to the number of dimensions" << std::endl;
    abort();
  }
  
  std::vector<int> seli;
  int k2=1;
  for (int i = 0; i < ki.size(); i++)
        k2*=ki[i];

  dichotomic_grouping(p, k2, part, G);
    
  int P=G.size();
  
  // The required wavelet level is now supplied as an input
  // If we supply a negative level, the code computes and uses the maximum number of levels
  if (J<=0)
    J=int(log2(double(P))+1.0);

  

  std::list<std::vector<int> >::iterator itp;
  itp=part.begin();
  
  std::list<boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> > Mi,Ui,Uti;
  std::list<boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> >::iterator itM, itU;

  // building a table for the vanishing moments orders in each direction
  std::vector< std::vector < int > > table;
  
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
  
  Uj.push_back(Ui);
  
  boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> Uii, Mii, Vii;
  int nj, nnj;
    
  for (int j = 2; j <= J; j++) {   // for each scale
    
    // at this scale, we have nj = P/2^(j-1) groups
    nnj=1;  
    for (int k = 1; k <= (j-1); k++)  nnj*=2; // n/(2^j*k);
    nj=floor(double(P)/double(nnj));
    
    // update each sub matrix
    itU=Ui.begin();
    itM=Mi.begin();
    std::list<boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> > MMi,UUi,UUti;
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
    Uj.push_back(Ui);
  }
}
