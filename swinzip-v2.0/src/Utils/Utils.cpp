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

#include "Utils.hpp"

using namespace std;


void dichotomic_grouping(vector<vector<double> >& p, int ptmax, list<vector<int> > &part, vector<int>& G) {
  
  int cpt=0;
  int depthmax=1000000;
  int pp, type, nP;
  int N=p.size(), d=(p[0]).size();
  type=d;
  
  list<vector<int> > part_prev, part1;
  list<vector<int> >::iterator it1, itp;
  
  vector<int> Ni, P, IP;
  for (int i = 0; i < N; i++)
    Ni.push_back(i);
  
  part.push_back(Ni);
  int done=0;
  
  while ((cpt<depthmax)&&(done==0)) {
    
    cpt+=1;
    part_prev = part;
    part1 = part;
    part.clear();
    
    it1 = part1.begin();
    
    for (int k=0;k<part1.size();k++) {
      
      P = *it1;
      nP=P.size();      
      pp = floor(nP/2);    // half # points

      int s = (cpt-1)%type;
      vector<double> tmp;
      vector<size_t> I(nP);
      for (int j=0;j<nP;j++) {
	tmp.push_back(p[P[j]][s]);
	I[j]=j;
      }
      sortandindices( tmp, I);      
      vector<int> ppart1(&I[0], &I[pp]);
      vector<int> ppart2(&I[pp], &I[nP]);

      // assign new partition
      IP.clear();
      for (int j=0;j<ppart1.size();j++)
	IP.push_back(P[ppart1[j]]);
      part.push_back(IP);
            
      IP.clear();
      for (int j=0;j<ppart2.size();j++)
	IP.push_back(P[ppart2[j]]);
      part.push_back(IP);
            
      ++it1;
    }
    
    itp=part.begin();
    for (int j=0;j<part.size();j++) {
      if ((*itp).size()<ptmax) {
	part = part_prev;
	done = 1;
	break;
      }
    ++itp;  
    }
  
  }

  itp=part.begin();
  for (int i=0;i<part.size();i++) {
    G.push_back((*itp).size());
    ++itp;
  }

}


void sortandindices(vector<double>& numbers, vector<size_t>& indexes) {

  sort(
      indexes.begin(),
      indexes.end(),
      [&numbers](size_t i1, size_t i2) {return numbers[i1] < numbers[i2];});
  
}

void prctile(double *f, int N, double p, double *c) {
  int i;
  vector<double> tmp(N);
  vector<size_t> I(N);
  
  for (i=0;i<N;i++) {
    tmp[i]=f[i];
    I[i]=i;
  }
  
  sortandindices(tmp, I);
  
  i=ceil(p*(double)N/100.0e0);
  if (i<0)
    i=0;
  if (i>=N)
    i=N-1;
  *c=tmp[I[i]];
}