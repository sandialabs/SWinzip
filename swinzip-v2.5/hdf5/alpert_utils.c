#include <stdio.h>
#include "alpert_utils.h"
#include <malloc.h>
#include <stdlib.h>
#include <math.h>

void transpose(double *M, int N) {
    int i,j;
    double t;
    for (i=0; i<N; i++) {
      for (j=0; j<i; j++) {
	t=M[i*N+j];
	M[i*N+j]=M[j*N+i];
	M[j*N+i]=t;
      }
    }
}

int compq (const void *a, const void *b) {
    struct pos *f, *s;
    f = (struct pos *) a;
    s = (struct pos *) b;
    if (f->posx > s->posx) return  1;
    if (f->posx < s->posx) return -1;
    return 0;
}

void prctile(double *f, int N, double p, double *c) {
  int i;
  struct pos *v;
  
  v=(struct pos *) malloc(sizeof(struct pos)*N);    
  for (i=0; i<N; i++) {
    v[i].posx=f[i];
    v[i].posId=i;
  }
  
  qsort(v,N,sizeof(struct pos),compq);
  
  i=ceil(p*(double)N/100.0e0);
  if (i<0)
    i=0;
  if (i>=N)
    i=N-1;
  *c=v[i].posx;
  
  free(v);
}

void fmean(double *f, int N, double mxmn, double *m) {
  int i;
  double d=0.0e0;
  
  for (i=0; i<N; i++) {
    d+=fabs(f[i]);
  }
  
  d/=(double)N;
  d/=mxmn;
  *m=d;
}

void fmeana(double *f, int N, double x, double mxmn, double *m) {
  int i;
  double d=0.0e0;
  
  for (i=0; i<N; i++) {
    d+=fabs(fabs(f[i])/mxmn-x);
  }
  
  d/=(double)N;
  *m=d;
}

void mxmn(double *f, int N, double *c) {
  int i;
  double mn,mx;
  
  mn=f[0];
  for (i=1; i<N; i++)
    if (f[i]<=mn)
      mn=f[i];
  
  mx=f[0];
  for (i=1; i<N; i++)
    if (f[i]>=mx)
      mx=f[i];
    
  *c=mx-mn;  
}
