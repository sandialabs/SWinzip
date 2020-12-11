#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include "alpert_transform.h"

void alpert_transform(double *(*Uj)[], double *f, double *r, int dire, int *(*part)[], int *G, int P, int J, int k2, int N) {
   int i,k,j,m,kk,twok2,nJj,nJ;
   int *si, *j_list, *seli, *selj;
   int *selio, *seljo;
   double *(*Ui)[];
   double *w;
   int offs,lo,mj,nj;
   double *U, *Md;
   double *Mdo, *Uo;
   
   twok2=2*k2;
   w = (double *) malloc(sizeof(double)*N);
   for (i=0; i<N; i++)
     w[i]=f[i];
     
   si = (int *) malloc(sizeof(int)*(P+1));
   si[0]=0;
   for (i=0; i<P; i++)
     si[i+1]=G[i]+si[i];
      
   j_list = (int *) malloc(sizeof(int)*J);
   for (j=0; j<J; j++)
     j_list[j]=j+1;
   
   if (dire==-1)
     for (j=0; j<J; j++)
       j_list[j]=J-j;
   
   selio = (int *) malloc(sizeof(int)*k2);
   seljo = (int *) malloc(sizeof(int)*twok2);
   Mdo = (double *) malloc(sizeof(double)*k2*twok2);
   
   nJ=1;
   for (i=0; i<J; i++)
     nJ*=2;
   nJ-=1;
   
   for (j = 0; j < J; j++) {
     
     nJj=pow(2,J-j_list[j]);
     
     Ui=(double*(*)[])malloc(nJj*sizeof(double*));
     for (k=0; k<nJj; k++)
       (*Ui)[k]=(*Uj)[nJ+1-2*nJj + k];
     
     if (j_list[j]==1) {
       
       for (k=0; k<N; k++)
	 r[k]=w[k];
       
       for (i=0; i<P; i++) {
	 
	 selj=(*part)[i];
	 offs = si[i]-i*k2;
         lo = G[i]-k2;
	 
	 seli = (int *) malloc(sizeof(int)*lo);
	 for (k=0; k<lo; k++)
	   seli[k]=k+offs;
	 
	 kk=G[i];
	 U = (*Ui)[i];
	 Md = (double *) malloc(sizeof(double)*lo*kk);
	 
	 for (k=0; k<kk; k++)
	   for (m=k2; m<kk; m++)
	     Md[(kk-k2)*k+m-k2]=U[kk*k+m];
	   
	 if (dire==1) {
	   for (k=0; k<lo; k++) {
	     r[seli[k]]=0.0e0;
	     for (m=0; m<kk; m++)
		r[seli[k]]+=Md[m*lo+k]*w[selj[m]];
	   }
	 }
	 else {
	   for (k=0; k<kk; k++) {
	     r[selj[k]]=0.0e0;
	     for (m=0; m<lo; m++)
		r[selj[k]]+=Md[k*lo+m]*w[seli[m]];
	   }
	 }
	   
	 offs = N - P*k2 + 1 + i*k2;
	 
	 free(Md);
	 
	 if (offs>=0) {
	   for (k=0; k<k2; k++)
	     selio[k]=k+offs-1;
	      
	   Md = (double *) malloc(sizeof(double)*k2*kk);	 
	   for (k=0; k<kk; k++)
	     for (m=0; m<k2; m++)
	       Md[k2*k+m]=U[kk*k+m];
	     
	   if (dire==1) {
	      for (k=0; k<k2; k++) {
		r[selio[k]]=0.0e0;
		for (m=0; m<kk; m++)
		   r[selio[k]]+=Md[m*k2+k]*w[selj[m]];
	      }
	    }
	    else
	      for (k=0; k<kk; k++)
		for (m=0; m<k2; m++)
		  r[selj[k]]+=Md[k*k2+m]*w[selio[m]];
		
	    free(Md);   
	 }
	 else
	   printf("Pbm empty bin.\n");
         
	 free(seli);
       }
       
       for (k=0; k<N; k++)
	 w[k]=r[k];
              
     }
     
     else {       
       
       nj=P/pow(2,j_list[j]-1);
       mj=nj*2*k2;
       
       for (k=0; k<N; k++)
	 r[k]=w[k];
       
       offs=N-mj;
       
       for (i=0; i<nj; i++) {
	 
	 for (k=1; k<=twok2; k++)
	   seljo[k-1]=offs+twok2*i+k-1;
	 for (k=1; k<=k2; k++)
	   selio[k-1]=offs+k2*i+k-1;
	 
	 Uo = (*Ui)[i];
	  
	 for (k=0; k<twok2; k++)
	   for (m=0; m<k2; m++)
	     Mdo[k2*k+m]=Uo[twok2*k+m+k2];
	 
	 if (dire==1) {     
	   for (k=0; k<k2; k++) {
	     r[selio[k]]=0.0e0;
	     for (m=0; m<twok2; m++)
	       r[selio[k]]+=Mdo[m*k2+k]*w[seljo[m]];
	    }     
	 }
	 else {
	   for (k=0; k<twok2; k++) {
	     r[seljo[k]]=0.0e0;
	     for (m=0; m<k2; m++)
	       r[seljo[k]]+=Mdo[k*k2+m]*w[selio[m]];	     
	    }
	 }
       
       for (k=1; k<=k2; k++)
	 selio[k-1]=offs+mj/2+k2*i+k-1;
       
       for (k=0; k<twok2; k++)
	 for (m=0; m<k2; m++)
	   Mdo[k2*k+m]=Uo[twok2*k+m];
       
       if (dire==1) {
	 for (k=0; k<k2; k++) {
	   r[selio[k]]=0.0e0;
	   for (m=0; m<twok2; m++)
	     r[selio[k]]+=Mdo[m*k2+k]*w[seljo[m]];
	 }
       }
       else
	 for (k=0; k<twok2; k++) 
	   for (m=0; m<k2; m++)
	     r[seljo[k]]+=Mdo[k*k2+m]*w[selio[m]];
	    
     }
     
     for (k=0; k<N; k++)
       w[k]=r[k];
     
     }
     
     free(Ui);
   }
   
   free(w);
   free(si);free(j_list);
   free(selio);free(seljo);
   free(Mdo);
}

void thresh(double *w, int N, double th0, double mxmn, double mxmng, double *th) {
  double c,d;
  
  fmean(w,N,mxmng,&d);
  fmeana(w,N,d,mxmng,&c);
  *th=th0*sqrt(mxmng/mxmn)/d*c/0.6745e0*sqrt(2.0e0*log((double)N));
}

void get_large_indices(double *w, int N, double th, double mxmn, int *M) {
  int i;
  int n=0;
  
  for (i=0; i<N; i++)
    if (fabs(w[i])/mxmn >= th)
      n+=1;
    
  *M=n;  
 }
 