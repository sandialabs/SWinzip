#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include "alpert_utils.h"
#include <stdlib.h>

void get_Nbins(int *ki, int N, int d,int *J, int *k2, int *P, int *nJ) {
   int ptmax,i,j;
   long int depthmax;
   
   depthmax=1000000;
   
   ptmax=1;
   for (i=0; i<d; i++)
     ptmax*=ki[i];
   
   /* Estimating the number of bins */
   for (i=0; i<depthmax; i++) {
     if (pow(2,i) >= N/ptmax) {
       *P=pow(2,i-1);
       *J=i;
       break;
     }
   }
   
   *k2=ptmax;
   
   j=1;
   for (i=0; i<*J; i++)
    j*=2;
   j-=1;
   
   *nJ=j;

}

int compqi (const void *a, const void *b) {
    struct pos *f, *s;
    f = (struct pos *) a;
    s = (struct pos *) b;
    if (f->posx > s->posx) return  1;
    if (f->posx < s->posx) return -1;
    return 0;
}

void dichotomic_grouping(double *p, int N, int d, int *G, int *(*part)[], int J, int k2, int P) {
   
   int i,j,k,cpt,done,lp,nP,pp,s;
   int *Ni, *nPi, *nPi_prev, *I, *PP;
   int *(*part_prev)[];
   int *ppart1, *ppart2;
   struct pos *tmp;
   
   cpt=0;   
      
   Ni = (int *) malloc(sizeof(int)*N);
   for (i=0; i<N; i++)
     Ni[i]=i;
   
   done=0;
   nPi = (int *) malloc(sizeof(int)*P);
   nPi_prev = (int *) malloc(sizeof(int)*P);
   
   part_prev=(int*(*)[])malloc(P*sizeof(int*));
   
   for (i=0; i<P; i++) {
    nPi[i]=0;
    nPi_prev[i]=0;
    G[i]=0;
    (*part)[i]=NULL;
    (*part_prev)[i]=NULL;
   }
   
   (*part)[0] = (int *) malloc(sizeof(int)*N);
   memcpy((*part)[0], Ni, sizeof(int)*N);
   free(Ni);
   
   
   lp=0;
   nP=N;
   nPi[0]=N;

   for (done=0; done<J-1; done++) {
     cpt+=1;
     s = (cpt-1)%d;
     
     for (i=0; i<pow(2,lp); i++) {
       (*part_prev)[i] = (int *) malloc(sizeof(int)*nPi[i]);
       memcpy((*part_prev)[i], (*part)[i], sizeof(int)*nPi[i]);
       free((*part)[i]);
       nPi_prev[i]=nPi[i];
     }
     
     for (k=0; k<pow(2,lp); k++) {       
       PP=(*part_prev)[k];
       nP=nPi_prev[k];

       pp = floor(nP/2);       
       I = (int *) malloc(sizeof(int)*nP);
       tmp = (struct pos *) malloc(sizeof(struct pos)*nP);
       
       for (j=0; j<nP; j++) {
	 tmp[j].posx=p[d*PP[j]+s];
	 tmp[j].posId=j;
       }
       
       /* sorting tmp */
       /*mergeSort(tmp, 0, nP-1);*/
       qsort(tmp,nP,sizeof(struct pos),compqi);
       for (j=0; j<nP; j++)
	 I[j]=tmp[j].posId;
       
       ppart1 = (int *) malloc(sizeof(int)*pp);
       ppart2 = (int *) malloc(sizeof(int)*(nP-pp));
       
       for (j=0; j<pp; j++)
	 ppart1[j]=PP[I[j]];
       
       for (j=pp; j<nP; j++)
	 ppart2[j-pp]=PP[I[j]];
       
       (*part)[2*k] = (int *) malloc(sizeof(int)*pp);
       (*part)[2*k+1] = (int *) malloc(sizeof(int)*(nP-pp));
       memcpy((*part)[2*k], ppart1, sizeof(int)*pp);
       memcpy((*part)[2*k+1], ppart2, sizeof(int)*(nP-pp));
       
       nPi[2*k]=pp;
       nPi[2*k+1]=nP-pp;
       
       free(I);free(tmp);
       free(ppart1);free(ppart2);
     }

     lp+=1;
     
     for (i=0; i<pow(2,lp-1); i++)
       free((*part_prev)[i]);        
   }
   
   for (i=0; i<P; i++)
     G[i]=nPi[i];
   
   free(nPi);free(nPi_prev);
   free(part_prev);
}


void compute_Uj(double *p, int *(*part)[], int *ki, int *G, int k2, int P, int J, int d, double *(*Uj)[]) {
   int i,j,k,m,cnt,twok2,kk,ii,nj,info,lwork;
   int *table, *iid, *seli;
   char trans='N';
   double alpha=1.0e0, beta=0.0e0;
   double *M, *Ut, *Mu, *Md, *tau, *work;
   double *(*Mi)[], *(*Ui)[];
   double *(*MMi)[], *(*UUi)[];
   double *Mo, *MMo;
   int nJ;
   
   twok2=2*k2;
   table = (int *) malloc(sizeof(int)*k2*d);
   iid = (int *) malloc(sizeof(int)*d);
   cnt=0;
   
   for (j=0; j<d; j++)
     iid[j]=0;
   
   while (cnt<k2) {
    for (j=0; j<d; j++)
      table[d*cnt+j] = iid[j];
    
    iid[0]+=1;
    i = 0;
    while ((i<(d-1)) && (iid[i]==ki[i])) {
      iid[i]=0;
      i+=1;
      iid[i]+=1;
    }
    
    cnt+=1;
   }
   
   Mi=(double*(*)[])malloc(P*sizeof(double*));
   Ui=(double*(*)[])malloc(P*sizeof(double*));
   Mo = (double *) malloc(sizeof(double)*twok2*twok2);
   MMo = (double *) malloc(sizeof(double)*twok2*twok2);
   for (i=0; i<P; i++) {
     seli=(*part)[i];
     kk=G[i];
     
     lwork=kk;
     M = (double *) malloc(sizeof(double)*kk*twok2);
     Ut = (double *) malloc(sizeof(double)*kk*kk);
     tau = (double *) malloc(sizeof(double)*kk);
     work = (double *) malloc(sizeof(double)*lwork);
     
     for (j=0; j<twok2; j++)
       for (ii=0; ii<kk; ii++)
	 M[j*kk+ii]=0.0e0;
       
     for (j=0; j<k2; j++) {
        for (ii=0; ii<kk; ii++) {
        M[j*kk+ii]=1.0e0;
        for (k=0; k<d; k++) 
          M[j*kk+ii]*=pow(p[d*seli[ii]+k],(double)table[j*d+d-k-1]);
      }
     }
     
     (*Mi)[i] = (double *) malloc(sizeof(double)*kk*twok2);
     memcpy((*Mi)[i], M, sizeof(double)*kk*twok2);
     
     for (j=0;j<kk*kk;j++)
       Ut[j]=M[j];
     
     dgeqrf_(&kk, &kk, &Ut[0], &kk, &tau[0], &work[0], &lwork, &info);
     dorgqr_(&kk, &kk, &kk, &Ut[0], &kk, &tau[0], &work[0], &lwork, &info);
     transpose(&Ut[0], kk);
     
     (*Ui)[i] = (double *) malloc(sizeof(double)*kk*kk);
     memcpy((*Ui)[i], Ut, sizeof(double)*kk*kk);
     
     (*Uj)[i] = (double *) malloc(sizeof(double)*kk*kk);
     memcpy((*Uj)[i], Ut, sizeof(double)*kk*kk);
     
     free(tau);free(work);
     free(M);free(Ut);
     
   }
   
   nJ=P;
   
   lwork=twok2;
   
   Md = (double *) malloc(sizeof(double)*twok2*k2);
   Mu = (double *) malloc(sizeof(double)*twok2*k2);
   tau = (double *) malloc(sizeof(double)*twok2);
   work = (double *) malloc(sizeof(double)*lwork);
       
   for (j=2; j<=J; j++) {
     nj = P/pow(2,j-1);
     
     MMi=(double*(*)[])malloc(nj*sizeof(double*));
     UUi=(double*(*)[])malloc(nj*sizeof(double*));
     
     for (i=0; i<nj; i++) {
 
       for (k=0; k<twok2*k2; k++) {
	 Mu[k]=0.0e0;
	 Md[k]=0.0e0;
       }
       
       if (j==2)
	 kk=G[2*i];
       else
	 kk=twok2;
       
       Ut = (double *) malloc(sizeof(double)*kk*k2);
       for (k=0; k<kk; k++)
	 for (m=0; m<k2; m++)
	   Ut[k2*k+m]=(*Ui)[2*i][kk*k+m];   
       dgemm_(&trans, &trans, &k2, &twok2, &kk, &alpha, &Ut[0], &k2, (*Mi)[2*i], &kk, &beta, &Mu[0], &k2);
	 
       if (j==2)
	 kk=G[2*i+1];
       else
	 kk=twok2;
       
       free(Ut);
       Ut = (double *) malloc(sizeof(double)*kk*k2);
       for (k=0; k<kk; k++)
	 for (m=0; m<k2; m++)
	   Ut[k2*k+m]=(*Ui)[2*i+1][kk*k+m];
       dgemm_(&trans, &trans, &k2, &twok2, &kk, &alpha, &Ut[0], &k2, (*Mi)[2*i+1], &kk, &beta, &Md[0], &k2);
       
       for (k=0; k<twok2; k++) {
	 for (m=0; m<k2; m++) {
	   Mo[k*twok2+m]=Mu[k*k2+m];
	   Mo[k*twok2+m+k2]=Md[k*k2+m];
	   MMo[k*twok2+m]=Mu[k*k2+m];
	   MMo[k*twok2+m+k2]=Md[k*k2+m];
	 }
       }
       
       (*MMi)[i] = (double *) malloc(sizeof(double)*twok2*twok2);
       memcpy((*MMi)[i], MMo, sizeof(double)*twok2*twok2);
       
       dgeqrf_(&twok2, &twok2, &Mo[0], &twok2, &tau[0], &work[0], &lwork, &info);
       dorgqr_(&twok2, &twok2, &twok2, &Mo[0], &twok2, &tau[0], &work[0], &lwork, &info);
       transpose(&Mo[0], twok2);
       
       (*UUi)[i] = (double *) malloc(sizeof(double)*twok2*twok2);
       memcpy((*UUi)[i], Mo, sizeof(double)*twok2*twok2);       
       free(Ut);
    
     }
     
     for (i=0; i<2*nj; i++) {
       free((*Mi)[i]);
       free((*Ui)[i]);
     }
     free(Mi);free(Ui);
     
     Mi=(double*(*)[])malloc(nj*sizeof(double*));
     Ui=(double*(*)[])malloc(nj*sizeof(double*));
     
     for (i=0; i<nj; i++) {
       (*Mi)[i] = (double *) malloc(sizeof(double)*twok2*twok2);
       memcpy((*Mi)[i], (*MMi)[i], sizeof(double)*twok2*twok2);
       free((*MMi)[i]);
       
       (*Ui)[i] = (double *) malloc(sizeof(double)*twok2*twok2);
       memcpy((*Ui)[i], (*UUi)[i], sizeof(double)*twok2*twok2);
       
       (*Uj)[nJ+i] = (double *) malloc(sizeof(double)*twok2*twok2);
       memcpy((*Uj)[nJ+i], (*UUi)[i], sizeof(double)*twok2*twok2);
       
       free((*UUi)[i]);
     }
     
     nJ+=nj;
     
     free(MMi);
     free(UUi);
     
   }
   
   if (J>=2) {
     for (i=0; i<nj; i++) {
	free((*Mi)[i]);
	free((*Ui)[i]);
     }
   }
   free(Mi);free(Ui);
   
   free(Mu);free(Md);
   free(tau);free(work);
   free(table);free(iid);
   free(Mo);free(MMo);
}
