/*
 * This file is an example of an HDF5 filter plugin.
 * The plugin can be used with the HDF5 library vesrion 1.8.11+ to read
 * HDF5 datasets compressed with ALPERT.
 * Written by Maher Salloum, October 2020
 */

#include <sys/types.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#if defined(_WIN32)
#include <Winsock2.h>
#endif
#include <H5PLextern.h>

#include "alpert_transform.h"
#include "init_alpert.h"
#include "alpert_utils.h"

static size_t H5Z_filter_ALPERT(unsigned int flags, size_t cd_nelmts,
        const unsigned int cd_values[], size_t nbytes,
        size_t *buf_size, void **buf);

/*#define H5Z_FILTER_ALPERT 32004*/

int *pp;
double *p;
int *ki, *G;
int *(*part)[];
double *(*Uj)[];
int N, kki, d, P,J,k2,nJ,Ng,nc;
int Nx,Ny,Nz;
int Ngx,Ngy,Ngz;

typedef struct {
    double w;
    int ii;
} doubleint;

const H5Z_class2_t H5Z_ALPERT[1] = {{
        H5Z_CLASS_T_VERS,       /* H5Z_class_t version */
        (H5Z_filter_t)H5Z_FILTER_ALPERT,         /* Filter id number             */
        1,              /* encoder_present flag (set to true) */
        1,              /* decoder_present flag (set to true) */
        "HDF5 ALPERT filter; see http://www.hdfgroup.org/services/contributions.html",
        /* Filter name for debugging    */
        NULL,                       /* The "can apply" callback     */
        NULL,                       /* The "set local" callback     */
        (H5Z_func_t)H5Z_filter_ALPERT,         /* The actual filter function   */
}};

/*H5PL_type_t   H5PLget_plugin_type(void) {return H5PL_TYPE_FILTER;}
const void *H5PLget_plugin_info(void) {return H5Z_ALPERT;}*/

static size_t H5Z_filter_ALPERT(unsigned int flags, size_t cd_nelmts,
        const unsigned int cd_values[], size_t nbytes,
        size_t *buf_size, void **buf)
{
    doubleint * outBufc = NULL;
    double * outBufd = NULL;
    size_t ret_value;
    
    int i,ii,M,n,j,k,l;
    double *w, *f;
    double th, mx, mn, mxmn, mxmng,th0;
    double x1,x2;
    
    ii=2;
    if (pp == NULL) {
      pp=&ii;
            
      d = cd_values[1];
      kki=cd_values[0];
      ki = (int *) malloc(sizeof(int)*d);
      for (i=0; i<d; i++)
	ki[i]=kki;
      
      Nx=cd_values[2];
      Ny=cd_values[3];
      Nz=cd_values[4];      
      N=Nx*Ny*Nz;
      
      Ngx=cd_values[11];
      Ngy=cd_values[12];
      Ngz=cd_values[13];
      
      p = (double *) malloc(sizeof(double)*d*N);
      if (d==1)
	for (i=0;i<N;i++)
	  p[i] = (double)i/(double)(Ngx-1);
      else if (d==2) {
	k=0;
	for (i=0;i<Nx;i++) {
	  for (j=0;j<Ny;j++) {
	    p[2*k]=(double)i/(double)(Ngx-1);
	    p[2*k+1]=(double)j/(double)(Ngy-1);
	    k+=1;
	  }
	}
      }
      else {
	l=0;
	for (i=0;i<Nx;i++) {
	  for (j=0;j<Ny;j++) {
	    for (k=0;k<Nz;k++) {
	      p[3*l]=(double)i/(double)(Ngx-1);
	      p[3*l+1]=(double)j/(double)(Ngy-1);
	      p[3*l+2]=(double)k/(double)(Ngz-1);
	      l+=1;
	    }
	  }
	}
      }
      
      get_Nbins(ki, N, d, &J, &k2, &P, &nJ);
      
      G = (int *) malloc(sizeof(int)*P);
      part=(int*(*)[])malloc(P*sizeof(int*));
      dichotomic_grouping(p, N, d, G, part, J, k2, P);
	
      Uj=(double*(*)[])malloc(nJ*sizeof(double*));
      compute_Uj(p, part, ki, G, k2, P, J, d, Uj);
      
      Ng=Ngx*Ngy*Ngz/N;
      nc=0;
    }
    
    x1=(double)cd_values[6];
    x2=pow(10.0e0,(double)cd_values[5]);
    mxmng = (cd_values[7]==1) ? x1*x2 : x1/x2;
    
    x1=(double)cd_values[9];
    x2=pow(10.0e0,(double)cd_values[8]);
    th0 = (cd_values[10]==1) ? x1*x2 : x1/x2;
        
    if (flags & H5Z_FLAG_REVERSE)
    {

        nc+=1;
    
        if (NULL==(outBufd = (double *) malloc(sizeof(double)*N)))
        {
            printf("cannot malloc\n");
            goto error;
        }
        
        outBufc=buf[0];
        w = (double *) malloc(sizeof(double)*N);
	M=buf_size[0]/sizeof(doubleint);
	
	for (i=0; i<N; i++)
	  w[i]=0.0e0;
	
	if ((M==1) && (outBufc[0].ii==-1))
	  for (i=0; i<N; i++)
	    outBufd[i]=outBufc[0].w;
	else {  
	  for (i=0; i<M; i++) {
	    w[outBufc[i].ii]=outBufc[i].w;            
	  }
	  alpert_transform(Uj, w, outBufd, -1, part, G, P, J, k2, N);
	}

        *buf = outBufd;
	/*memcpy(*buf, outBufd, sizeof(double)*N);
	free(outBufd);*/
        outBufd = NULL;
        ret_value = M*sizeof(doubleint);
	
	free(w);
    }
    else /* forward filter */
    {
        nc+=1;
    
	f=buf[0];

        if (nbytes > INT32_MAX)
        {
            /* can only compress chunks up to 2GB */
            goto error;
        }        
	
	w = (double *) malloc(sizeof(double)*N);
	alpert_transform(Uj, f, w, 1, part, G, P, J, k2, N);
	
	prctile(f, N, 99.0e0, &mx);
	prctile(f, N, 1.0e0, &mn);
	mxmn = mx - mn;
	thresh(w, N, th0, mxmn, mxmng, &th);
	get_large_indices(w, N, th, mxmng, &M);
	
	if (mxmn==0.0e0)
	  M=1;
	
	outBufc = (doubleint *) malloc(sizeof(doubleint)*M);
	
	if (mxmn==0.0e0) {
	  outBufc[0].w=mx;
	  outBufc[0].ii=-1;
	}
	else {  
	  n=0;
	  for (i=0; i<N; i++) {
	    if (fabs(w[i])/mxmng >= th) {
	      outBufc[n].w=w[i];
	      outBufc[n].ii=i;
	      n+=1;
	    }
	  }
	}

        /**buf = outBufc;*/
	memcpy(*buf, outBufc, sizeof(doubleint)*M);
	free(outBufc);
        *buf_size = M*sizeof(doubleint);
        outBufc = NULL;
        ret_value = M*sizeof(doubleint);
	
	free(w);
    }
    
    done:
    if (outBufd)
        free(outBufd);
    
    if (nc==Ng) {
      nc=0;
      pp = NULL;
      free(p);free(G);free(ki);
      
      for (i=0; i<P; i++)
	free((*part)[i]);
      free(part);
      
      for (i=0; i<nJ; i++)
	free((*Uj)[i]);
      free(Uj);
    }
    
    return ret_value;
    
    error:
    if (outBufc)
        free(outBufc);
    outBufc = NULL;
    if (outBufd)
        free(outBufd);
    outBufd = NULL;
    
    if (nc==Ng) {
      nc=0;
      pp = NULL;
      free(p);free(G);free(ki);
      
      for (i=0; i<P; i++)
	free((*part)[i]);
      free(part);
      
      for (i=0; i<nJ; i++)
	free((*Uj)[i]);
      free(Uj);
    }
    
    return 0;
}
