#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "hdf5.h"

#define H5FILE_NAME        "SDS_3d.h5"

int main (int argc, char *argv[])
{
    hid_t       file, dataset;         /* file and dataset handles */
    hid_t       plist;
    hid_t       datatype, dataspace;   /* handles */
    hid_t       grp;
    
    hsize_t     dimsd[3];              /* dataset dimensions */
    hsize_t     cdimd[3];
    herr_t      status;
    int         i, j,l,c, ex,k;
    int rank;
    double mx, mn, mxmnt, th,x,y,z;
    double pi=4.0e0*atan(1.0e0);
    
    size_t      size;
    int         numfilt;
    unsigned flags, filter_info;
    H5Z_filter_t filter_type;
    size_t   nelmts;
    
    htri_t avail;
    unsigned filter_config;
    
    double *d, *dr;
    
    int cd_nelmts;
    unsigned int *cd_values;
    
    int nx,ny,nz;
    
    /* Define a 3D regular grid */
    nx=240;
    ny=180;
    nz=200;
    
    /* Define a toy functional field comprised of sines and cosines */
    d = (double *) malloc(sizeof(double)*nx*ny*nz);
    k=0;
    for (i=0;i<nx;i++) {
      for (j=0;j<ny;j++) {
	for (l=0;l<nz;l++) {
	  x=(double)i/(double)(nx-1);
	  y=(double)j/(double)(ny-1);
	  z=(double)l/(double)(nz-1);
	  d[k]=48.0e0*(sin(2.0e0*pi*x) - sin(2.0e0*pi*y) + cos(3.0e0*pi*z))*sin(2.0e0*pi*x)*cos(3.0e0*pi*z)-52.0e0;
	  k+=1;
	}
      }
    }
    
    /* Find the max and min of the field to be used in the compression inout */
    mx=d[0];
    mn=d[0];
    for (i=0; i<nx*ny*nz; i++) {
      printf("%g\n",d[i]);          /* Output the original field data */
      
      if (d[i]>=mx)
	mx=d[i];
      
      if (d[i]<=mn)
	mn=d[i];
    }
    
    dimsd[0]=nx;
    dimsd[1]=ny;
    dimsd[2]=nz;
    
    /* Writing Alpert compressed H5 file */
    file = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    dataspace = H5Screate_simple(3, dimsd, NULL);
    plist  = H5Pcreate (H5P_DATASET_CREATE);
    
    /* The size of each chunk in all directions */
    cdimd[0] = 40;
    cdimd[1] = 90;
    cdimd[2] = 50;
    
    /* cd_values is the parameters input to Alpert compressed, it contains 14 elements */
    cd_nelmts=14;
    cd_values = (unsigned int *) malloc(sizeof(unsigned int)*cd_nelmts);
    
    /* These are the elements of the parameters cd_values */
    cd_values[0]=3;                  /* wavelet order */
    cd_values[1]=3;                  /* number of dimensions */
    cd_values[2]=cdimd[0];           /* number of points per chunk in the 1st dimension */
    cd_values[3]=cdimd[1];           /* number of points per chunk in the 2nd dimension */
    cd_values[4]=cdimd[2];           /* number of points per chunk in the 3rd dimension */    
    
    /* some way to convert a double to two integers, here we convert the max-min */
    ex=ceil(log10(mx-mn))-1;
    cd_values[5]=abs(ex);
    cd_values[6]=ceil((mx-mn)/pow(10.0e0,(double)ex));
        if (ex>=0)
      cd_values[7]=1;
    else 
      cd_values[7]=0;
    
    /* some way to convert a double to two integers, here we convert the truncation threshold */
    th=0.005e0;
    ex=ceil(log10(th))-1;
    cd_values[8]=abs(ex);
    cd_values[9]=ceil(th/pow(10.0e0,(double)ex));
        if (ex>=0)
      cd_values[10]=1;
    else 
      cd_values[10]=0;
    
    cd_values[11]=nx;
    cd_values[12]=ny;
    cd_values[13]=nz;
    
    for (i=0; i<cd_values[1]; i++) {
      if (dimsd[i]%cdimd[i]!=0)
	printf("Warning: uneven chunks in dimension %d, results may be erroneous\n",i+1);
    }
      
    status = H5Pset_chunk (plist, 3, cdimd);
    
    /* Set the Alpert filter to perform the compression upun writing */
    status = H5Pset_filter (plist, H5Z_FILTER_ALPERT, H5Z_FLAG_MANDATORY, (size_t)cd_nelmts, cd_values);
    if (status < 0) exit(1);
        
    datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    status = H5Tset_order(datatype, H5T_ORDER_LE);
         
    dataset = H5Dcreate2(file, "data", datatype, dataspace, H5P_DEFAULT, plist, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, d);
    
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Pclose(plist);
    H5Fclose(file);
    
    free(d);
    
    /* Done with writing, now reading Alpert compressed H5 file */
    file = H5Fopen(H5FILE_NAME, H5F_ACC_RDONLY, H5P_DEFAULT);
    dataset = H5Dopen2(file, "data", H5P_DEFAULT);
    datatype  = H5Dget_type(dataset);
    size  = H5Tget_size(datatype);    
    plist = H5Dget_create_plist (dataset);
    numfilt = H5Pget_nfilters (plist);
    nelmts = 0;
    filter_type = H5Pget_filter2 (plist, 0, &flags, &nelmts, NULL, 0, NULL, &filter_info);
    dataspace = H5Dget_space(dataset);    /* dataspace handle */
    rank      = H5Sget_simple_extent_ndims(dataspace);
    status  = H5Sget_simple_extent_dims(dataspace, dimsd, NULL);
    dr = (double *) malloc(sizeof(double)*nx*ny*nz);
    status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dr);
    
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Pclose(plist);
    H5Fclose(file);
    
    /* Output the decompressed field data */
    for (i=0; i<nx*ny*nz; i++)
      printf("%g\n",dr[i]);
    
    free(dr);
    free(cd_values);
    
    return 0;
}

