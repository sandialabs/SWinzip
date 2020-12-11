#include <stdio.h>

void alpert_transform(double *(*Uj)[], double *f, double *r, int dire, int *(*part)[], int *G, int P, int J, int k2, int N);
void thresh(double *w, int N, double th0, double mxmn, double mxmng, double *th);
void get_large_indices(double *w, int N, double th, double mxmn, int *M);

