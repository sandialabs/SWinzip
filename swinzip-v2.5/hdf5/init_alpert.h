#include <stdio.h>

void get_Nbins(int *ki, int N, int d,int *J, int *k2, int *P, int *nJ);
void dichotomic_grouping(double *p, int N, int d, int *G, int *(*part)[], int J, int k2, int P);
void compute_Uj(double *p, int *(*part)[], int *ki, int *G, int k2, int P, int J, int d, double *(*Uj)[]);
int compqi (const void *a, const void *b);
