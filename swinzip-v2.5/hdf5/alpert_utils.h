#include <stdio.h>

struct pos {
    int posId;
    double posx;
};

void dgeqrf_(int *rows, int *cols, double *matA, int *LDA, double *TAU, double *WORK, int *LWORK, int *INFO);
void dorgqr_(int *M, int *N, int *K, double *A, int *LDA, double *TAU, double *WORK, int *LWORK, int *INFO);
int dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc); 

void transpose(double *M, int N);
int compq (const void *a, const void *b);

void prctile(double *f, int N, double p, double *c);
void fmean(double *f, int N, double mxmn, double *m);
void fmeana(double *f, int N, double x, double mxmn, double *m);
void mxmn(double *f, int N, double *c);