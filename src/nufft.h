#ifndef _H_NUFFT_H
#define	_H_NUFFT_H

#ifndef _H_BCV_BASE_TYPE_H
typedef struct _complex8{
	float r;
	float i;
}COMPLEX8, *PCOMPLEX8;

typedef struct _complex16{
	double r;
	double i;
}COMPLEX16, *PCOMPLEX16;
#endif

//return the number next base which is 2^p*3^q*5^r, (p,q,r>=0)
int next235_(double *base);

/*
 * single precision, float,
 *
 */
void dirft1d1f_(int *nj, float *xj, PCOMPLEX8 cj, int *iflag,
				int *ms, PCOMPLEX8 fk);
void dirft1d2f_(int *nj, float *xj, PCOMPLEX8 cj, int *iflag,
				int *ms, PCOMPLEX8 fk);
void dirft1d3f_(int *nj, float *xj, PCOMPLEX8 cj, int *iflag,
				int *nk, float *sk, PCOMPLEX8 fk);

/* 1D type1 transformations */
void nufft1d1ff90_ffte_(int *nj, float *xj, PCOMPLEX8 cj, int *iflag, float *eps,
					   int *ms, PCOMPLEX8 fk, int *ier);
/* 1D type2 transformations */
void nufft1d2ff90_ffte_(int *nj, float *xj, PCOMPLEX8 cj, int *iflag, float *eps,
					   int *ms, PCOMPLEX8 fk, int *ier);
/* 1D type3 transformations */
void nufft1d3ff90_ffte_(int *nj, float *xj, PCOMPLEX8 cj, int *iflag, float *eps,
					   int *nk, float *sk, PCOMPLEX8 fk, int *ier);

/*
 * double precision
 *
 */
void dirft1d1_(int *nj, double *xj, PCOMPLEX16 cj, int *iflag, int *ms,PCOMPLEX16 fk);
void dirft1d2_(int *nj, double *xj, PCOMPLEX16 cj, int *iflag, int *ms,PCOMPLEX16 fk);
void dirft1d3_(int *nj, double *xj, PCOMPLEX16 cj, int *iflag, int *nk,
			  double *sk, PCOMPLEX16 fk);

void nufft1d1f90_ffte_(int *nj, double *xj, PCOMPLEX16 cj,int * iflag,
					  double *eps, int *ms, PCOMPLEX16 fk, int *ier);
void nufft1d2f90_ffte_(int *nj, double *xj, PCOMPLEX16 cj, int *iflag, double *eps,
					   int *ms, PCOMPLEX16 fk, int *ier);
void nufft1d3f90_ffte_(int *nj, double *xj, PCOMPLEX16 cj, int *iflag, double *eps,
					   int *nk, double *sk, PCOMPLEX16 fk, int *ier);
#endif
