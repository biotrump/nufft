#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include "nufft.h"

#ifndef	M_PI
# define M_PI       3.14159265358979323846  /* pi */
#endif
#ifndef M_PI_2
# define M_PI_2     1.57079632679489661923  /* pi/2 */
#endif
#ifndef M_PI_4
# define M_PI_4     0.78539816339744830962  /* pi/4 */
#endif

double cabs_sq(COMPLEX16 c)
{
	return c.r * c.r + c.i * c. i;
}

float cabsf_sq(COMPLEX8 c)
{
	return c.r * c.r + c.i * c. i;
}

double cabs_(COMPLEX16 c)
{
	return sqrt(c.r * c.r + c.i * c. i);
}

float cabsf_(COMPLEX8 c)
{
	return sqrtf(c.r * c.r + c.i * c. i);
}

COMPLEX8 cmplxf(float r, float i)
{
	COMPLEX8 ctemp;
	ctemp.r=r;
	ctemp.i=i;
	return ctemp;
}

COMPLEX16 cmplx(double r, double i)
{
	COMPLEX16 ctemp;
	ctemp.r=r;
	ctemp.i=i;
	return ctemp;
}

void errcompf(PCOMPLEX8 fk0, PCOMPLEX8 fk1, int n, double *err)
{
	double salg=0.0f,ealg=0.0f;

	for(int k = 0 ; k < n; k ++){
		COMPLEX16 ctemp;
		ctemp.r = fk1[k].r-fk0[k].r;
		ctemp.i = fk1[k].i-fk0[k].i;
		ealg += cabs_sq(ctemp);
		ctemp.r=(double)fk0[k].r;
		ctemp.i=(double)fk0[k].i;
		salg += cabs_sq(ctemp);
		//ealg += ( (fk1[k].r-fk0[k].r)*(fk1[k].r-fk0[k].r) +
		//		(fk1[k].i-fk0[k].i)*(fk1[k].i-fk0[k].i) );
		//salg +=  (fk0[k].r*fk0[k].r + fk0[k].i*fk0[k].i);
		printf("(%f,%f) (%f,%f)\n",fk1[k].r,fk1[k].i,fk0[k].r, fk0[k].i);
	}
	*err =sqrt(ealg/salg);
	printf("ealg/salg:%g/%g err=%e\n",ealg, salg, *err);
}

#define	MX	(100000)

int main(int argc, char **argv)
{
      int i,ier,iflag,j,k1,ms,nj;
      float xj[MX], sk[MX];
      float eps;
      double err;
      COMPLEX8 cj[MX],cj0[MX],cj1[MX];
      COMPLEX8 fk0[MX],fk1[MX];
//
//     --------------------------------------------------
//     create some test data
//     --------------------------------------------------
      ms = 90;
      nj = 128;
	  //fortran array index from "1", but c is from "0"
      for(k1 = -nj/2; k1 <= (nj-1)/2;k1++){
         j = k1+nj/2;
         xj[j] = M_PI * cos(-M_PI*(j+1)/nj);
         cj[j] = cmplxf( sin(M_PI*(j+1)/nj), cos(M_PI*(j+1)/nj));
		 printf("%f, (%f,%f)\n", xj[j], cj[j].r, cj[j].i);
	  }
//
//     --------------------------------------------------
//     start tests
//     --------------------------------------------------
//
      iflag = 1;
      printf("Start 1D testing: nj = %d  ms =%d\n",nj, ms);
      for(i = 0; i < 5; i ++){
        if (i==0) eps=1e-3;
        if (i==1) eps=1e-4;
        if (i==2) eps=1e-5;
        if (i==3) eps=1e-6;
        if (i==4) eps=1e-7;
		printf("Requested precision eps =%E\n",eps);
//
//     -----------------------
//    1D Type1 method
//     -----------------------
//
        dirft1d1f_(&nj, xj, cj, &iflag, &ms, fk0);
        nufft1d1ff90_ffte_(&nj, xj, cj, &iflag, &eps, &ms, fk1, &ier);
        errcompf(fk0, fk1, ms, &err);
        printf(" ier = %d  type 1 error = %E\n", ier, err);
//
//     -----------------------
//    1D Type2 method
//     -----------------------
//
        dirft1d2f_(&nj, xj, cj0, &iflag, &ms, fk0);
        nufft1d2ff90_ffte_(&nj,xj,cj1,&iflag, &eps, &ms, fk0, &ier);
        errcompf(cj0, cj1, nj, &err);
        printf(" ier = %d  type 2 error = %E\n",ier, err);
//
//     -----------------------
//    1D Type3 method
//     -----------------------
		for(k1 = 0; k1 < ms; k1++){
			sk[k1] = 48.0*cos(k1*M_PI/ms);
		}
        dirft1d3f_(&nj, xj, cj, &iflag, &ms, sk, fk0);
        nufft1d3ff90_ffte_(&nj, xj, cj, &iflag, &eps, &ms, sk, fk1, &ier);
        errcompf(cj0, cj1, nj, &err);
		printf(" ier = %d  type 3 error = %E\n",ier, err);
	}
}
