#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include "nufft.h"
#include "r110_m156.h"

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
		printf("%d:fk0=(%f,%f) fk1=(%f,%f)\n",k, fk0[k].r, fk0[k].i, fk1[k].r,fk1[k].i);
	}
	*err =sqrt(ealg/salg);
	printf("ealg/salg:%g/%g err=%e\n",ealg, salg, *err);
}

#define	MX	(512)

void nufft_hr(void)
{
	int nj = 51;
	int ier=0, iflag=-1;
	int ms= 32;
	float eps=1e-5;
	double err=0.0f;
	double observingT=0.0f;
	double deltaF=0.0f,deltaBPM=0.0f;
	COMPLEX8 cj[MX];
    COMPLEX8 fk0[MX],fk1[MX];
	float xj[MX];
/*
 * xj : non-equspaced x position
 * cj : y position
 */
	int size_x=sizeof(uxj)/sizeof(unsigned);
	int size_raw=sizeof(rcj)/sizeof(float);
	int size_det=sizeof(dcj)/sizeof(float);
	printf("size_x=%d,size_raw=%d, size_det=%d\n", size_x, size_raw, size_det);
	nj=size_x;
	observingT = uxj[size_x-1];
	printf("observingT=%f us\n", observingT);

	/* xj[] is normalized to [-PI,+PI]
	 * (2*PI / T) * (t - T/2)
	 */
	for(int i = 0; i < size_x; i ++){
		xj[i]= (uxj[i]- observingT/2.0f) * (2.0f * M_PI / observingT);
		printf("%f\n", xj[i]);
	}

	/* nufft does not change the observing period, but the sampling number
	 * from N->2M, so only the highest working frequency changes.
	 */
	observingT /= 1000000.0f;//us->1second
	deltaF=1.0f/observingT;
	deltaBPM=60.0f * deltaF;
	printf("observingT=%f s, deltaF=%f, deltaBPM=%fbpm\n", observingT, deltaF, deltaBPM);

	/* ms points in frequency domain, so at least
	 * 180.0bpm/deltaBPM are needed.
	 * ms is oversampling 2*ms in time domain, but
	 * ms points in frequency domain.
	 * half of ms, ms/2, are mirrored, so 2*ms are needed.
	 */
	//double dms=ceil(180.0f/deltaBPM)*2.0f;
	//ms = next235_(&dms);
	ms=60;

	/*
	 * raw trace
	 */
	printf(">>>>raw trace data\n");
	printf("nj=%d, ms=%d, size_raw=%d\n", nj, ms, size_raw);
	printf("observingT=%f, dF=%f, dBPM=%f\n", observingT, deltaF, deltaBPM);
	for(int i = 0; i < size_raw; i ++){
		cj[i].r= rcj[i];
		cj[i].i=0.0f;
	}

	iflag=-1;//forward
	dirft1d1f_(&nj, xj, cj, &iflag, &ms, fk0);
	iflag=-1;	//forward
	nufft1d1ff90_ffte_(&nj, xj, cj, &iflag, &eps, &ms, fk1, &ier);
	errcompf(fk0, fk1, ms, &err);
	printf("ier=%d, err=%e\n", ier, err);

	int mI=0;
	double dMax=0.0f;
	printf("nj=%d, ms=%d\n", nj,ms);
	printf("observingT=%f, dF=%f, dBPM=%f\n", observingT, deltaF, deltaBPM);
	for(int i = ms/2; i < ms; i ++){
		float bpm;
		bpm=(i-ms/2)*deltaBPM;
		printf("\n%d:%f,%fbpm\n", i , (i-ms/2)*deltaF, bpm);

		printf("fk0=(%f,%f),fk1=(%f,%f)\n", fk0[i].r,fk0[i].i,
			fk1[i].r,fk1[i].i);
		printf("(%f,%f, %f)\n",cabsf_sq(fk0[i]), cabsf_sq(fk1[i]), dMax);

		if( (bpm < 180.0f) && (bpm > 50.0f) && (cabsf_sq(fk1[i]) > dMax) ){
			dMax =cabsf_sq(fk1[i]);
			mI = i;
		}
	}
	printf("*********** mI=%d,bpm=%f, max=%f\n", mI, (mI - ms/2)*deltaBPM, dMax);

	/*
	 * detrend trace
	 */
	printf("\n\n>>>>detrend data\n");
	printf("nj=%d, ms=%d, size_det=%d\n", nj, ms, size_det);
	printf("observingT=%f, dF=%f, dBPM=%f\n", observingT, deltaF, deltaBPM);
	for(int i = 0; i < size_det; i ++){
		cj[i].r= dcj[i];
		cj[i].i=0.0f;
	}

	iflag=-1;//forward
	dirft1d1f_(&nj, xj, cj, &iflag, &ms, fk0);
	iflag=-1;	//forward
	nufft1d1ff90_ffte_(&nj, xj, cj, &iflag, &eps, &ms, fk1, &ier);
	errcompf(fk0, fk1, ms, &err);
	printf("ier=%d, err=%e\n", ier, err);

	mI=0;
	dMax=0.0f;
	printf("nj=%d, ms=%d\n", nj,ms);
	printf("observingT=%f, dF=%f, dBPM=%f\n", observingT, deltaF, deltaBPM);
	for(int i = ms/2; i < ms; i ++){
		float bpm;
		bpm=(i-ms/2)*deltaBPM;
		printf("\n%d:%f,%fbpm\n", i , (i-ms/2)*deltaF, bpm);

		printf("fk0=(%f,%f),fk1=(%f,%f)\n", fk0[i].r,fk0[i].i,
			fk1[i].r,fk1[i].i);
		printf("(%f,%f, %f)\n",cabsf_sq(fk0[i]), cabsf_sq(fk1[i]), dMax);

		if( (bpm < 180.0f) && (bpm > 50.0f) && (cabsf_sq(fk1[i]) > dMax) ){
			dMax =cabsf_sq(fk1[i]);
			mI = i;
		}
	}
	printf("************ mI=%d,bpm=%f, max=%f\n", mI, (mI - ms/2)*deltaBPM, dMax);
}

int main(int argc, char **argv)
{
	srand(time(NULL));
	nufft_hr();
}