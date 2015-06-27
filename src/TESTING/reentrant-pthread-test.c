#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
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
		//printf("%d:fk0=(%f,%f) fk1=(%f,%f)\n",k, fk0[k].r, fk0[k].i, fk1[k].r,fk1[k].i);
	}
	*err =sqrt(ealg/salg);
	//printf("ealg/salg:%g/%g err=%e\n",ealg, salg, *err);
}

#define	MX	(512)

//http://forum.cvapp.org/viewtopic.php?f=21&t=431&sid=4e8c86381ce6b3a43a53d7d8f8b22058
void *nufft_sine(void *ptr)
{
	int ier,iflag,j,k1,ms,nj;
	float xj[MX];
	float eps=1e-5;
	double err;
	COMPLEX8 cj[MX];
	COMPLEX8 fk0[MX],fk1[MX];
	double observingT=0.0f;
	int w=*(int *)ptr;
	char filename[80];
	FILE *file;
	if(ptr) free(ptr);
	sprintf(filename,"log-%d",w);
	printf("+%s:w=%f\n",__func__,(float)w);
	file = fopen(filename, "w");
	if(file == NULL){
		printf("%s opening failure\n", filename);
		return NULL;
	}
	fprintf(file,"+%s:w=%f\n",__func__,(float)w);
#if 1
//
//     --------------------------------------------------
//     create some test data
//     --------------------------------------------------
	ms = 36;	//oversampling 2*ms in time domain
	nj = 57;
	//fortran array index from "1", but c is from "0"
	printf("input xj, nj and cj\n");
	printf("nj bin=%f\n", M_PI*2.0/(nj));
	for(k1 = -nj/2; k1 <= (nj-1)/2;k1++){
		j = k1+nj/2;
		//xj[j] = M_PI * cos(-M_PI*(j+1)/nj);
		//float r_s=M_PI*2.0/(nj*2) * ((rand()%10 - 5 )/5.0f);
		float r_s=M_PI*2.0/(nj) * ((rand()%10 - 5.0 )/5.0f);
		xj[j] = M_PI * k1/(nj/2) + r_s;//generate a random non-equispaced
		//cj[j] = cmplxf(cos(5.0*xj[j])+cos(1.0*xj[j]),
		//			   sin(2.0*xj[j]));
		//cj[j].r = cos(5.0*xj[j])+cos(1.0*xj[j]);
		//cj[j].i = sin(2.0*xj[j]);
		cj[j].r = cos(w*xj[j]);

		//printf("%d,%d:%f, %f,(%f,%f)\n", k1, j, r_s, xj[j], cj[j].r, cj[j].i);
	}
	printf("dirft1d1f_\n");
	iflag=-1;	//forward
	dirft1d1f_(&nj, xj, cj, &iflag, &ms, fk0);
	printf("nufft1d1ff90_ffte_\n");
	iflag=-1;	//forward
	nufft1d1ff90_ffte_(&nj, xj, cj, &iflag, &eps, &ms, fk1, &ier);
	errcompf(fk0, fk1, ms, &err);
	printf("ier=%d, err=%e\n\n", ier, err);
	printf("nj=%d, ms=%d\n",nj,ms);
	float fmax0=0.0f, fmax1=0.0f;
	int I0=0,I1=0;
	for(int i = 0; i < ms; i ++){
		/*printf("%d:fk0=(%f,%f),fk1=(%f,%f)\n", i, fk0[i].r,fk0[i].i,
			fk1[i].r,fk1[i].i);
		printf("(%f,%f)\n",cabsf_sq(fk0[i]), cabsf_sq(fk1[i]));*/
		if(cabsf_sq(fk0[i]) > fmax0){
			fmax0=cabsf_sq(fk0[i]);
			I0=i;
		}
		if(cabsf_sq(fk1[i]) > fmax1){
			fmax1=cabsf_sq(fk1[i]);
			I1=i;
		}
	}
	printf("****%s[%d]:I0=%d/%d, max=%f\n", __func__, w, ms/2-I0, ms, fmax0);
	printf("****%s[%d]:I1=%d/%d, max=%f\n", __func__, w, ms/2-I1, ms, fmax1);

	fprintf(file, "****%s[%d]:I0=%d/%d, max=%f\n", __func__, w, ms/2-I0, ms, fmax0);
	fprintf(file, "****%s[%d]:I1=%d/%d, max=%f\n", __func__, w, ms/2-I1, ms, fmax1);
#endif
	printf("-%s:%d\n",__func__, w);

	fprintf(file, "-%s:%d\n",__func__, w);

	fclose(file);
	pthread_exit(NULL);
	//return NULL;
}

#define	MAX_THREADS	(10)
/*
 * I want to verify if ffte and nufft are reentrant.
 * !!!After the verifications, nufft and ffte are indeed "REENTRANT".
 */
int main(int argc, char **argv)
{
	srand(time(NULL));

	pthread_t thread_id[MAX_THREADS];

	for(int i = 0; i < MAX_THREADS ; i ++){
		int ret;
		/*
		 * arg passed to pthead should be "local" to the thread,
		 * so arg won't be changed in the main thread before the pthread executes.
		 * if i is used, and &i is passed to pthread,
		 * it could be that i is changed in the iteration, but the pthread which is created
		 * has not been executed yet. When the created pthread is scheduled to execute,
		 * the i has been changed.
		 */
		int *id=malloc(sizeof(int));
		*id=i;
		ret=pthread_create (thread_id+i, NULL, &nufft_sine, id);
		printf("pthread_create:%d, ret=%d\n",i,ret);
	}

	/* Make sure
	the first thread has finished. */
	for(int i = 0 ; i < MAX_THREADS ; i ++){
		printf(">>pthread_join:%d\n",i);
		pthread_join(thread_id[i], NULL);
		printf("<<pthread_join:%d\n",i);
	}

	//nufft_hr();
}
