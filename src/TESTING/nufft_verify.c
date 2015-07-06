#include <stdio.h>
#include <stdlib.h>
#define _GNU_SOURCE         /* See feature_test_macros(7) */
#include <string.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include "nufft.h"

#define	FORMAT_ERROR_TRACE		(-1)
#define	FORMAT_EOF				(-2)
#define	FORMAT_RAW_TRACE			(1)
#define	FORMAT_DETREND_TRACE		(2)
#define	FORMAT_ELAPSED_TIME		(3)
#define	FORMAT_MS_ASSIGN			(4)

#define	MAX_CHARS_PER_LINE		(256)
#define	MAX_ROWS		(10000)


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

/*
 * word counts in a line
 */
int word_count(const char sentence[ ])
{
    int counted = 0; // result

    // state:
    const char* it = sentence;
    int inword = 0;

    do switch(*it) {
        case '\0':
        case ' ': case '\t': case '\n': case '\r': // TODO others?
            if (inword) { inword = 0; counted++; }
            break;
        default: inword = 1;
    } while(*it++);

    return counted;
}

int line_count(FILE *f)
{
	int c=0,n=0;

	// Extract characters from file and store in character c
	for (c = getc(f); c != EOF; c = getc(f))
		if (c == '\n') // Increment count if this character is newline
			n++;
	rewind(f);
	printf("lines=%d\n", n);
	return n;
}

int read_line(FILE *file, char *line)
{
	char *s;
	int len,i;

	while( (s=fgets(line, MAX_CHARS_PER_LINE, file)) != NULL){
		//puts(s);
		if(line[0] == '\n') continue;//skip empty line
		if(line[0] != '#') break;	//not a comment line
	}

	if(s==NULL)
		return FORMAT_EOF;

	len = strlen(s);
	if(strstr(line,"uxj")!=NULL)
		return FORMAT_ELAPSED_TIME;
	if(strstr(line,"raw")!=NULL)
		return FORMAT_RAW_TRACE;
	if(strstr(line,"detrend")!=NULL)
		return FORMAT_DETREND_TRACE;
	if(strstr(line,"ms")!=NULL){
		i=2;	//skip "ms"
		while( line[i] != '=') i++;
		if(i == len)
			return FORMAT_ERROR_TRACE;
		while(i >= 0)	 line[i--]=' ';
		printf("line=%s\n", line);
		return FORMAT_MS_ASSIGN;
	}
	return FORMAT_ERROR_TRACE;
}

int import_file(FILE *file, int *ms,  int *nj, unsigned *uxj, float *cjr, float *dcjr)
{
	int ret=-1;
	int i,nx=-1, nc=-1, nd=-1;
	if(file){
		char line[MAX_CHARS_PER_LINE];
		do{
			ret=read_line(file, line);
			if(ret==FORMAT_EOF){
				printf("EOF.\n");
				break;
			}else if(ret == FORMAT_ERROR_TRACE){
				printf("FORMAT_ERROR_TRACE\n");
				break;
			}else if(ret == FORMAT_MS_ASSIGN){
				*ms=atoi(line);
				printf("ms=%d",*ms);
			}else if(ret == FORMAT_RAW_TRACE){
				i = 0;
				printf("\nFORMAT_RAW_TRACE\n");
				while( fgets(line, MAX_CHARS_PER_LINE, file) != NULL){
					if(line[0]=='\n') {
						printf("empty string, break section\n");
						break;//empty string to break current section
					}
					cjr[i]=atof(line);
					//printf("%d:%f %s\n",i,cjr[i], line);
					i++;
				}
				nc=i;
			}else if(ret == FORMAT_DETREND_TRACE){
				i = 0;
				printf("\nFORMAT_DETREND_TRACE\n");
				while( fgets(line, MAX_CHARS_PER_LINE, file) != NULL){
					if(line[0]=='\n'){
						printf("empty string, break section\n");
						break;//empty string to break current section
					}
					dcjr[i]=atof(line);
					//printf("%d:%f %s\n",i,dcjr[i], line);
					i++;
				}
				nd=i;
			}else if(ret == FORMAT_ELAPSED_TIME){
				i = 0;
				printf("\nFORMAT_ELAPSED_TIME\n");
				while( fgets(line, MAX_CHARS_PER_LINE, file) != NULL){
					//printf("%d:%s\n", i, line);
					if(line[0]=='\n'){
						printf("empty string, break section\n");
						break;//empty string to break current section
					}
					uxj[i]=atoi(line);
					//printf("%d:%u %s\n",i, uxj[i], line);
					i++;
				}
				nx=i;
			}
		}while( (nx ==-1) || (nd ==-1) || (nc ==-1));
		fclose(file);
	}
	if( (nx != nd) || (nx != nc) || (nd != nc)){
		printf("error:nx=%d, nc=%d, nd=%d\n", nx, nc,nd);
	}else {
		*nj=nx;
		ret = 0;
		printf("nj=%d\n", *nj);
	}
	return ret;
}

int verify_trace(const char *filename)
{
	int ms, nj, i, nx=-1, nc=-1, nd=-1;
	int ret=-1;

	int ier=0, iflag=-1;
	float eps=1e-5;
	double err=0.0f;
	double deltaF=0.0f,deltaBPM=0.0f;
	unsigned uxj[MAX_ROWS];
	float xj[MAX_ROWS];
	float cjr[MAX_ROWS], dcjr[MAX_ROWS];
	COMPLEX8 cj[MAX_ROWS];
    COMPLEX8 fk0[MAX_ROWS],fk1[MAX_ROWS];

	FILE *file=fopen(filename,"r");
	if(file){
		if(ret = import_file(file, &ms,  &nj, uxj, cjr, dcjr)){
			printf("import_file fail=%d\n",ret);
			return ret;
		}
	}

	double observingT = uxj[nj-1];
	printf("observingT=%f us\n", observingT);

	/* xj[] is normalized to [-PI,+PI]
	 * (2*PI / T) * (t - T/2)
	 */
	for(int i = 0; i < nj; i ++){
		xj[i]= (uxj[i]- observingT/2.0f) * (2.0f * M_PI / observingT);
		//printf("%f\n", xj[i]);
	}

	/* nufft does not change the observing period, but the sampling number
	 * from N->2M, so only the highest working frequency changes.
	 */
	observingT /= 1000000.0f;//us->1second
	deltaF=1.0f/observingT;
	deltaBPM=60.0f * deltaF;
	printf("observingT=%f s, deltaF=%f, deltaBPM=%fbpm\n", observingT, deltaF, deltaBPM);

	/*
	 * raw trace
	 */
	printf(">>>>raw trace data\n");
	printf("nj=%d, ms=%d, nj=%d\n", nj, ms, nj);
	printf("observingT=%f, dF=%f, dBPM=%f\n", observingT, deltaF, deltaBPM);
	for(int i = 0; i < nj; i ++){
		cj[i].r= cjr[i];
		cj[i].i=0.0f;
	}

	iflag=-1;//forward
	dirft1d1f_(&nj, xj, cj, &iflag, &ms, fk0);
	iflag=-1;	//forward
	nufft1d1ff90_ffte_(&nj, xj, cj, &iflag, &eps, &ms, fk1, &ier);
	errcompf(fk0, fk1, ms, &err);
	printf("*** ier=%d, err=%e\n", ier, err);

	/*
	 * find the max strong Heart beat
	 */
	int mI=0;
	double dMax=0.0f;
	printf("\nnj=%d, ms=%d\n", nj,ms);
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
	printf("nj=%d, ms=%d, nj=%d\n", nj, ms, nj);
	printf("observingT=%f, dF=%f, dBPM=%f\n", observingT, deltaF, deltaBPM);
	for(int i = 0; i < nj; i ++){
		cj[i].r= dcjr[i];
		cj[i].i=0.0f;
	}

	iflag=-1;//forward
	dirft1d1f_(&nj, xj, cj, &iflag, &ms, fk0);
	iflag=-1;	//forward
	nufft1d1ff90_ffte_(&nj, xj, cj, &iflag, &eps, &ms, fk1, &ier);
	errcompf(fk0, fk1, ms, &err);
	printf("*** ier=%d, err=%e\n", ier, err);

	/*
	 * find the max strong Heart beat
	 */
	mI=0;
	dMax=0.0f;
	printf("\nnj=%d, ms=%d\n", nj,ms);
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

	return ret;
}

/*
 *
 */
int main(int argc, char **argv)
{

	if(argc >=2){
		verify_trace(argv[1]);
	}else{
		printf("%s export.txt\n",argv[0]);
	}




}
#if 0

		if(lc >0){
			char str[100];
			if(fgets(str, 99, f) != NULL){
				int channels;
				rewind(f);//reset file location to start after a trial read!!!
				channels = word_count(str);
				if(channels >=1 && channels <= MAX_CHANNELS){//limite 100 channels
					if(raw_buf=	malloc(lc*channels*sizeof(float))){
						int lda;
						int order=CblasColMajor;
						pr_debug(DSP_INFO, "channels=%d, lc=%d\n",channels,lc);
						//row-major stored : channels * N
						//col-major stored : N * channels
						lda=lc;
						if(CblasColMajor==order){
							*M=lc;
							*N=channels;
						}else{//row-major
							*M=channels;
							*N=lc;
						}
						pr_debug(DSP_INFO, "m=%d, n=%d, lda=%d\n", *M , *N, lda);
						int i,j,ret;
						//read data by major or by col
						for(i=0;i < lc; i ++){
							for(j = 0; j < channels; j ++){
								if(CblasColMajor==order){
									ret=fscanf(f,"%f",&raw_buf[j*lda+i]);
									//pr_debug(DSP_INFO, "%7.3f ",raw_buf[j*lda+i]);
								}else{
									ret=fscanf(f,"%f",&raw_buf[i*lda + j]);
									//pr_debug(DSP_INFO, "%7.3f ",raw_buf[i*lda + j]);
								}
								if(EOF==ret)
									break;
							}
							//putchar('\n');
							//pr_debug(DSP_INFO,"\n");
						}
					}else{
						pr_debug(DSP_ERROR,"%s:raw_buf malloc error!\n", __func__);
					}
				}else{
					pr_debug(DSP_ERROR,"%s:channels [%d] exceeds limit %d or error!\n", \
					__func__, channels, MAX_CHANNELS);
				}
			}else{
				pr_debug(DSP_ERROR,"%s:fgets error!\n", __func__);
			}
		}else{
			pr_debug(DSP_ERROR,"%s:line counts error!\n", __func__);
		}
#endif
