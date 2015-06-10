#ifndef DFFTPACK_H
#define DFFTPACK_H

#include <complex.h>

extern void dcfftb_(int *, double *, double *);
extern void dcfftf_(int *, double *, double *);
extern void dcffti_(int *, double *);

extern void dzffti_(int *, double *);
extern void dzfftf_(int *, double *, double *, double *, double *, double *);
extern void dzfftb_(int *, double *, double *, double *, double *, double *);

extern void dsinti_(int *, double *);
extern void dsint_(int *, double *, double *);

extern void dcosti_(int *, double *);
extern void dcost_(int *, double *, double *);

extern void dsinqi_(int *, double *);
extern void dsinqf_(int *, double *, double *);
extern void dsinqb_(int *, double *, double *) ;

extern void dcosqi_(int *, double *);
extern void dcosqf_(int *, double *, double *);
extern void dcosqb_(int *, double *, double *);

extern void zffti_(int *, double *);
extern void zfftf_(int *, double complex *, double *);
extern void zfftb_(int *, double complex *, double *);

#endif
