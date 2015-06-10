#ifndef DFFTPACK_HPP
#define DFFTPACK_HPP

#include <complex>

extern "C" void dfftb_(int *, double *, double *);
extern "C" void dfftf_(int *, double *, double *);
extern "C" void dffti_(int *, double *);

extern "C" void dzffti_(int *, double *);
extern "C" void dzfftf_(int *, double *, double *, double *, double *, double *);
extern "C" void dzfftb_(int *, double *, double *, double *, double *, double *);

extern "C" void dsinti_(int *, double *);
extern "C" void dsint_(int *, double *, double *);

extern "C" void dcosti_(int *, double *);
extern "C" void dcost_(int *, double *, double *);

extern "C" void dsinqi_(int *, double *);
extern "C" void dsinqf_(int *, double *, double *);
extern "C" void dsinqb_(int *, double *, double *) ;

extern "C" void dcosqi_(int *, double *);
extern "C" void dcosqf_(int *, double *, double *);
extern "C" void dcosqb_(int *, double *, double *);

extern "C" void zffti_(int *, double *);
extern "C" void zfftf_(int *, std::complex<double> *, double *);
extern "C" void zfftb_(int *, std::complex<double> *, double *);

#endif
