#######################################################################
#  This is the makefile to create a library for the FFTBASE.
#  The files are grouped as follows:
#
#       SFFTBASE1 -- Single precision real FFTBASE routines
#       CFFTBASE1 -- Single precision complex FFTBASE routines
#       DFFTBASE1 -- Double precision real FFTBASE routines
#       ZFFTBASE1 -- Double precision complex FFTBASE routines
#
#       CB1AUX -- Real FFTBASE routines called by complex routines
#       ZB1AUX -- D.P. real FFTBASE routines called by d.p. complex
#                 routines
#
#      ALLFFTBASE -- Auxiliary routines for Level 2 and 3 FFTBASE
#
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  To create or add to the library, enter make
#  followed by one or more of the precisions desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Note that these commands are not safe for parallel builds.
#
#  Alternatively, the commands
#       make all
#  or
#       make
#  without any arguments creates a library of all four precisions.
#  The name of the library is held in FFTBASELIB, which is set in the
#  top-level make.inc
#
#  To remove the object files after the library is created, enter
#       make clean
#  To force the source files to be recompiled, enter, for example,
#       make single FRC=FRC
#
#---------------------------------------------------------------------
#
#  Edward Anderson, University of Tennessee
#  March 26, 1990
#  Susan Ostrouchov, Last updated September 30, 1994
#  ejr, May 2006.
#
#######################################################################

#---------------------------------------------------------
#  Comment out the next 6 definitions if you already have
#  the Level 1 FFTBASE.
#---------------------------------------------------------
#set(SBLAS1 isamax.f sasum.f saxpy.f scopy.f sdot.f snrm2.f
#	srot.f srotg.f sscal.f sswap.f sdsdot.f srotmg.f srotm.f)

#set(CBLAS1 scabs1.f scasum.f scnrm2.f icamax.f caxpy.f ccopy.f
#	cdotc.f cdotu.f csscal.f crotg.f cscal.f cswap.f csrot.f)

#set(DBLAS1 idamax.f dasum.f daxpy.f dcopy.f ddot.f dnrm2.f
#	drot.f drotg.f dscal.f dsdot.f dswap.f drotmg.f drotm.f)
#nufft*d_demof90.f nufft*df90.f dirft*d.f dfftpack.f next235.f
set(ZFFTF90 ../dfftpack.f)

#nufft*d_demof77.f nufft*df77.f dirft*d.f dfft.f next235.f
set(ZFFTF77 ../dfft.f)

# default build all of it
set(ALLOBJF90 ${ZFFTF90})
set(ALLOBJF77 ${ZFFTF77})

#option(FFTBASE_F90 "Build FFTBASE_F90" ON)

#option(FFTBASE_SINGLE "Build FFTBASE SINGLE" ON)
#option(FFTBASE_DOUBLE "Build FFTBASE DOUBLE" ON)
#option(FFTBASE_COMPLEX "Build FFTBASE COMPLEX" ON)
#option(FFTBASE_COMPLEX16 "Build FFTBASE COMPLEX16" ON)

if(FFTBASE_SINGLE)
  set(ALLOBJF90 ${ZFFTF90})
  set(ALLOBJF77 ${ZFFTF77})
endif()
if(FFTBASE_DOUBLE)
  set(ALLOBJF90 ${ZFFTF90})
  set(ALLOBJF77 ${ZFFTF77})
endif()
if(FFTBASE_COMPLEX)
  set(ALLOBJF90 ${ZFFTF90})
  set(ALLOBJF77 ${ZFFTF77})
endif()
if(FFTBASE_COMPLEX16)
  set(ALLOBJF90 ${ZFFTF90})
  set(ALLOBJF77 ${ZFFTF77})
endif()

if(NUFFT_F90)
add_library(fftbase ${ALLOBJF90})
else()
add_library(fftbase ${ALLOBJF77})
endif()

if(UNIX)
  target_link_libraries(fftbase m)
endif()

if(ANDROID_NDK)
  target_link_libraries(fftbase m)
endif(ANDROID_NDK)

#target_link_libraries(fftbase)
nufft_install_library(fftbase)
