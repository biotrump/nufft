#######################################################################
#  This is the makefile to create a library for the NUFFT.
#  The files are grouped as follows:
#
#       SNUFFT1 -- Single precision real NUFFT routines
#       CNUFFT1 -- Single precision complex NUFFT routines
#       DNUFFT1 -- Double precision real NUFFT routines
#       ZNUFFT1 -- Double precision complex NUFFT routines
#
#       CB1AUX -- Real NUFFT routines called by complex routines
#       ZB1AUX -- D.P. real NUFFT routines called by d.p. complex
#                 routines
#
#      ALLNUFFT -- Auxiliary routines for Level 2 and 3 NUFFT
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
#  The name of the library is held in NUFFTLIB, which is set in the
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
#  the Level 1 NUFFT.
#---------------------------------------------------------
#set(SBLAS1 isamax.f sasum.f saxpy.f scopy.f sdot.f snrm2.f
#	srot.f srotg.f sscal.f sswap.f sdsdot.f srotmg.f srotm.f)

#set(CBLAS1 scabs1.f scasum.f scnrm2.f icamax.f caxpy.f ccopy.f
#	cdotc.f cdotu.f csscal.f crotg.f cscal.f cswap.f csrot.f)

#set(DBLAS1 idamax.f dasum.f daxpy.f dcopy.f ddot.f dnrm2.f
#	drot.f drotg.f dscal.f dsdot.f dswap.f drotmg.f drotm.f)
#nufft*d_demof90.f nufft*df90.f dirft*d.f dfftpack.f next235.f
set(ZNUFFTF90 next235.f dirft1d.f dirft2d.f dirft3d.f
	nufft1df90.f nufft2df90.f nufft3df90.f)

#nufft*d_demof77.f nufft*df77.f dirft*d.f dfft.f next235.f
set(ZNUFFTF77 next235.f dirft1d.f dirft2d.f dirft3d.f
	nufft1df77.f nufft2df77.f nufft3df77.f)

#set(CB1AUX  isamax.f sasum.f saxpy.f scopy.f snrm2.f sscal.f)
#set(ZB1AUX idamax.f dasum.f daxpy.f dcopy.f dnrm2.f dscal.f)

#---------------------------------------------------------------------
#  The following line defines auxiliary routines needed by both the
#  Level 2 and Level 3 BLAS.  Comment it out only if you already have
#  both the Level 2 and 3 BLAS.
#---------------------------------------------------------------------
#set(ALLBLAS  lsame.f xerbla.f xerbla_array.f)

# default build all of it
set(ALLOBJF90 ${ZNUFFTF90})
set(ALLOBJF77 ${ZNUFFTF77})

#option(NUFFT_SINGLE "Build NUFFT SINGLE" ON)
#option(NUFFT_DOUBLE "Build NUFFT DOUBLE" ON)
#option(NUFFT_COMPLEX "Build NUFFT COMPLEX" ON)
#option(NUFFT_COMPLEX16 "Build NUFFT COMPLEX16" ON)

if(NUFFT_SINGLE)
  set(ALLOBJF90 ${SNUFFTF90})
  set(ALLOBJF77 ${SNUFFTF77})
endif()
if(NUFFT_DOUBLE)
  set(ALLOBJF90 ${DNUFFTF90})
  set(ALLOBJF77 ${DNUFFTF77})
endif()
if(NUFFT_COMPLEX)
  set(ALLOBJF90 ${CNUFFTF90})
  set(ALLOBJF77 ${CNUFFTF77})
endif()
if(NUFFT_COMPLEX16)
  set(ALLOBJF90 ${ZNUFFTF90})
  set(ALLOBJF77 ${ZNUFFTF77})
endif()

if(NUFFT_F90)
add_library(nufft ${ALLOBJF90})
else()
add_library(nufft ${ALLOBJF77})
endif()

if(UNIX)
  target_link_libraries(nufft m)
endif()

if(ANDROID_NDK)
  target_link_libraries(nufft m)
endif(ANDROID_NDK)

#target_link_libraries(nufft)
nufft_install_library(nufft)


############################
#option(BUILD_TESTING "Build BLAS" ON)

if(BUILD_TESTING)
add_subdirectory(TESTING)
endif(BUILD_TESTING)
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/blas.pc.in ${CMAKE_CURRENT_BINARY_DIR}/blas.pc)
#install(FILES
#  ${CMAKE_CURRENT_BINARY_DIR}/blas.pc
#  DESTINATION ${PKG_CONFIG_DIR}
#  )
