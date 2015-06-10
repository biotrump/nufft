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
#libnufft.a : nufft
#libfftbase.a : Netlib fftpack
macro(add_nufft_test name src)
#  get_filename_component(baseNAME ${src} NAME_WE)
  set(baseNAME "demo")
  set(TEST_INPUT "${NUFFT_SOURCE_DIR}/${baseNAME}.in")
  add_executable(${name} ${src})
  target_link_libraries(${name} nufft fftbase)
#  if(EXISTS "${TEST_INPUT}")
#    add_test(NAME NUFFT-${name} COMMAND "${CMAKE_COMMAND}"
#      -DTEST=$<TARGET_FILE:${name}>
#      -DINPUT=${TEST_INPUT}
#      -DINTDIR=${CMAKE_CFG_INTDIR}
#      -P "${NUFFT_SOURCE_DIR}/TESTING/runtest.cmake")
#    else()
#      add_test(NAME NUFFT-${name} COMMAND "${CMAKE_COMMAND}"
#        -DTEST=$<TARGET_FILE:${name}>
#        -DINTDIR=${CMAKE_CFG_INTDIR}
#        -P "${NUFFT_SOURCE_DIR}/TESTING/runtest.cmake")
#    endif()
endmacro(add_nufft_test)

#option(BUILD_SINGLE "Build NUFFT SINGLE" ON)
#option(BUILD_DOUBLE "Build NUFFT DOUBLE" ON)
#option(BUILD_COMPLEX "Build NUFFT COMPLEX" ON)
#option(BUILD_COMPLEX16 "Build NUFFT COMPLEX16" ON)

if(BUILD_SINGLE)
if(NUFFT_F90)
  add_nufft_test(nufft1d_demof90 nufft1d_demof90.f nufft90)
  add_nufft_test(nufft2d_demof90 nufft2d_demof90.f nufft90)
  add_nufft_test(nufft3d_demof90 nufft3d_demof90.f nufft90)
else(NUFFT_F90)
  add_nufft_test(nufft1d_demof77 nufft1d_demof77.f nufft77)
  add_nufft_test(nufft2d_demof77 nufft2d_demof77.f nufft77)
  add_nufft_test(nufft3d_demof77 nufft3d_demof77.f nufft77)
endif(NUFFT_F90)
endif()

if(BUILD_DOUBLE)
if(NUFFT_F90)
  add_nufft_test(nufft1d_demof90 nufft1d_demof90.f nufft90)
  add_nufft_test(nufft2d_demof90 nufft2d_demof90.f nufft90)
  add_nufft_test(nufft3d_demof90 nufft3d_demof90.f nufft90)
else(NUFFT_F90)
  add_nufft_test(nufft1d_demof77 nufft1d_demof77.f nufft77)
  add_nufft_test(nufft2d_demof77 nufft2d_demof77.f nufft77)
  add_nufft_test(nufft3d_demof77 nufft3d_demof77.f nufft77)
endif(NUFFT_F90)
endif()

if(BUILD_COMPLEX)
if(NUFFT_F90)
  add_nufft_test(nufft1d_demof90 nufft1d_demof90.f)
  add_nufft_test(nufft2d_demof90 nufft2d_demof90.f)
  add_nufft_test(nufft3d_demof90 nufft3d_demof90.f)
else(NUFFT_F90)
  add_nufft_test(nufft1d_demof77 nufft1d_demof77.f)
  add_nufft_test(nufft2d_demof77 nufft2d_demof77.f)
  add_nufft_test(nufft3d_demof77 nufft3d_demof77.f)
endif(NUFFT_F90)
endif()

if(BUILD_COMPLEX16)
if(NUFFT_F90)
  add_nufft_test(nufft1d_demof90 nufft1d_demof90.f)
  add_nufft_test(nufft2d_demof90 nufft2d_demof90.f)
  add_nufft_test(nufft3d_demof90 nufft3d_demof90.f)
else(NUFFT_F90)
  add_nufft_test(nufft1d_demof77 nufft1d_demof77.f)
  add_nufft_test(nufft2d_demof77 nufft2d_demof77.f)
  add_nufft_test(nufft3d_demof77 nufft3d_demof77.f)
endif(NUFFT_F90)
endif()

#add_nufft_test(nufft1d_demof90-ffte nufft1d_demof90.f)
add_executable(nufft1d_demof90-ffte nufft1d_demof90-ffte.f)
target_link_libraries(nufft1d_demof90-ffte
	nufft
	${FFTELIB} )
