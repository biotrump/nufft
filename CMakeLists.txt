cmake_minimum_required(VERSION 2.8.10)
project(NUFFT Fortran)

set(NUFFT_VERSION 1.3.2)

# Configure the warning and code coverage suppression file
configure_file(
  "${NUFFT_SOURCE_DIR}/CTestCustom.cmake.in"
  "${NUFFT_BINARY_DIR}/CTestConfig.cmake"
  COPYONLY
)

# Add the CMake directory for custon CMake modules
set(CMAKE_MODULE_PATH "${NUFFT_SOURCE_DIR}/CMAKE" ${CMAKE_MODULE_PATH})

#F95, F90, F77
option(NUFFT_F90 "Build NUFFT by fortran F90" ON)

if(ANDROID_NDK)
	# export FC=/path/to/android-ndk-r9/toolchains/arm-linux-androideabi-4.8.0/prebuilt/linux-x86/bin/arm-linux-androideabi-gfortran
	message ("**** CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
	get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
	if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
	# gfortran
	message ("*** Fortran compiler: " ${Fortran_COMPILER_NAME})
#KitKat4.4.x and Lollip 5.x :executable
#	set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fPIE -pie -fPIC -fno-f2c -O3")
#	set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -fPIE -pie -fPIC -g")
#ICS4.0.x :executable
	set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fPIC -fno-f2c -O3")
	set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -fPIC -g")
	else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
	message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
	message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
	message ("No optimized Fortran compiler flags are known, we just try -O2...")
	set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
	set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
	endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")
else(ANDROID_NDK)
	if (UNIX)
	if ( "${CMAKE_Fortran_COMPILER}" MATCHES "ifort" )
	set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fltconsistency -fp_port" )
	endif ()
	if ( "${CMAKE_Fortran_COMPILER}" MATCHES "xlf" )
	set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qnosave -qstrict=none" )
	endif ()
	# Delete libmtsk in linking sequence for Sun/Oracle Fortran Compiler.
	# This library is not present in the Sun package SolarisStudio12.3-linux-x86-bin
	STRING(REPLACE \;mtsk\; \; CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")
	endif ()
	message ("*** Fortran compiler RELEASE FLAGS: " ${CMAKE_Fortran_FLAGS_RELEASE})
	message ("*** Fortran compiler DEBUG FLAGS: " ${CMAKE_Fortran_FLAGS_DEBUG})
endif(ANDROID_NDK)

# Get Python
IF(WIN32 AND NOT UNIX)
find_package(PythonInterp)
message(STATUS "Looking for Python found - ${PYTHONINTERP_FOUND}")
if (PYTHONINTERP_FOUND)
   message(STATUS "Using Python version ${PYTHON_VERSION_STRING}")
endif()
ENDIF(WIN32 AND NOT UNIX)

# --------------------------------------------------
#install all .a, .so to nufft/lib
set(NUFFT_INSTALL_EXPORT_NAME nufft-targets)

macro(nufft_install_library lib)
  install(TARGETS ${lib}
    EXPORT ${NUFFT_INSTALL_EXPORT_NAME}
    ARCHIVE DESTINATION lib${LIB_SUFFIX}
    LIBRARY DESTINATION lib${LIB_SUFFIX}
    RUNTIME DESTINATION bin
  )
endmacro()

# --------------------------------------------------
# Testing

enable_testing()
#BUILD_TESTING is turned ON by enable_testing()
include(CTest)
enable_testing()
# --------------------------------------------------

# Organize output files.  On Windows this also keeps .dll files next
# to the .exe files that need them, making tests easy to run.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${NUFFT_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${NUFFT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${NUFFT_BINARY_DIR}/lib)

# --------------------------------------------------
# Check for any necessary platform specific compiler flags
include( CheckNUFFTCompilerFlags )
CheckNUFFTCompilerFlags()

# --------------------------------------------------
# Precision to build
# By default all precisions are generated
if(ANDROID_NDK)
	set( CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER )
	set( CMAKE_FIND_ROOT_PATH_MODE_LIBRARY NEVER )
	set( CMAKE_FIND_ROOT_PATH_MODE_INCLUDE NEVER )
endif(ANDROID_NDK)
# --------------------------------------------------
# Subdirectories that need to be processed

# --------------------------------------------------
# FFTE
#option(FFTE "Build FFTE" OFF)

#if(FFTE)
#  add_subdirectory(ffte)
#endif(FFTE)

# --------------------------------------------------
#find ffte
message("FFTE_DIR=${FFTE_DIR}")
find_path(FFTE_INC ffte.h
			HINTS ${FFTE_DIR}/include )
message("FFTE_INC=${FFTE_INC}" )
if(FFTE_INC)
	include_directories("${FFTE_INC}")
endif(FFTE_INC)

message("FFTE_OUT=${FFTE_OUT} ")
if(ANDROID_NDK)
find_library (FFTELIB
			NAMES ${FFTE_LIB_NAME}
			PATHS ${FFTE_OUT}
			PATH_SUFFIXES lib
			NO_DEFAULT_PATH)
else()
find_library (FFTELIB
			NAMES ${FFTE_LIB_NAME}
			PATHS ${FFTE_OUT}
			PATH_SUFFIXES lib
			NO_DEFAULT_PATH)
endif()
message("**FFTELIB=${FFTELIB} ")
# --------------------------------------------------
if(ANDROID_NDK)
	set( CMAKE_FIND_ROOT_PATH_MODE_PROGRAM ONLY )
	set( CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY )
	set( CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY )
endif(ANDROID_NDK)
# NUFFT
option(USE_OPTIMIZED_FFT "Whether or not to use an optimized fft library instead of included netlib fftpack" OFF)

# Check the usage of the user provided BLAS libraries
if(FFT_LIBRARIES)
  include(CheckFortranFunctionExists)
  set(CMAKE_REQUIRED_LIBRARIES ${FFT_LIBRARIES})
  CHECK_FORTRAN_FUNCTION_EXISTS("dgemm" FFT_FOUND)
  unset( CMAKE_REQUIRED_LIBRARIES )
  if(FFT_FOUND)
    message(STATUS "--> BLAS supplied by user is WORKING, will use ${FFT_LIBRARIES}.")
  else(FFT_FOUND)
    message(ERROR "--> BLAS supplied by user is not WORKING, CANNOT USE ${FFT_LIBRARIES}.")
    message(ERROR "-->     Will use REFERENCE BLAS (by default)")
    message(ERROR "-->     Or Correct your FFT_LIBRARIES entry ")
    message(ERROR "-->     Or Consider checking USE_OPTIMIZED_FFT")
  endif(FFT_FOUND)

# User did not provide a BLAS Library but specified to search for one
elseif( USE_OPTIMIZED_FFT )
  find_package( FFT
  PATHS ../ffte)
endif (FFT_LIBRARIES)

# Neither user specified or optimized BLAS libraries can be used
if(NOT FFT_FOUND)
  message(STATUS "Using supplied NETLIB fftpack implementation")
  add_subdirectory(src/fftpack)
  set( FFT_LIBRARIES fftbase )
else()
  set( CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}"
    CACHE STRING "Linker flags for executables" FORCE)
  set( CMAKE_MODULE_LINKER_FLAGS
    "${CMAKE_MODULE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}"
    CACHE STRING "Linker flags for modules" FORCE)
  set( CMAKE_SHARED_LINKER_FLAGS
    "${CMAKE_SHARED_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}"
    CACHE STRING "Linker flags for shared libs" FORCE)
endif( NOT FFT_FOUND )

# --------------------------------------------------
# CNUFFT
option(CNUFFT "Build CNUFFT" OFF)

if(CNUFFT)
  add_subdirectory(CNUFFT)
endif(CNUFFT)

# --------------------------------------------------
# NUFFT
# User did not provide a LAPACK Library but specified to search for one
if( USE_OPTIMIZED_NUFFT )
  find_package( NUFFT )
endif (USE_OPTIMIZED_NUFFT)

# Check the usage of the user provided or automatically found NUFFT libraries
if(NUFFT_LIBRARIES)
  include(CheckFortranFunctionExists)
  set(CMAKE_REQUIRED_LIBRARIES ${NUFFT_LIBRARIES})
  # Check if new routine of 3.4.0 is in NUFFT_LIBRARIES
  CHECK_FORTRAN_FUNCTION_EXISTS("next235" LATESTNUFFT_FOUND)
  unset( CMAKE_REQUIRED_LIBRARIES )
  if(LATESTNUFFT_FOUND)
    message(STATUS "--> NUFFT supplied by user is WORKING, will use ${NUFFT_LIBRARIES}.")
  else(NUFFT_FOUND)
    message(ERROR "--> NUFFT supplied by user is not WORKING or is older than NUFFT 1.3.2, 	CANNOT USE ${NUFFT_LIBRARIES}.")
    message(ERROR "-->     Will use REFERENCE NUFFT (by default)")
    message(ERROR "-->     Or Correct your NUFFT_LIBRARIES entry ")
    message(ERROR "-->     Or Consider checking USE_OPTIMIZED_NUFFT")
  endif(LATESTNUFFT_FOUND)
endif (NUFFT_LIBRARIES)

# Neither user specified or optimized NUFFT libraries can be used
if(NOT LATESTNUFFT_FOUND)
  message(STATUS "Using supplied NUFFT implementation")
  set( NUFFT_LIBRARIES nufft )
  option(BUILD_SINGLE "Build NUFFT Single Precision" OFF)
  option(BUILD_DOUBLE "Build NUFFT Double Precision" OFF)
  option(BUILD_COMPLEX "Build NUFFT Complex Precision" OFF)
  option(BUILD_COMPLEX16 "Build NUFFT Double Complex Precision" ON)
  add_subdirectory(src)
else()
  set( CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} ${NUFFT_LINKER_FLAGS}"
    CACHE STRING "Linker flags for executables" FORCE)
  set( CMAKE_MODULE_LINKER_FLAGS
    "${CMAKE_MODULE_LINKER_FLAGS} ${NUFFT_LINKER_FLAGS}"
    CACHE STRING "Linker flags for modules" FORCE)
  set( CMAKE_SHARED_LINKER_FLAGS
    "${CMAKE_SHARED_LINKER_FLAGS} ${NUFFT_LINKER_FLAGS}"
    CACHE STRING "Linker flags for shared libs" FORCE)
endif( NOT LATESTNUFFT_FOUND )

message(STATUS "BUILD TESTING : ${BUILD_TESTING}" )
if(BUILD_TESTING)
#  add_subdirectory(TESTING)
endif(BUILD_TESTING)

# --------------------------------------------------
# CPACK Packaging

SET(CPACK_PACKAGE_NAME "NUFFT")
SET(CPACK_PACKAGE_VENDOR "The Courant Mathematics and Computing Laboratory (CMCL) is a research center \
						at the Courant Institute of New York University.")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "NUFFT- Linear Algebra Package")
set(CPACK_PACKAGE_VERSION_MAJOR 1)
set(CPACK_PACKAGE_VERSION_MINOR 3)
set(CPACK_PACKAGE_VERSION_PATCH 2)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/license.txt")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "NUFFT")
IF(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
  SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\http://www.cims.nyu.edu/cmcl/nufft/nufft.html")
  SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.cims.nyu.edu/cmcl/nufft/nufft.html")
  SET(CPACK_NSIS_CONTACT "thomas@biotrump.com")
  SET(CPACK_NSIS_MODIFY_PATH ON)
  SET(CPACK_NSIS_DISPLAY_NAME "NUFFT-${NUFFT_VERSION}")
  set(CPACK_PACKAGE_RELOCATABLE "true")
ELSE(WIN32 AND NOT UNIX)
  SET(CPACK_GENERATOR "TGZ")
  SET(CPACK_SOURCE_GENERATOR TGZ)
  SET(CPACK_SOURCE_PACKAGE_FILE_NAME "nufft-${NUFFT_VERSION}" )
  SET(CPACK_SOURCE_IGNORE_FILES ~$ .gitignore ${CPACK_SOURCE_IGNORE_FILES} )
ENDIF(WIN32 AND NOT UNIX)
INCLUDE(CPack)

# --------------------------------------------------
# By default static library
OPTION(BUILD_SHARED_LIBS "Build shared libraries" OFF )
OPTION(BUILD_STATIC_LIBS "Build static libraries" ON )
#OPTION(BUILD_SHARED_LIBS "Build shared libraries" ON )

if(NOT FFT_FOUND)
  set(ALL_TARGETS ${ALL_TARGETS} fftbase)
endif(NOT FFT_FOUND)

if(NOT LATESTNUFFT_FOUND)
  set(ALL_TARGETS ${ALL_TARGETS} nufft)
endif(NOT LATESTNUFFT_FOUND)

# Export nufft targets, not including nuffte, from the
# install tree, if any.
set(_nufft_config_install_guard_target "")
if(ALL_TARGETS)
  install(EXPORT nufft-targets
    DESTINATION lib/cmake/nufft-${NUFFT_VERSION})

  # Choose one of the nufft targets to use as a guard for
  # nufft-config.cmake to load targets from the install tree.
  list(GET ALL_TARGETS 0 _nufft_config_install_guard_target)
endif()

# Include cnufft in targets exported from the build tree.
if(CNUFFT)
  set(ALL_TARGETS ${ALL_TARGETS} cnufft)
endif(CNUFFT)

# Export nufft and nuffte targets from the build tree, if any.
set(_nufft_config_build_guard_target "")
if(ALL_TARGETS)
  export(TARGETS ${ALL_TARGETS} FILE nufft-targets.cmake)

  # Choose one of the nufft or nuffte targets to use as a guard
  # for nufft-config.cmake to load targets from the build tree.
  list(GET ALL_TARGETS 0 _nufft_config_build_guard_target)
endif()

configure_file(${NUFFT_SOURCE_DIR}/CMAKE/nufft-config-version.cmake.in
  ${NUFFT_BINARY_DIR}/nufft-config-version.cmake @ONLY)
configure_file(${NUFFT_SOURCE_DIR}/CMAKE/nufft-config-build.cmake.in
  ${NUFFT_BINARY_DIR}/nufft-config.cmake @ONLY)


#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/nufft.pc.in ${CMAKE_CURRENT_BINARY_DIR}/nufft.pc)
#  install(FILES
#  ${CMAKE_CURRENT_BINARY_DIR}/nufft.pc
#  DESTINATION ${PKG_CONFIG_DIR}
#   )

configure_file(${NUFFT_SOURCE_DIR}/CMAKE/nufft-config-install.cmake.in
  ${NUFFT_BINARY_DIR}/CMakeFiles/nufft-config.cmake @ONLY)
install(FILES
  ${NUFFT_BINARY_DIR}/CMakeFiles/nufft-config.cmake
  ${NUFFT_BINARY_DIR}/nufft-config-version.cmake
  DESTINATION lib/cmake/nufft-${NUFFT_VERSION}
  )
