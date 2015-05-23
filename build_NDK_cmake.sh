#!/bin/bash
#export CMAKE_BUILD_TYPE "Debug"
export CMAKE_BUILD_TYPE="Release"

if [ -z "${NDK_ROOT}"  ]; then
	export NDK_ROOT=${HOME}/NDK/android-ndk-r10d
	#export NDK_ROOT=${HOME}/NDK/android-ndk-r9
fi
export ANDROID_NDK=${NDK_ROOT}

while [ $# -ge 1 ]; do
	case $1 in
	-ABI|-abi)
		echo "\$1=-abi"
		shift
		APP_ABI=$1
		shift
		;;
	-clean|-c|-C) #
		echo "\$1=-c,-C,-clean"
		clean_build=1
		shift
		;;
	-l|-L)
		echo "\$1=-l,-L"
		local_build=1
		;;
	--help|-h|-H)
		# The main case statement will give a usage message.
		echo "$0 -c|-clean -abi=[armeabi, armeabi-v7a, armv8-64,mips,mips64el, x86,x86_64]"
		exit 1
		break
		;;
	-*)
		echo "$0: unrecognized option $1" >&2
		exit 1
		;;
	*)
		break
		;;
	esac
done
echo APP_ABI=$APP_ABI
export APP_ABI

if [[ ${NDK_ROOT} =~ .*"-r9".* ]]
then
#ANDROID_APIVER=android-8
#ANDROID_APIVER=android-9
#android 4.0.1 ICS and above
ANDROID_APIVER=android-14
#TOOL_VER="4.6"
#gfortran is in r9d V4.8.0
TOOL_VER="4.8.0"
else
#r10d : android 4.0.1 ICS and above
if [ "$APP_ABI" = "arm64-v8a" -o \
	"$APP_ABI" = "x86_64" -o \
	"$APP_ABI" = "mips64" ]; then
	ANDROID_APIVER=android-21
else
	ANDROID_APIVER=android-14
fi
TOOL_VER="4.9"
fi
#echo ANDROID_APIVER=$ANDROID_APIVER
#read
#default is arm
#export PATH="$NDK_ROOT/toolchains/${TARGPLAT}-${TOOL_VER}/prebuilt/${HOSTPLAT}/bin/:\
#$NDK_ROOT/toolchains/${TARGPLAT}-${TOOL_VER}/prebuilt/${HOSTPLAT}/${TARGPLAT}/bin/:$PATH"
case $APP_ABI in
  armeabi)
    TARGPLAT=arm-linux-androideabi
    TOOLCHAINS=arm-linux-androideabi
    ARCH=arm

	export BLIS_LIB_NAME=${BLIS_LIB_NAME:-blis-NDK-arm}
	FFTS_LIB_NAME=ffts
	#libpico-NDK-arm.a
	PICORT_LIB_NAME=picort
	#enable VFP only
  ;;
  armeabi-v7a)
    TARGPLAT=arm-linux-androideabi
    TOOLCHAINS=arm-linux-androideabi
    ARCH=arm
	#enable NEON
	export BLIS_LIB_NAME=${BLIS_LIB_NAME:-blis-NDK-arm}
	FFTS_LIB_NAME=ffts
	PICORT_LIB_NAME=picort
  ;;
  arm64-v8a)
    TARGPLAT=aarch64-linux-android
    TOOLCHAINS=aarch64-linux-android
    ARCH=arm64

	export BLIS_LIB_NAME=${BLIS_LIB_NAME:-blis-NDK-arm}
	FFTS_LIB_NAME=ffts
	PICORT_LIB_NAME=picort
  ;;
  x86)#atom-32
    TARGPLAT=i686-linux-android
    TOOLCHAINS=x86
    ARCH=x86
	#specify assembler for x86 SSE3, but ffts's sse.s needs 64bit x86.
	#intel atom z2xxx and the old atoms are 32bit
	#http://forum.cvapp.org/viewtopic.php?f=13&t=423&sid=4c47343b1de899f9e1b0d157d04d0af1
	export  CCAS="${TARGPLAT}-as"
	export  CCASFLAGS="--32 -march=i686+sse3"

	export BLIS_LIB_NAME=${BLIS_LIB_NAME:-blis-NDK-x86}
	PICORT_LIB_NAME=picort

	echo "$APP_ABI is not supported in dsplib yet!!!"
  ;;
  x86_64)
    TARGPLAT=x86_64-linux-android
    TOOLCHAINS=x86_64
    ARCH=x86_64
    #specify assembler for x86 SSE3, but ffts's sse.s needs 64bit x86.
	#atom-64 or x86-64 devices only.
	#http://forum.cvapp.org/viewtopic.php?f=13&t=423&sid=4c47343b1de899f9e1b0d157d04d0af1
	export  CCAS="${TARGPLAT}-as"
#	export  CCASFLAGS="--64 -march=i686+sse3"
	export  CCASFLAGS="--64"

	export BLIS_LIB_NAME=${BLIS_LIB_NAME:-blis-NDK-x86}
	PICORT_LIB_NAME=picort
  ;;
  mips)
	## probably wrong
	TARGPLAT=mipsel-linux-android
	TOOLCHAINS=mipsel-linux-android
	ARCH=mips
	echo "$APP_ABI is not supported in dsplib yet!!!"
  ;;
  mips64)
	## probably wrong
	TARGPLAT=mips64el-linux-android
	TOOLCHAINS=mips64el-linux-android
	ARCH=mips64
	echo "$APP_ABI is not supported in FFTS yet!!!"
  ;;
  *) echo $0: Unknown target; exit
esac

echo "Using: $NDK_ROOT/toolchains/${TOOLCHAINS}-${TOOL_VER}/prebuilt/${HOSTPLAT}/bin"
export PATH="${NDK_ROOT}/toolchains/${TOOLCHAINS}-${TOOL_VER}/prebuilt/${HOSTPLAT}/bin/:$PATH"

export SYS_ROOT="${NDK_ROOT}/platforms/${ANDROID_APIVER}/arch-${ARCH}/"
export CC="${TARGPLAT}-gcc --sysroot=$SYS_ROOT"
export LD="${TARGPLAT}-ld"
export AR="${TARGPLAT}-ar"
export RANLIB="${TARGPLAT}-ranlib"
export STRIP="${TARGPLAT}-strip"
#export CFLAGS="-Os -fPIE"
export CFLAGS="-Os -fPIE --sysroot=$SYS_ROOT"
export CXXFLAGS="-fPIE --sysroot=$SYS_ROOT"
export FORTRAN="${TARGPLAT}-gfortran --sysroot=$SYS_ROOT"

#!!! quite importnat for cmake to define the NDK's fortran compiler.!!!
#Don't let cmake decide it.
export FC=${FORTRAN}
export AM_ANDROID_EXTRA="-llog -fPIE -pie"

#include path :
#platforms/android-21/arch-arm/usr/include/

if [ -z "$PICO_DIR" ]; then
	export DSP_HOME=${DSP_HOME:-`pwd`/../../../../dsp}
	export PICO_DIR=${PICO_DIR:-`pwd`/../..}
	export BLIS_DIR=${DSP_HOME}/blis
	export BLIS_DIR=${DSP_HOME}/blis
	export BLISLIB_DIR=${BLISLIB_DIR:-${BLIS_DIR}/lib}
	export FFTS_DIR=${DSP_HOME}/ffts
	export LAPACK_SRC=${LAPACK_SRC:-${DSP_HOME}/LAPACK}
else
	export BLISLIB_DIR=${BLISLIB_DIR:-$BLIS_OUT/lib}
fi
export LAPACKE_SRC=${LAPACKE_SRC:-${LAPACK_SRC}/LAPACKE}
export CBLAS_SRC=${CBLAS_SRC:-${LAPACK_SRC}/CBLAS}
export PICORT_DIR=${PICORT_DIR:-`pwd`}

if [ -z "$PICO_OUT" ]; then
	export PICO_OUT=${PICO_OUT:-$PICO_DIR/build}
	export local_build=1
fi

#check if it needs a clean build?
if [ -d "$PICO_OUT/$APP_ABI" ]; then
	if [ -n "$clean_build" ]; then
		rm -rf $PICO_OUT/$APP_ABI/*
	fi
else
	mkdir -p $PICO_OUT/$APP_ABI
fi

pushd ${PICO_OUT}/$APP_ABI

case $APP_ABI in
	armeabi)
	cmake -DCMAKE_TOOLCHAIN_FILE=${PICORT_DIR}/android.toolchain.cmake \
	-DANDROID_NATIVE_API_LEVEL=${ANDROID_APIVER} \
	-DANDROID_NDK=${ANDROID_NDK} -DANDROID_TOOLCHAIN_NAME=${TOOLCHAINS}-${TOOL_VER} \
	-DCMAKE_BUILD_TYPE=Release -DANDROID_ABI=$APP_ABI \
	-DFFTS_DIR:FILEPATH=${FFTS_DIR} -DFFTS_LIB_NAME=${FFTS_LIB_NAME} \
	-DFFTS_OUT:FILEPATH=${FFTS_OUT} \
	-DLAPACK_SRC:FILEPATH=${LAPACK_SRC} -DLAPACKE_SRC:FILEPATH=${LAPACKE_SRC} \
	-DLAPACK_OUT:FILEPATH=${LAPACK_OUT} \
	-DUSE_BLIS=1 -DBLIS_DIR:FILEPATH=${BLIS_DIR} -DBLISLIB_DIR:FILEPATH=${BLISLIB_DIR} \
	-DBLIS_LIB_NAME=${BLIS_LIB_NAME} \
	-DAPP_ABI=${APP_ABI} \
	-DCBLAS_SRC:FILEPATH=${CBLAS_SRC} ${PICORT_DIR}
	;;
	armeabi-v7a)
	cmake -DCMAKE_TOOLCHAIN_FILE=${PICORT_DIR}/android.toolchain.cmake \
	-DANDROID_NATIVE_API_LEVEL=${ANDROID_APIVER} \
	-DANDROID_NDK=${ANDROID_NDK} -DANDROID_TOOLCHAIN_NAME=${TOOLCHAINS}-${TOOL_VER} \
	-DCMAKE_BUILD_TYPE=Release -DANDROID_ABI="armeabi-v7a with NEON" \
	-DFFTS_DIR:FILEPATH=${FFTS_DIR} -DFFTS_LIB_NAME=${FFTS_LIB_NAME} \
	-DFFTS_OUT:FILEPATH=${FFTS_OUT} \
	-DLAPACK_SRC:FILEPATH=${LAPACK_SRC} -DLAPACKE_SRC:FILEPATH=${LAPACKE_SRC} \
	-DLAPACK_OUT:FILEPATH=${LAPACK_OUT} \
	-DUSE_BLIS=1 -DBLIS_DIR:FILEPATH=${BLIS_DIR} -DBLISLIB_DIR:FILEPATH=${BLISLIB_DIR} \
	-DBLIS_LIB_NAME=${BLIS_LIB_NAME} \
	-DAPP_ABI=${APP_ABI} \
	-DCBLAS_SRC:FILEPATH=${CBLAS_SRC} ${PICORT_DIR}
	;;
	arm64-v8a)
	cmake -DCMAKE_TOOLCHAIN_FILE=${PICORT_DIR}/android.toolchain.cmake \
	-DANDROID_NATIVE_API_LEVEL=${ANDROID_APIVER} \
	-DANDROID_NDK=${ANDROID_NDK} -DANDROID_TOOLCHAIN_NAME=${TOOLCHAINS}-${TOOL_VER} \
	-DCMAKE_BUILD_TYPE=Release -DANDROID_ABI=$APP_ABI \
	-DFFTS_DIR:FILEPATH=${FFTS_DIR} -DFFTS_LIB_NAME=${FFTS_LIB_NAME} \
	-DFFTS_OUT:FILEPATH=${FFTS_OUT} \
	-DLAPACK_SRC:FILEPATH=${LAPACK_SRC} -DLAPACKE_SRC:FILEPATH=${LAPACKE_SRC} \
	-DLAPACK_OUT:FILEPATH=${LAPACK_OUT} \
	-DUSE_BLIS=1 -DBLIS_DIR:FILEPATH=${BLIS_DIR} -DBLISLIB_DIR:FILEPATH=${BLISLIB_DIR} \
	-DBLIS_LIB_NAME=${BLIS_LIB_NAME} \
	-DAPP_ABI=${APP_ABI} \
	-DCBLAS_SRC:FILEPATH=${CBLAS_SRC} ${PICORT_DIR}
	;;
  x86)
	cmake -DCMAKE_TOOLCHAIN_FILE=${PICORT_DIR}/android.toolchain.cmake \
	-DANDROID_NATIVE_API_LEVEL=${ANDROID_APIVER} \
	-DANDROID_NDK=${NDK_ROOT} -DANDROID_TOOLCHAIN_NAME=${TOOLCHAINS}-${TOOL_VER} \
	-DCMAKE_BUILD_TYPE=Release -DANDROID_ABI=$APP_ABI \
	-DFFTS_DIR:FILEPATH=${FFTS_DIR} -DFFTS_LIB_NAME=${FFTS_LIB_NAME} \
	-DFFTS_OUT:FILEPATH=${FFTS_OUT} \
	-DLAPACK_SRC:FILEPATH=${LAPACK_SRC} -DLAPACKE_SRC:FILEPATH=${LAPACKE_SRC} \
	-DLAPACK_OUT:FILEPATH=${LAPACK_OUT} \
	-DUSE_BLIS=1 -DBLIS_DIR:FILEPATH=${BLIS_DIR} -DBLISLIB_DIR:FILEPATH=${BLISLIB_DIR} \
	-DBLIS_LIB_NAME=${BLIS_LIB_NAME} \
	-DAPP_ABI=${APP_ABI} \
	-DCBLAS_SRC:FILEPATH=${CBLAS_SRC} ${PICORT_DIR}
  ;;
  x86_64)
	cmake -DCMAKE_TOOLCHAIN_FILE=${PICORT_DIR}/android.toolchain.cmake \
	-DANDROID_NATIVE_API_LEVEL=${ANDROID_APIVER} \
	-DANDROID_NDK=${NDK_ROOT} -DANDROID_TOOLCHAIN_NAME=${TOOLCHAINS}-${TOOL_VER} \
	-DCMAKE_BUILD_TYPE=Release -DANDROID_ABI=$APP_ABI \
	-DFFTS_DIR:FILEPATH=${FFTS_DIR} -DFFTS_LIB_NAME=${FFTS_LIB_NAME} \
	-DFFTS_OUT:FILEPATH=${FFTS_OUT} \
	-DLAPACK_SRC:FILEPATH=${LAPACK_SRC} -DLAPACKE_SRC:FILEPATH=${LAPACKE_SRC} \
	-DLAPACK_OUT:FILEPATH=${LAPACK_OUT} \
	-DUSE_BLIS=1 -DBLIS_DIR:FILEPATH=${BLIS_DIR} -DBLISLIB_DIR:FILEPATH=${BLISLIB_DIR} \
	-DBLIS_LIB_NAME=${BLIS_LIB_NAME} \
	-DAPP_ABI=${APP_ABI} \
	-DCBLAS_SRC:FILEPATH=${CBLAS_SRC} ${PICORT_DIR}
  ;;
  mips)
	cmake -DCMAKE_TOOLCHAIN_FILE=${PICORT_DIR}/android.toolchain.cmake \
	-DANDROID_NATIVE_API_LEVEL=${ANDROID_APIVER} \
	-DANDROID_NDK=${NDK_ROOT} -DANDROID_TOOLCHAIN_NAME=${TOOLCHAINS}-${TOOL_VER} \
	-DCMAKE_BUILD_TYPE=Release -DANDROID_ABI=$APP_ABI \
	-DFFTS_DIR:FILEPATH=${FFTS_DIR} -DFFTS_LIB_NAME=${FFTS_LIB_NAME} \
	-DFFTS_OUT:FILEPATH=${FFTS_OUT} \
	-DLAPACK_SRC:FILEPATH=${LAPACK_SRC} -DLAPACKE_SRC:FILEPATH=${LAPACKE_SRC} \
	-DLAPACK_OUT:FILEPATH=${LAPACK_OUT} \
	-DUSE_BLIS=1 -DBLIS_DIR:FILEPATH=${BLIS_DIR} -DBLISLIB_DIR:FILEPATH=${BLISLIB_DIR} \
	-DBLIS_LIB_NAME=${BLIS_LIB_NAME} \
	-DAPP_ABI=${APP_ABI} \
	-DCBLAS_SRC:FILEPATH=${CBLAS_SRC} ${PICORT_DIR}
  ;;
  mips64)
	cmake -DCMAKE_TOOLCHAIN_FILE=${PICORT_DIR}/android.toolchain.cmake \
	-DANDROID_NATIVE_API_LEVEL=${ANDROID_APIVER} \
	-DANDROID_NDK=${NDK_ROOT} -DANDROID_TOOLCHAIN_NAME=${TOOLCHAINS}-${TOOL_VER} \
	-DCMAKE_BUILD_TYPE=Release -DANDROID_ABI=$APP_ABI \
	-DFFTS_DIR:FILEPATH=${FFTS_DIR} -DFFTS_LIB_NAME=${FFTS_LIB_NAME} \
	-DFFTS_OUT:FILEPATH=${FFTS_OUT} \
	-DLAPACK_SRC:FILEPATH=${LAPACK_SRC} -DLAPACKE_SRC:FILEPATH=${LAPACKE_SRC} \
	-DLAPACK_OUT:FILEPATH=${LAPACK_OUT} \
	-DUSE_BLIS=1 -DBLIS_DIR:FILEPATH=${BLIS_DIR} -DBLISLIB_DIR:FILEPATH=${BLISLIB_DIR} \
	-DBLIS_LIB_NAME=${BLIS_LIB_NAME} \
	-DAPP_ABI=${APP_ABI} \
	-DCBLAS_SRC:FILEPATH=${CBLAS_SRC} ${PICORT_DIR}
  ;;
  *) echo $0: Unknown target; exit
esac

ret=$?
echo "ret=$ret"
if [ "$ret" != '0' ]; then
echo "$0 cmake error!!!!"
exit -1
fi

make -j${CORE_COUNT}

ret=$?
echo "ret=$ret"
if [ "$ret" != '0' ]; then
echo "$0 make error!!!!"
exit -1
fi

popd

pushd ${PICO_OUT}
mkdir -p libs/$APP_ABI
#cp ${PICORT_DIR}/Android.mk ${PICO_OUT}/libs
rm -rf libs/$APP_ABI/*
ln -s ${PICO_OUT}/$APP_ABI/lib/libpicort.a libs/$APP_ABI/libpicort.a
popd

exit 0
