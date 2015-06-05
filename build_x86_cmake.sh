#!/bin/bash
# Make sure you have NDK_ROOT defined in .bashrc or .bash_profile

#export CMAKE_BUILD_TYPE "Debug"
export CMAKE_BUILD_TYPE="Release"

#get cpu counts
case $(uname -s) in
  Darwin)
    CONFBUILD=i386-apple-darwin`uname -r`
    HOSTPLAT=darwin-x86
    CORE_COUNT=`sysctl -n hw.ncpu`
  ;;
  Linux)
    CONFBUILD=x86-unknown-linux
    HOSTPLAT=linux-`uname -m`
    CORE_COUNT=`grep processor /proc/cpuinfo | wc -l`
  ;;
CYGWIN*)
	CORE_COUNT=`grep processor /proc/cpuinfo | wc -l`
	;;
  *) echo $0: Unknown platform; exit
esac

if [ -z "${NUFFT_SRC}" ];then
export NUFFT_SRC=`pwd`
fi

BUILD_NUFFT_GENERIC=1
BUILD_NUFFT_VEC=0
BUILD_NUFFT_NEON=0
BUILD_NUFFT_CUDA=0

while [ $# -ge 1 ]; do
	case $1 in
	-acc|-ACC)
		shift
		case $1 in
		GENERIC|generic)
			BUILD_NUFFT_GENERIC=1
			;;
		VEC|vec)
			BUILD_NUFFT_VEC=1
			;;

		NEON|neon)
			BUILD_NUFFT_NEON=1
			;;

		CUDA|cuda)
			BUILD_NUFFT_CUDA=1
			;;

		*)
			BUILD_NUFFT_GENERIC=1
			;;
		esac
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
if [ -z "$NUFFT_OUT" ]; then
	export NUFFT_DIR=`pwd`
	export NUFFT_OUT=${NUFFT_OUT:-$NUFFT_DIR/build}
fi

#check if it needs a clean build?
if [ -d "$NUFFT_OUT/$TARGET_ARCH" ]; then
	if [ -n "$clean_build" ]; then
		rm -rf $NUFFT_OUT/$TARGET_ARCH/*
	fi
else
	mkdir -p $NUFFT_OUT/$TARGET_ARCH
fi

pushd ${NUFFT_OUT}/$TARGET_ARCH

cmake ${NUFFT_DIR}
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
pushd ${NUFFT_OUT}

mkdir -p libs/$TARGET_ARCH
rm -rf libs/$TARGET_ARCH/*

ln -s ${NUFFT_OUT}/$TARGET_ARCH/lib/libnufft.a ${NUFFT_OUT}/libs/$TARGET_ARCH/libnufft.a
ln -s ${NUFFT_OUT}/$TARGET_ARCH/lib/libfftbase.a ${NUFFT_OUT}/libs/$TARGET_ARCH/libfftbase.a

popd
exit 0
