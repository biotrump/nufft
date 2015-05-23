#!/bin/bash

case `uname` in
"Darwin")
	# Should also work on other BSDs
	CORE_COUNT=`sysctl -n hw.ncpu`
	;;
"Linux")
	CORE_COUNT=`grep processor /proc/cpuinfo | wc -l`
	;;
CYGWIN*)
	CORE_COUNT=`grep processor /proc/cpuinfo | wc -l`
	;;
*)
	echo Unsupported platform: `uname`
	exit -1
esac

if [ ! -d build_x86 ]; then
mkdir -p build_x86
else
rm -rf build_x86/*
fi
pushd build_x86
cmake ..
make -j${CORE_COUNT}
popd
