#!/bin/bash
#fortran compiler is not supported in all platorms
if [ $TARGET_ARCH == "all" ]; then
	./build_NDK_cmake.sh -abi armeabi -c
	if [ "$?" != "0" ]; then
		exit -1
	fi

	./build_NDK_cmake.sh -abi armeabi-v7a -c
	if [ "$?" != "0" ]; then
		exit -1
	fi
#./build_NDK_cmake.sh -abi arm64-v8a -c
#	if [ "$?" != "0"]; then
#		exit -1
#	fi

	./build_NDK_cmake.sh -abi x86 -c
	if [ "$?" != "0" ]; then
		exit -1
	fi
#./build_NDK_cmake.sh -abi x86_64 -c


#./build_NDK_cmake.sh -abi mips -c

#./build_NDK_cmake.sh -abi mips64 -c
else
	./build_NDK_cmake.sh -abi $TARGET_ARCH -c
	ret1=$?
	#echo "ret1=$ret1"
	#read
	if [ "$ret1" != "0" ]; then
		exit -1
	fi
fi
#echo "00000"
exit 0
