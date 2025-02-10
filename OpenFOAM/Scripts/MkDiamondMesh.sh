#! /bin/bash -e
#=============================================================================#
#                      "Scripts/RunDiamondCasesAll":                          #
#=============================================================================#
[[ "$1" == "" ]] || cd $1

# This script must be started from the directory containing the "diamondCase":
if [ ! -d diamondCase ]
then
	echo "Run this script from a directory containing the 'diamondCase' sub-dir"
	exit 1
fi

cd diamondCase
blockMesh
checkMesh -allGeometry -allRegions -allTopology

