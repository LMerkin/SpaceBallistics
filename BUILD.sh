#! /bin/bash -e
# vim:ts=2:et
#-----------------------------------------------------------------------------#
# Paths:                                                                      #
#-----------------------------------------------------------------------------#
AbsPath0=$(realpath $0)
TopDir=$(dirname $AbsPath0)
Project=$(basename $TopDir)

#-----------------------------------------------------------------------------#
# Command-Line Params:                                                        #
#-----------------------------------------------------------------------------#
ConfigMake=0
CleanUp=0
Build=0
Verbose=0
DebugMode=0
ReleaseMode=0
UnCheckedMode=0
ToolChain="GCC"
Jobs=$(nproc)
CountSrcLines=0

function usage
{
  echo "ERROR: Invalid option: $1"
  echo "Available options:"
  echo "-t ToolChain, currently supported: GCC (default), CLang"
  echo "-c       : Configure"
  echo "-C       : as above, but with full clean-up"
  echo "-b       : Build"
  echo "-j Jobs  : Number of concurrent jobs in make (default:  auto)"
  echo "-d: (use with -c|-C): Configure in the Debug     mode"
  echo "-r: (use with -c|-C): Configure in the Release   mode (default)"
  echo "-u: (use with -c|-C): Configure in the UnChecked mode"
  echo "-v: Verbose output"
  echo "-q: Quiet mode (no stdout output)"
  echo "-l: Count the Eligible Source Code Lines"
  exit 1
}

while getopts ":t:j:cCbdruvlq" opt
do
  case $opt in
    t) ToolChain="$OPTARG";;
    c) ConfigMake=1;;
    C) ConfigMake=1;    CleanUp=1;;
    b) Build=1;;
    d) DebugMode=1;;
    r) ReleaseMode=1;;
    u) UnCheckedMode=1;;
    q) Verbose=0;;
    v) Verbose=1;;
    j) Jobs="$OPTARG";;
    l) CountSrcLines=1;;
    *) usage $OPTARG;;
  esac
done

#-----------------------------------------------------------------------------#
# Set the BuildType and BuildMode:                                            #
#-----------------------------------------------------------------------------#
# NB:
# (*) UnCheckedMode is the highest-optimisation mode (higher than ReleaseMode);
# (*) UnCheckedMode cannot be combined with DebugMode;
# (*) Release and Debug Modes can be combined, but if the DebugMode is set along
#     with ReleaseMode, then UnCheckedMode cannot be set (see above);
# (*) BuildType and UnCheckedMode defs are passed to CMake, whereas BuildMode is
#     only for creating dirs (and then passed to CMake via {BIN,LIB}_DIR):
#
if [ $UnCheckedMode -eq 1 ]
then
  if [ $DebugMode -eq 1 ]
  then
    echo "UnCheckedMode is incompatible with DebugMode"
    exit 1
  fi
  # BuildType is "Release" but BuildMode is "UnChecked":
  BuildType="Release"
  BuildMode="UnChecked"
else
  # Set the BuildType. If it is not specified at all, "Release" is assumed:
  BuildType="Release"
  [ $ReleaseMode -eq 0 -a $DebugMode -eq 1 ] && BuildType="Debug"
  [ $ReleaseMode -eq 1 -a $DebugMode -eq 1 ] && BuildType="RelWithDebInfo"
  # Here BuildMode is the same as BuildType:
  BuildMode="$BuildType"
fi

#-----------------------------------------------------------------------------#
# Verify the ToolChain and set the C and C++ Compilers:                       #
#-----------------------------------------------------------------------------#
# XXX: We use the compilers of the appropriate ToolChain just as they appear in
# the current PATH:
#
case "$ToolChain" in
  "GCC")
    CXX=$(which g++)
    ;;
  "CLang")
    CXX=$(which clang++)
    ;;
  *) echo "ERROR: Invalid ToolChain=$ToolChain (must be: GCC|CLang)";
     exit 1
esac

#-----------------------------------------------------------------------------#
# Go Ahead:                                                                   #
#-----------------------------------------------------------------------------#
BldTop="$TopDir/__BUILD__/$ToolChain-$BuildMode"
BldDir="$BldTop/build"
BinDir="$BldTop/bin"
LibDir="$BldTop/lib"

# Create dirs if they don't exist:
mkdir -p $BldDir
mkdir -p $BinDir
mkdir -p $LibDir

#-----------------------------------------------------------------------------#
# Count the Lines in Non-Auto-Generated Src Files:                            #
#-----------------------------------------------------------------------------#
if [ $CountSrcLines -eq 1 ]
then
  EligFiles=$(find  . -type f              | \
              grep -v Docs                 | \
              grep -v .git                 | \
              grep -v __BUILD__            | \
              grep -v DE440T-1650-1750.hpp | \
              grep -v DE440T-1750-1850.hpp | \
              grep -v DE440T-1850-1950.hpp | \
              grep -v DE440T-1950-2050.hpp | \
              grep -v DE440T-2050-2150.hpp | \
              grep -v GravityFld-Earth.hpp | \
              grep -v GravityFld-Moon.hpp )
  echo $(wc -l $EligFiles | grep total)
fi

#-----------------------------------------------------------------------------#
# Configure:                                                                  #
#-----------------------------------------------------------------------------#
# Generate Makefiles if requested or build directory is empty:
if [ $ConfigMake -eq 1 ] || [ ! "$(ls -A $BldDir)" ]
then
  # Remove all files in {Bld,Bin,Lib}Dir if requested:
  [ $CleanUp -eq 1 ] && rm -fr $BldDir/*
  [ $CleanUp -eq 1 ] && rm -fr $BinDir/*
  [ $CleanUp -eq 1 ] && rm -fr $LibDir/*

  echo "Generating files in $BldDir..."

  # Run CMake:
  # {CMAKE,PROJECT}_SOURCE_DIR will be "$TopDir" (passed via the -S arg),
  # {CMAKE,PROJECT}_BINARY_DIR will be "$BldDir" (passed via the -B arg):
  cmake \
    -G "Unix Makefiles"   \
    -D CMAKE_CXX_COMPILER="$CXX" \
    -D TOOL_CHAIN="$ToolChain"   \
    -D CMAKE_BUILD_TYPE="$BuildType"      \
    -D UNCHECKED_MODE="$UnCheckedMode"    \
    -D PROJ_NAME="$Project"    \
    -D LIB_DIR="$LibDir"  \
    -D BIN_DIR="$BinDir"  \
    -S "$TopDir" \
    -B "$BldDir"
fi

#-----------------------------------------------------------------------------#
# Build:                                                                      #
#-----------------------------------------------------------------------------#
# NB: Only Makefile-based build can be done in batch mode, so the postfix after
# "BinBase" is always empty:
#
if [ $Build -eq 1 ]
then
  if [ $Verbose -eq 1 ]; then MVerbose="VERBOSE=1"; else MVerbose=""; fi

  cmake --build $BldDir -- -j $Jobs $MVerbose
fi
exit 0
