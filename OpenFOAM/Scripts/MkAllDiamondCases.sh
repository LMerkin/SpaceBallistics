#! /bin/bash -e
#=============================================================================#
#                      "Scripts/RunDiamondCasesAll":                          #
#=============================================================================#
AbsPath0=$(realpath $0)
TopDir=$(dirname $AbsPath0)
cd $TopDir/../diamondCase

# Clean everything up:
dcDir=/var/tmp/DiamondCases
rm -fr $dcDir
mkdir  $dcDir

# Generate all cases under "$dcDir":
# Angle of Attack, in deg:
for alphaDeg in {0..30..1}
do
  # Air speed at infinity:
  for Uoo in {10..990..10}
  do
    # Create the CaseDir:
    caseDir=$(printf $dcDir/%02d/%03d $alphaDeg $Uoo)
    mkdir -p $caseDir

    # Copy the prototype files over:
    cp -a 0 constant system VarParams $caseDir

    # Create the "AU" file:
    echo -e "__alphaDeg__ $alphaDeg;\n__Uoo__ $Uoo;" > $caseDir/AU
  done
done
# Generation Done!
