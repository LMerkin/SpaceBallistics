#! /bin/bash -e
#=============================================================================#
#                      "Scripts/RunDiamondCasesAll":                          #
#=============================================================================#
# This script must be started from the directory containing the "diamondCase":
if [ ! -d diamondCase ]
then
	echo "Run this script from a directory containing the 'diamondCase' sub-dir"
	exit 1
fi
cd diamondCase

# Speed of Sound (XXX: under the standard atmospheric conditions):
a=340.294

# Transonic Mach Limits:
Ul=$(echo | awk "{print int($a * 0.8)}")
Uu=$(echo | awk "{print int($a * 1.2)}")

# Clean everything up:
dcDir=/tmp/DiamondCases
rm -fr $dcDir
mkdir  $dcDir

# Generate all cases under /tmp:
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

		# Override the 0/U file depending on the curr Mach number:
		if   [[ $Uoo -lt $Ul ]]
		then
			cp templates/p-subsonic   $caseDir/0/p
		elif [[ $Uoo -lt $Uu ]]
		then
			cp templates/p-transonic  $caseDir/0/p
		else
			cp templates/p-supersonic $caseDir/0/p
		fi

		# Create the "AU" file:
		echo -e "__alphaDeg__ $alphaDeg;\n__Uoo__ $Uoo;" > $caseDir/AU
	done
done
# Generation Done!
