#! /bin/bash -e
#=============================================================================#
#                        "Scipts/RunAllDiamondCases.sh":                      #
#=============================================================================#
ewDir=~/Tmp/EW
cd $ewDir

# Run "make" for all cases:
make >& make.log

# Collect the Coeffs for each AngleOfAttack, for all M vals:
shopt -s nullglob
AllCoeffs=(*/*/Coeffs)
[[ ${#AllCoeffs[@]} -eq 0 ]] || cat ${AllCoeffs[@]} > Coeffs

