#! /bin/bash -e
#=============================================================================#
#                        "Scipts/RunAllDiamondCases.sh":                      #
#=============================================================================#
# How many processes to run in parallel:
NPar=$(( $(nproc) / 2 ))

ewDir=/var/tmp/EllipticWing
cd $ewDir

# Run "make" for all cases, creating at most $NPar parallel processes:
make -j $NPar

# Collect the Coeffs for each AngleOfAttack, for all M vals:
shopt -s nullglob
AllCoeffs=(*/*/Coeffs)
[[ ${#AllCoeffs[@]} -eq 0 ]] || cat ${AllCoeffs[@]} > Coeffs

