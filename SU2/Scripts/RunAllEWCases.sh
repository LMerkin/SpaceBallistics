#! /bin/bash -e
#=============================================================================#
#                        "Scipts/RunAllDiamondCases.sh":                      #
#=============================================================================#
# How many processes to run in parallel:
NPar=$(( $(nproc) / 2 ))

AbsPath0=$(realpath $0)
TopDir=$(dirname $AbsPath0)

shopt -s nullglob
ewDir=~/Tmp/EllipticWing
cd $ewDir

#=============================================================================#
AllCases=(*/*)
#=============================================================================#
AllTargs=${AllCases[@]/%/\/Coeffs}

# Generate the Makefile:
echo "all: ${AllTargs[@]}" > Makefile
for c  in  ${AllCases[@]}
do
  alpha=$(dirname $c)
  M=$(basename $c)
  echo -e \
    "\n$c/Coeffs: $c/Config.cfg\n\tcd $c &&" \
    "SU2_CFD Config.cfg >& cfd.log &&" \
    "awk -v alpha=$alpha -v M=$M -f $TopDir/EWCoeffs.awk cfd.log > Coeffs" \
    >> Makefile
done

# Run "make" for all cases, creating at most $NPar parallel processes:
make -j $NPar

# Collect the Coeffs for each AngleOfAttack, for all M vals:
AllCoeffs=(*/*/Coeffs)
[[ ${#AllCoeffs[@]} -eq 0 ]] || cat ${AllCoeffs[@]} > Coeffs

