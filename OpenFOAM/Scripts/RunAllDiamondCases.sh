#! /bin/bash -e
#=============================================================================#
#                        "Scipts/RunAllDiamondCases.sh":                      #
#=============================================================================#
# How manu processes to run in parallel:
NPar=8

AbsPath0=$(realpath $0)
TopDir=$(dirname $AbsPath0)
a=$(awk '/__a0__/{print ($2 + 0)}' $TopDir/../diamondCase/VarParams)

shopt -s nullglob
dcDir=/var/tmp/DiamondCases
cd $dcDir

#=============================================================================#
#AllCases=($*/*)
AllCases=(05/{{0..5}?0,6{0..6}0})
#=============================================================================#
AllTargs=${AllCases[@]/%/\/Coeffs}

# Generate the Makefile cotaining dependenies on the srcs likely to change:
echo "all: ${AllTargs[@]}" > Makefile
for c  in  ${AllCases[@]}
do
  V=$(basename $c)
  echo -e \
    "\n$c/Coeffs: $c/system/fvSchemes $c/system/fvSolution" \
    "$c/constant/thermophysicalProperties"   \
    "$c/constant/turbulenceProperties\n\tcd $c &&" \
    "rhoSimpleFoam >& RSF.log &&" \
    "rhoSimpleFoam -postProcess -func forceCoeffs -latestTime >& PP.log &&" \
    "awk -v V=$V -v a=$a -f $TopDir/DiamondCoeffs.awk PP.log > Coeffs" \
    >> Makefile
done
exit 0

# Run "make" for all cases, creating at most $NPar parallel processes:
make -j $NPar

# Collect the Coeffs for each AngleOfAttack, for all M vals:
AllCoeffs=(*/*/Coeffs)
[[ ${#AllCoeffs[@]} -eq 0 ]] || cat ${AllCoeffs[@]} > Coeffs

