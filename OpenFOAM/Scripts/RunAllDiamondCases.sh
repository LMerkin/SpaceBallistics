#! /bin/bash -e
#=============================================================================#
#                        "Scipts/RunAllDiamondCases.sh":                      #
#=============================================================================#
# How manu processes to run in parallel:
NPar=8

shopt -s nullglob
#AllCases=(/tmp/DiamondCases/*/*)
AllCases=(/tmp/DiamondCases/10/2*)

NC=${#AllCases[@]}
from=0
while [[ $from -lt $NC ]]
do
	# Take at most $NPar list elements:
	to=$((from + NPar - 1))
	[[ $to -lt $NC ]] || to=$((NC - 1))

	# Run the batch {$from .. $to} of cases:
	for j in $(seq $from $to)
	do
		# Run the "j"th case:
		caseDir=${AllCases[$j]}
		V=$(basename $caseDir)
		( \
		 cd $caseDir && \
 	   rhoSimpleFoam >& RSF.log && \
     rhoSimpleFoam -postProcess -func forceCoeffs -latestTime >& PP.log && \
		 awk -v V=$V '/Cd:/{Cd=$2}/Cl:/{Cl=$2}END{print V, V/340.294, Cd, Cl}' \
			 PP.log > Coeffs \
		) &
	done
	# Wait for all parallel processes to complete:
	wait
	# To the next batch:
	from=$((to + 1))
done

# Collect the Coeffs for each AngleOfAttack, for all M vals:
for a in /tmp/DiamondCases/*
do
	cd $a
	cat */Coeffs > Coeffs
done
