#! /bin/bash -e
#=============================================================================#
#                      "SU2/Scripts/MkAllEWCases.sh":                         #
#=============================================================================#
AbsPath0=$(realpath $0)
TopDir=$(dirname $AbsPath0)
cd $TopDir/../EllipticWing

NProc=$(nproc)

#-----------------------------------------------------------------------------#
# Top Matter:																																	#
#-----------------------------------------------------------------------------#
# Clean everything up:
ewDir=~/Tmp/EW
rm -fr $ewDir
mkdir  $ewDir

# Generate the Mesh File:
$TopDir/../../__BUILD__/GCC-Debug/bin/MkZhukovskyMesh \
	-A 10 -N 256 -o $ewDir/MeshEW.su2

#-----------------------------------------------------------------------------#
# Generate All Cases under "$ewDir":                                          #
#-----------------------------------------------------------------------------#
# Angle of Attack, in deg:
for alphaDeg in {0..30..1}
do
	# NB: The Convergence Criterion Field depends on "alphaDeg":
	if [ $alphaDeg -eq 0 ]
	then
		ConvFld=DRAG
	else
		ConvFld=LIFT
	fi

  # Mach Number at infinity (as a multiple 0.01):
  for M100 in {10..300..5}
  do
    M=$(echo $M100 | awk '{printf("%.02f", $1/100)}')

    # Create the CaseDir:
    caseDir=$(printf $ewDir/%02d/%s $alphaDeg $M)
    mkdir -p $caseDir

    # Copy the Config-Euler.cfg over:
    cp Config-Euler.cfg $caseDir/C.cfg

    # Install the variable params.
    # NB: The MaxCFL depends on the Mach Number:
    if [ $M100 -le 30 -o $M100 -ge 170 ]
		then
			MaxCFL=1
		else
			MaxCFL=0.1
		fi
    cat >> $caseDir/C.cfg <<EOF
MACH_NUMBER= $M
AOA= $alphaDeg
CFL_ADAPT_PARAM= (0.618, 1.618, 0.01, $MaxCFL, 1E-5)
CONV_FIELD= $ConvFld
EOF

    # There is no need to copy the Mesh file; symlink will do:
    ln -s $ewDir/MeshEW.su2 $caseDir/MeshEW.su2
  done
done

#-----------------------------------------------------------------------------#
# Generate the Makefile:                                                      #
#-----------------------------------------------------------------------------#
cd $ewDir
shopt -s nullglob
AllCases=(*/*)
AllTargs=${AllCases[@]/%/\/Coeffs}

echo "all: ${AllTargs[@]}" > Makefile
for c  in  ${AllCases[@]}
do
  alpha=$(dirname $c)
  M=$(basename $c)
  echo -e \
    "\n$c/Coeffs: $c/C.cfg\n\tcd $c &&" \
    "SU2_CFD -t $NProc C.cfg >& cfd.log &&" \
    "awk -v alpha=$alpha -v M=$M -f $TopDir/EWCoeffs.awk cfd.log > Coeffs" \
    >> Makefile
done
# Generation Done!
