#! /bin/bash -e
#=============================================================================#
#                   "SU2/Scripts/MkAllEWCases-Euler.sh":                      #
#=============================================================================#
AbsPath0=$(realpath $0)
TopDir=$(dirname $AbsPath0)
cd $TopDir/../EllipticWing
ewDir=~/Tmp/EW-Euler

#-----------------------------------------------------------------------------#
# Top Matter:                                                                 #
#-----------------------------------------------------------------------------#
# Generate the Mesh File with the default params:
$TopDir/../../__BUILD__/GCC-Debug/bin/MkZhukovskyMesh \
	-N 512 \
	-o $ewDir/MeshEW.su2

#-----------------------------------------------------------------------------#
# Generate All Cases under "$ewDir":                                          #
#-----------------------------------------------------------------------------#
# Angle of Attack, in deg:
for alphaDeg in {0..21..3}
do
# NB: The Convergence Criterion Field depends on "alphaDeg":
  if [ $alphaDeg -eq 0 ]
  then
    ConvFld=DRAG
  else
    ConvFld=LIFT
  fi

  # Mach Number at infinity (as a multiple of 0.01):
  for M100 in {10..300..5}
  do
    M=$(echo $M100 | awk '{printf("%.02f", $1/100)}')

    # Create the CaseDir:
    caseDir=$(printf $ewDir/%02d/%s $alphaDeg $M)
    mkdir -p $caseDir

    # Copy the Config-Euler.cfg over into the corresp "Cfg":
    cp Config-Euler-Proto.cfg $caseDir/Cfg

		# Guess a reasonable initial CFL_NUMBER value:
		if   [ $M100 -le 75 -o $M100 -ge 95 ]
		then
			InitCFL=1
			Iter=250000
		else
			InitCFL=0.1
			Iter=1000000
		fi

    # Install the variable params (some of them may be adjusted manually
    # later):
    cat >> $caseDir/Cfg <<EOF
MACH_NUMBER= $M
AOA= $alphaDeg
CONV_FIELD= $ConvFld
RESTART_SOL= NO
ITER= $Iter
CFL_NUMBER= $InitCFL
CFL_ADAPT= NO
%CFL_ADAPT_PARAM= (0.8, 1.25, 0.01, 1.0, 1E-5)
EOF
    # There is no need to copy the Mesh file; a symlink will do:
    [[ -s $caseDir/MeshEW.su2 ]] || ln -s $ewDir/MeshEW.su2 $caseDir/MeshEW.su2
  done
done

#-----------------------------------------------------------------------------#
# Generate the Makefile:                                                      #
#-----------------------------------------------------------------------------#
cd $ewDir
shopt -s nullglob
AllCases=(??/?.??)
AllTargs=${AllCases[@]/%/\/DL}

echo "all: ${AllTargs[@]}" > Makefile
for c  in  ${AllCases[@]}
do
  alpha=$(dirname $c)
  M=$(basename $c)
  cat >> Makefile <<EOF

$c/DL: $c/Log
	cd $c && awk -v alpha=$alpha -v M=$M -f $TopDir/EWCoeffs.awk Log > DL
$c/Log: $c/Cfg
	cd $c && SU2_CFD Cfg >& Log
EOF
done
# Generation Done!
