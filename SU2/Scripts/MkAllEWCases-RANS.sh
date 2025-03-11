#! /bin/bash -e
#=============================================================================#
#                    "SU2/Scripts/MkAllEWCases-RANS.sh":                      #
#=============================================================================#
AbsPath0=$(realpath $0)
TopDir=$(dirname $AbsPath0)
cd $TopDir/../EllipticWing
ewDir=~/Tmp/EW-RANS

#-----------------------------------------------------------------------------#
# Generate All Cases under "$ewDir":                                          #
#-----------------------------------------------------------------------------#
# Angle of Attack, in deg:
for alphaDeg in {0..21..3}
do
  # Mach Number at infinity (as a multiple of 0.01):
  for M100 in {10..300..5}
  do
		# The actual Mach number, the FreeStream Velocity and the Raynolds Number:
		# XXX: BEWARE: We assume the Speed of Sound to be @ 15C, ie 340.294 m/sec:
		Params=($(echo $M100 | awk '{printf("%.02f %.6f", $1/100, $1*231492.517)}'))
		M=${Params[0]}
		Re=${Params[1]}

    # Create the CaseDir:
    caseDir=$(printf $ewDir/%02d/%s $alphaDeg $M)
		caseDir0=$ewDir/00/$M
    mkdir -p $caseDir

		# Generate the Mesh File (or link an existing one):
		if [ $alphaDeg -eq 0 ]
		then
			$TopDir/../../__BUILD__/GCC-Debug/bin/MkZhukovskyMesh \
					-N 512 -y 70 -M $M -o $caseDir/MeshEW.su2
		else
			ln -s $caseDir0/MeshEW.su2 $caseDir/MeshEW.su2
		fi

    # Copy the Config-RANS.cfg over into the corresp "Cfg":
    cp Config-RANS-Proto.cfg $caseDir/Cfg

    # Install the variable params (some of them may be adjusted manually later).
    # XXX: For the moment, we assume Iter=25000 and CFLN=1.
    #      For Convergence, use Drag:
    cat >> $caseDir/Cfg <<EOF
ITER= 35000
CFL_NUMBER= 1
MACH_NUMBER= $M
AOA= $alphaDeg
REYNOLDS_NUMBER= $Re
RESTART_SOL= NO
EOF
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
