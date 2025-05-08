#! /bin/bash -e
#=============================================================================#
#                    "SU2/Scripts/MkAllEWCases-RANS.sh":                      #
#=============================================================================#
AbsPath0=$(realpath $0)
TopDir=$(dirname $AbsPath0)
cd $TopDir/../EllipticWing/Configs
ewDir=~/Tmp/EW-RANS

#-----------------------------------------------------------------------------#
# Generate All Cases under "$ewDir":                                          #
#-----------------------------------------------------------------------------#
# Angle of Attack, in deg:
for alphaDeg in {0..12..3}
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

    # We currently DO NOT use Wall Functions -- they may apparently do more harm
    # than good from the convergence point of view.  Without them, use y+ = 1:
    # The Meshes are the same for all AoA; so we can generate them for AoA=0,
    # and use symlinks otherwise;
    # Convergence is measured in CD for alphaDeg=0, and in CL otherwise;
    # "yPlus" in the Mesg and the number of iterations in the Cfg file depend on
    # the Mach number:
    yPlus=1
		CFLInit=5
		CFLMax=500
    if [ $M100 -ge 150 ]
    then
     # In "sufficiently high super-sonic modes", a more coarse "yPlus" is OK:
     yPlus=10
    fi
		if [ $M100 -ge 240 ]
		then
			CFLInit=1
		  CFLMax=100
		fi

    # Generate (or link) the Mesh File:
    if [ $alphaDeg -eq 0 ]
    then
      $TopDir/../../__BUILD__/GCC-Debug/bin/MkZhukovskyMesh \
          -N 256 -y $yPlus -M $M -s 1000 -R 1.1 -n 1 -o $caseDir/MeshEW.su2
      ConvField=DRAG
    else
      ln -sf $caseDir0/MeshEW.su2 $caseDir/MeshEW.su2
      ConvField=LIFT
    fi

    # Copy the config into the corresp "Cfg":
    cp Config-RANS-NoWF-Proto.cfg $caseDir/Cfg

    # Install the variable params (some of them may be adjusted manually later):
    cat >> $caseDir/Cfg <<EOF
CFL_ADAPT=              YES
CFL_NUMBER=             $CFLInit
CFL_ADAPT_PARAM=        ( 0.1, 1.2, 0.1, $CFLMax )
ITER=                   100000
CONV_FIELD=             $ConvField
MACH_NUMBER=            $M
AOA=                    $alphaDeg
REYNOLDS_NUMBER=        $Re
RESTART_SOL=            NO
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
$c/Log: $c/Cfg $c/MeshEW.su2
	cd $c && SU2_CFD Cfg >& Log
EOF
done
# Generation Done!
