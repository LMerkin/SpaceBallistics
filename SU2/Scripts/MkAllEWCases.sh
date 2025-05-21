#! /bin/bash -e
#=============================================================================#
#                         "SU2/Scripts/MkAllEWCases.sh":                      #
#=============================================================================#
AbsPath0=$(realpath $0)
TopDir=$(dirname $AbsPath0)
cd $TopDir/../EllipticWing/Configs

# The Model is "RANS" by default, unless overridden by "$1". The following
# Models are supported:
Model="RANS"
[[ "$1" != "" ]] && Model="$1"

# XXX: For now, only the following Models are supported:
if [ "$Model" != "RANS" -a "$Model" != "RANS-WF" -a "$Model" != "NS" -a \
     "$Model" == "Euler" ]
then
  echo "Invalid Model: $Model"
  exit 1
fi

ewDir=~/Tmp/EW-"$Model"
echo $ewDir
exit 0

#-----------------------------------------------------------------------------#
# Generate All Cases under "$ewDir":                                          #
#-----------------------------------------------------------------------------#
# Angle of Attack, in deg:
for alpha10 in {0..90..15}
do
  alphaDeg=$(echo $alpha10 | awk '{printf("%.01f", $1/10)}')

  # Mach Number at infinity (as a multiple of 0.01):
  for M100 in {10..300..5}
  do
    # The actual Mach number, the FreeStream Velocity and the Raynolds Number:
    # XXX: BEWARE: We assume the Speed of Sound to be @ 15C, ie 340.294 m/sec:
    Params=($(echo $M100 | awk '{printf("%.02f %.6f", $1/100, $1*231492.517)}'))
    M=${Params[0]}
    Re=${Params[1]}

    # Create the CaseDir:
    caseDir=$(printf $ewDir/%s/%s $alphaDeg $M)
    caseDir0=$ewDir/0.0/$M
    mkdir -p $caseDir

    # We currently DO NOT use Wall Functions -- they may apparently do more harm
    # than good from the convergence point of view.  Without them, use y+ = 1:
    # The Meshes are the same for all AoA; so we can generate them for AoA=0,
    # and use symlinks otherwise;
    # Convergence is measured in CD for alphaDeg=0, and in CL otherwise;
    # "yPlus" in the Mesg and the number of iterations in the Cfg file depend on
    # the Mach number:
    CFLInit=5
    CFLMax=500
    if [ $M100 -ge 240 ]
    then
      CFLInit=1
      CFLMax=100
    fi

    case "$Model" in
      "RANS")    yPlus=1;   [[ $M100 -ge 150 ]] && yPlus=10;;
      "RANS-WF") yPlus=100;;
      "NS")      yPlus=10;;    # Not tested much
      "Euler")   yPlus=0;;     # Not used!
    esac

    # Generate (or link) the Mesh File:
    if [ $alpha10 -eq 0 ]
    then
      if [ "$Model" != "Euler" ]
      then
        $TopDir/../../__BUILD__/GCC-Debug/bin/MkZhukovskyMesh \
          -N 256 -y $yPlus -M $M -s 1000 -R 1.1 -n 1 -o $caseDir/MeshEW.su2
      else
        $TopDir/../../__BUILD__/GCC-Debug/bin/MkZhukovskyMesh \
          -N 512 -o $caseDir/MeshEW.su2
      fi
      ConvField=DRAG
    else
      ln -sf $caseDir0/MeshEW.su2 $caseDir/MeshEW.su2
      ConvField=LIFT
    fi

    # Copy the config into the corresp "Cfg":
    cp Config-"$Model"-Proto.cfg $caseDir/Cfg

    # Install the variable params (some of them may be adjusted manually later):
    cat >> $caseDir/Cfg <<EOF
CFL_ADAPT=              YES
CFL_NUMBER=             $CFLInit
CFL_ADAPT_PARAM=        ( 0.1, 1.2, 0.1, $CFLMax )
ITER=                   250000
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
AllCases=(?.?/?.??)
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
