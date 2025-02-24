#! /bin/bash -e
#=============================================================================#
#                      "SU2/Scripts/MkAllEWCases.sh":                         #
#=============================================================================#
AbsPath0=$(realpath $0)
TopDir=$(dirname $AbsPath0)
cd $TopDir/../EllipticWing

# Clean everything up:
ewDir=/var/tmp/EllipticWing
rm -fr $ewDir
mkdir  $ewDir

#-----------------------------------------------------------------------------#
# Generate All Cases under "$ewDir":                                          #
#-----------------------------------------------------------------------------#
# Angle of Attack, in deg:
for alphaDeg in {0..15..1}
do
  # Mach Number at infinity (as a multiple 0.01):
  for M100 in {10..300..5}
  do
    M=$(echo $M100 | awk '{printf("%.02f", $1/100)}')

    # Create the CaseDir:
    caseDir=$(printf $ewDir/%02d/%s $alphaDeg $M)
    mkdir -p $caseDir

    # Copy the Config.cfg over and install the variable params in it:
    cp Config.cfg $caseDir
    echo -e "MACH_NUMBER= $M\nAOA= $alphaDeg\n" >> $caseDir/Config.cfg

    # There is no need to copy the Mesh file; symlink will do:
    ln -s $PWD/Mesh.su2 $caseDir/Mesh.su2
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
    "\n$c/Coeffs: $c/Config.cfg\n\tcd $c &&" \
    "SU2_CFD Config.cfg >& cfd.log &&" \
    "awk -v alpha=$alpha -v M=$M -f $TopDir/EWCoeffs.awk cfd.log > Coeffs" \
    >> Makefile
done
# Generation Done!
