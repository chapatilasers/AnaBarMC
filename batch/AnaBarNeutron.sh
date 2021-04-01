#!/bin/bash
#PBS -N CDetOptical
#PBS -m n
#PBS -M lydia.lorenti.15@cnu.edu
#PBS -l walltime=40:00:00
#PBS -V

export CENTOS_VERSION=`rpm --eval '%{centos_ver}'`

if [ $CENTOS_VERSION == "7" ]
then
	source /home/brash/geant4_C7/G4setup_batch.sh
	export G4BINARY=/home/brash/geant4_C7/bin/Linux-g++/AnaBarMC
else
	source /home/brash/geant4_C8/G4setup_batch.sh
	export G4BINARY=/home/brash/geant4_C8/bin/Linux-g++/AnaBarMC
fi

export nevents=100
export tempdir=/home/brash/CDetOptical/batch

export MACRO_PATH=/home/brash/CDetOptical/macros/
export MCMACRO=$tempdir/AnaBarMC_$RUN_NUMBER.mac

echo "/control/macroPath $MACRO_PATH"	 	                         >   $MCMACRO
echo "/AnaBarMC/physics/addPhysics standard_opt3"                        >>   $MCMACRO
echo "/AnaBarMC/physics/optical 1"	                                 >>  $MCMACRO
echo "/AnaBarMC/detector/AnaBarXpos 0.00"	                         >>  $MCMACRO
echo "/run/initialize"                                                   >>  $MCMACRO
echo "/AnaBarMC/generator/Mode 1"                              >>  $MCMACRO
echo "/AnaBarMC/generator/InputFile $tempdir/data/AnaBarMC_Gen_$RUN_NUMBER.root" >>  $MCMACRO
echo "/AnaBarMC/analysis/setOutputFile $tempdir/rootfiles/AnaBarMC_$RUN_NUMBER.root" >>  $MCMACRO

cd $tempdir
export ROOTSYS=/cern/root/pro
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
export PATH=$ROOTSYS/bin:$PATH
export DISPLAY=jlabanalysis.pcs.cnu.edu:0.0
#nohup root -l -q GenCosmics.C++\($nevents,$RUN_NUMBER\) #>& /dev/null
nohup $G4BINARY $MCMACRO #>& /dev/null
echo "****************** AnaBarMC Finished"

cp    ${tempdir}/rootfiles/"AnaBarMC_$RUN_NUMBER.root"   ${OUTPUT_DIR}/
rm -f ${tempdir}/rootfiles/"AnaBarMC_$RUN_NUMBER.root"
rm -f $MCMACRO
