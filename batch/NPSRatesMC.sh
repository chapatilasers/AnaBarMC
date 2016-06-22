#!/bin/tcsh
#PBS -N NPSRatesMC
#PBS -m n
#PBS -M david.j.hamilton@glasgow.ac.uk
#PBS -l walltime=40:00:00
#PBS -V

setenv nevents 5000000

setenv tempdir "/scratch/david"
mkdir -p ${tempdir}

set MCMACRO = "/home/david/JLab/C/NPSRateMC/batchfarm/farm_$RUN_NUMBER.mac"

echo "/NPSRatesMC/physics/addPhysics QGSP_BIC_EMY"                       >   $MCMACRO
echo "/NPSRatesMC/physics/setCuts 0.1 mm"                                >>  $MCMACRO
echo "/run/initialize"                                                   >>  $MCMACRO
echo "/NPSRatesMC/generator/Mode 0"                                      >>  $MCMACRO
echo "/NPSRatesMC/generator/Seed $RUN_NUMBER"                            >>  $MCMACRO
echo "/NPSRatesMC/analysis/setOutputFile $tempdir/farm_$RUN_NUMBER.root" >>  $MCMACRO
echo "/gps/particle e-"                                                  >>  $MCMACRO
echo "/gps/pos/centre 0. 0. -200. mm"                                    >>  $MCMACRO
echo "/gps/pos/type Beam"                                                >>  $MCMACRO
echo "/gps/pos/shape Square"                                             >>  $MCMACRO
echo "/gps/pos/sigma_x 3 mm"                                             >>  $MCMACRO
echo "/gps/pos/sigma_y 3 mm"                                             >>  $MCMACRO
echo "/gps/ang/type beam2d"                                              >>  $MCMACRO
echo "/gps/ang/sigma_x 0.0 deg"                                          >>  $MCMACRO
echo "/gps/ang/sigma_y 0.0 deg"                                          >>  $MCMACRO
echo "/gps/ang/rot1 -1 0 0"                                              >>  $MCMACRO
echo "/gps/energy 6.6 GeV"                                               >>  $MCMACRO
echo "/run/beamOn $nevents"                                              >>  $MCMACRO

cd /home/david/JLab/C/NPSRateMC/
source /home/david/Misc/geant4/G4setup.csh
nohup /home/david/geant4/bin/Linux-g++/NPSRatesMC $MCMACRO #>& /dev/null
echo "****************** NPSRates MC Finished"

cp    ${tempdir}/"farm_$RUN_NUMBER.root"   ${OUTPUT_DIR}/
rm -f ${tempdir}/"farm_$RUN_NUMBER.root"
rm -f $MCMACRO
