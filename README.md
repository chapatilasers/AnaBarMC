# AnaBarMC
To compile the code type:
make or gmake

To run the event generator type:
root -q GenCosmics.C++

This creates a root file in the data/ directory

To run the simulation interactively type:

AnaBarMC [return]
then at the prompt type:
/control/execute macros/vis.mac

To run the simaultion in batch mode type:

AnaBarMC macros/batch.mac




Tanner's Comments:
Change Primary Particle in ~/CDetOptical/batch/GenCosmics.C at line 100
