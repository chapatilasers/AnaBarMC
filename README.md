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

---Tanner's Comments---

To Run Simulation:

-Go into /CDetOptical/batch directory and type: "root generate_80.C" to create the primary particle rootfiles that are then store in /CDetOptical/batch/data
                                                                                  
-Type: "./sendjobs.sh runlist80" to send the primary particle rootfiles that are   listed in the according runlist file that is sent as an argumant. The output data will be stored in /CDetOptical/data
                                                                                  
Changing the Primary Particle:

-Switch the Generator Function in ~/CDetOptical/batch/GenCosmics.C at line 100

-Change PrimaryParticleID in ~/CDetOptical/AnalyseSignals.C at line 28
