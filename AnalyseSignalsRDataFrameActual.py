#!/usr/bin/env python

import ROOT as root
import numpy as np
import random
from timer import Timer

t = Timer()
t.start()

Analyse_Secondaries = 1
Theta_min_cut = 2.524
ThetaVerticalCut = 3.02

MaxHits = 50000
MaxPMTNo = 20
MaxPMTHits = 5000
Finger_Edep_Max = 10.0
AnaBar_Edep_Max = 10.0
pedastel_sigma = 2.9
Detector_Offset = 0

Finger_NPhotons_Max = 250
AnaBar_NPhotons_Max = 200

NUMPADDLE = 14

NMaxPMT = 14

Detector_pdg = np.array([0 for i in range(MaxHits)],dtype=int)
Detector_id = np.array([0 for i in range(MaxHits)],dtype=int)

Detector_x = np.array([0 for i in range(MaxHits)],dtype=float)
Detector_y = np.array([0 for i in range(MaxHits)],dtype=float)
Detector_z = np.array([0 for i in range(MaxHits)],dtype=float)
Detector_t = np.array([0 for i in range(MaxHits)],dtype=float)
Detector_Ed = np.array([0 for i in range(MaxHits)],dtype=float)

PMT_Nphotons = np.array([0 for i in range(MaxPMTNo)],dtype=int)
PMT_Nphotons_Noise = np.array([0 for i in range(MaxPMTNo)],dtype=int)

PMT_KineticEnergy = np.array([[0 for i in range(MaxPMTNo)] for j in range(MaxPMTHits)],dtype=float)

fileName = "data/AnaBarMC_777777.root"
treeName = "T"

f = root.TFile(fileName)
myTree = f.Get(treeName)
myTree.GetEntry(0)
particleType = myTree.Prim_pdg

if (particleType == 11):
    fMass = 0.511
else:
    if (particleType == 13):
        fMass = 105.7
    else:
        if (particleType == 2212):
            fMass = 938.28
        else:
            fMass = 939.5

print ("Incident particle mass = ",fMass)

d = root.RDataFrame(treeName,fileName)

cut1 = '






