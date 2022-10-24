#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ROOT as root
import numpy as np
import random
from timer import Timer

t = Timer()
t.start()


# In[2]:


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


# In[ ]:


f = root.TFile("data/AnaBarMC_77777.root")
f.ls()

myTree = f.Get("T")
myTree.Show(0)

myTree.Print()

entries = myTree.GetEntriesFast()

hPrimE = root.TH1F("PrimE","Primary Energy", 100, 0, 25000)
hPrimTh = root.TH1F("PrimTh","Primary Theta", 100, 0, np.pi)
hPrimPh = root.TH1F("PrimPh","Primary Phi", 100, 0, 2.0*np.pi)
hPrimPdg = root.TH1F("PrimPdg","Primary PDG ID", 30, 0, 30)


# In[ ]:

print (entries)
for i in range(entries):
    
    if (i%100 == 0):
        print("Event Number: %d" % i)
    
    anabar_hit_paddle = np.array([False for i in range(NMaxPMT)],dtype=bool)
    edeptot = np.array([0.0 for i in range(NMaxPMT)],dtype=float)
    edep0tot = 0.0
    
    myTree.GetEntry(i)
    
    trigger = False
    finger_hit = False
    anabar_hit = False
    anabar_top_hit = False
    anabar_bottom_hit = False
    j_finger = 0
    j_anabar = 0
    
    if (myTree.Prim_pdg == 11):
        fMass = 0.511
    else:
        if (myTree.Prim_pdg == 13):
            fMass = 105.7
        else:
            if (myTree.Prim_pdg == 2212):
                fMass = 938.28
            else:
                fMass = 939.5
    
    trigger2 = False
    trigger3 = False
    
    fMomentum = np.sqrt(myTree.Prim_E**2 - fMass**2)
    fPx = fMomentum*np.sin(myTree.Prim_Th)*np.cos(myTree.Prim_Ph)
    fPy = fMomentum*np.sin(myTree.Prim_Th)*np.sin(myTree.Prim_Ph)
    fPz = fMomentum*np.cos(myTree.Prim_Th)
    fNewTheta = np.arccos(fPy/fMomentum)
    if (fPx<0):
        fNewPhi = np.arctan(fPz/fPx)+np.pi
    else:
        if (fPx>0 and fPz<0):
            fNewPhi = np.arctan(fPz/fPx)+2.0*np.pi
        else:
            fNewPhi=np.arctan(fPz/fPx)
    
    for j in range(myTree.Detector_Nhits):
        if(myTree.Detector_id[j] == Detector_Offset and not finger_hit):
            finger_hit = True
            j_finger = j
            
        for ibar in range(1,15):
            if(myTree.Detector_id[j+Detector_Offset] == (ibar+Detector_Offset)):
                anabar_hit = True
                anabar_hit_paddle[ibar-1]=True
                j_anabar = j
            
    if (finger_hit and anabar_hit):
        trigger = True
        
    if trigger:
        hPrimE.Fill(myTree.Prim_E)
        hPrimTh.Fill(fNewTheta)
        hPrimPh.Fill(fNewPhi)
        hPrimPdg.Fill(myTree.Prim_pdg)


# In[ ]:


c2 = root.TCanvas("c1","c1",800,800)
c2.Divide(2,2,0.01,0.01,0)

c2.cd(1)
hPrimE.Draw()
c2.cd(2)
hPrimTh.Draw()
c2.cd(3)
hPrimPh.Draw()
c2.cd(4)
hPrimPdg.Draw()


# In[ ]:


c2.Draw()
c2.Print("plots/c2R.pdf")

t.stop()
