#!/usr/bin/env python
# coding: utf-8

# In[2]:


import ROOT as root
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random


# In[3]:


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


# In[4]:


# In[12]:


f = root.TFile("data/AnaBarMC_777777.root")
f.ls()

myTree = f.Get("T")
myTree.Show(0)

myTree.Print()

entries = myTree.GetEntriesFast()

hFingerX = root.TH1F("FingerX","Finger X Position",100,-120,120)
hFingerY = root.TH1F("FingerY","Finger Y Position",100,30,80)
hFingerZ = root.TH1F("FingerZ","Finger Z Position",100,-140,60)
hFingerT = root.TH1F("FingerT","Finger Time",100,0.0,0.4)

hPrimE = root.TH1F("PrimE","Primary Energy", 100, 0, 25000)
hPrimTh = root.TH1F("PrimTh","Primary Theta", 100, 0, np.pi)
hPrimPh = root.TH1F("PrimPh","Primary Phi", 100, 0, 2.0*np.pi)
hPrimPdg = root.TH1F("PrimPdg","Primary PDG ID", 30, 0, 30)

hDetectorNhits = root.TH1F("DetectorNhits","Detector Number of Hits", 100, 0, 400)
hDetectorPdg = root.TH1F("DetectorPdg","Detector PDG ID", 50, -20, 30)
hDetectorID = root.TH1F("DetectorID","Detector ID Number", 30, 0, 30)
hPMTID = root.TH1F("PMTID","PMT ID Number", 15, 0, 15)

hFingerPMTNphot = root.TH1F("FingerPMTNphot","Finger PMT Number of Photons", Finger_NPhotons_Max+10, -10, Finger_NPhotons_Max)
hFingerEd = root.TH1F("FingerEd","Finger Energy Deposited", 100, 0.01, Finger_Edep_Max)

hAnaBarPMTNphot = []
hAnaBarPMTNphotCut = []

for i in range(NUMPADDLE):  
    name = ("AnaBarPMTNphotA%d" % i)
    title = ("AnaBar_PMT_Number_of_Photons_A%d" % i)
    print (name,title)
    hAnaBarPMTNphot.append(root.TH1F(name, title, int(AnaBar_NPhotons_Max*0.9+20), -20, int(AnaBar_NPhotons_Max*0.9)))
    
for i in range(NUMPADDLE):  
    name = ("AnaBarPMTNphotA%dCut" % i)
    title = ("AnaBar_PMT_Number_of_Photons_A%d_Cut" % i)
    print (name,title)
    hAnaBarPMTNphotCut.append(root.TH1F(name, title, int(AnaBar_NPhotons_Max*0.9+20), -20, int(AnaBar_NPhotons_Max*0.9)))

hAnaBarEd = root.TH1F("AnaBarEd","AnaBar Energy Deposited", 100, 0.01, AnaBar_Edep_Max)
hAnaBarEdCut = root.TH1F("AnaBarEdCut","AnaBar Energy Deposited", 100, 0.01, AnaBar_Edep_Max)

hAnaBarX = root.TH1F("AnaBarX","AnaBar X Position", 100, -120, 120)
hAnaBarY = root.TH1F("AnaBarY","AnaBar Y Position", 100, -30, 30)
hAnaBarZ = root.TH1F("AnaBarZ","AnaBar Z Position", 100, -30, 30)
hAnaBarT = root.TH1F("AnaBarT","AnaBar Time", 100, 0, .4)

hE1vsE2 = root.TH2F("E1vsE2", "AnaBar Edep vs. Finger Edep", 100, 0.01, Finger_Edep_Max, 100, 0.01, AnaBar_Edep_Max)

hFinger_Edep_vs_Nphot = root.TH2F("FingerEdepVsNphot", "Finger Edep vs. Number of Photons", Finger_NPhotons_Max, 0, Finger_NPhotons_Max, 100, 0.01, Finger_Edep_Max)
hNphot0_vs_Nphot1 = root.TH2F("AnaBarVsFingerNphot", "AnaBar vs. Finger Number of Photons", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max, Finger_NPhotons_Max, 0, Finger_NPhotons_Max)
hAnaBar_Edep_vs_Nphot = root.TH2F("AnaBarEdepVsNphot", "AnaBar Edep vs. Number of Photons", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max, 100, 0.01, AnaBar_Edep_Max)
hAnaBar_Edep_vs_NphotCut = root.TH2F("AnaBarEdepVsNphotCut", "AnaBar Edep vs. Number of Photons", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max, 100, 0.01, AnaBar_Edep_Max)

PMT_Nphotons_Total = 0

hAnaBarEdAll = []

for i in range(NUMPADDLE):
    name = ("AnaBarEd%d" % i)
    title = ("AnaBar Energy Deposited A%d" % i)
    hAnaBarEdAll.append(root.TH1F(name, title, 100, 0.01, AnaBar_Edep_Max))
    
hAnaBarEdAllCut = []

for i in range(NUMPADDLE):
    name = ("AnaBarEd%dCut" % i)
    title = ("AnaBar Energy Deposited A%d" % i)
    hAnaBarEdAllCut.append(root.TH1F(name, title, 100, 0.01, AnaBar_Edep_Max))

hAnaBarPMTNoiseCutNphot = []

for i in range(NUMPADDLE):
    name = ("AnaBarPMTNoiseCutNphotA%d" % i)
    title = ("AnaBar PMT Number of Photons A%d" % i)
    hAnaBarPMTNoiseCutNphot.append(root.TH1F(name, title, AnaBar_NPhotons_Max+20, -20, AnaBar_NPhotons_Max))
   # hAnaBarPMTNoiseCutNphot[i].SetLineColor("red")

hFingerPMTKE = root.TH1F("FingerPMTKE","Photon Wavelength Production Spectrum", 400, 300.0, 700.0)
hAnaBarPMTKEA1 = root.TH1F("AnaBarPMTKEA1","Photon Wavelength in WLS at PMT", 400, 300.0, 700.0)

hAnaBarMult = root.TH1F("AnaBarMult","Anabar PMT Multiplicity",12,0,12)

hPrimPx = root.TH1F("PrimPx","Primary Px", 100, -2000, 2000)
hPrimPz = root.TH1F("PrimPz","Primary Pz", 100, -2000, 2000)
hPrimPy = root.TH1F("PrimPy","Primary Py", 100, -2000, 0)

counter = 0

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
        
    if (finger_hit and anabar_hit and fNewTheta > Theta_min_cut):
        trigger2 = True
        
    if (finger_hit and anabar_hit and fNewTheta > ThetaVerticalCut):
        trigger3 = True
        
    if trigger:
        for icount in range(15):
            PMT_Nphotons_Noise[icount]=myTree.PMT_Nphotons[icount]+np.random.normal(loc=0.0,scale=pedastel_sigma)
                
        for i in range(NUMPADDLE):
            hAnaBarPMTNphot[i].Fill(PMT_Nphotons_Noise[i])
                
        hFingerPMTNphot.Fill(PMT_Nphotons_Noise[14])
    
    for j in range(myTree.Detector_Nhits):
        if trigger:
            if (myTree.Detector_id[j] == Detector_Offset and myTree.Detector_pdg[j] == myTree.Prim_pdg):
                counter = counter + 1
                hFingerX.Fill(myTree.Detector_x[j])
                hFingerY.Fill(myTree.Detector_y[j])
                hFingerZ.Fill(myTree.Detector_z[j])
                hFingerT.Fill(myTree.Detector_t[j])
                
            if (myTree.Detector_id[j] == Detector_Offset):
                hDetectorPdg.Fill(myTree.Detector_pdg[j])
                hDetectorID.Fill(myTree.Detector_id[j])
                
            if (myTree.Detector_id[j] == (1 + Detector_Offset)):
                hDetectorPdg.Fill(myTree.Detector_pdg[j])
                hDetectorID.Fill(myTree.Detector_id[j])
            
            if (myTree.Detector_id[j] == (1 + Detector_Offset) and myTree.Detector_pdg[j] == myTree.Prim_pdg):
                hAnaBarX.Fill(myTree.Detector_x[j])
                hAnaBarY.Fill(myTree.Detector_y[j])
                hAnaBarZ.Fill(myTree.Detector_z[j])
                hAnaBarT.Fill(myTree.Detector_t[j])

    #print ("Counter = ",counter)
    
    for j in range(myTree.Detector_Nhits):

        if (trigger):
            if (myTree.Detector_id[j] == Detector_Offset ):
                if (Analyse_Secondaries == 1 and fNewTheta > Theta_min_cut):
                    edep0tot += myTree.Detector_Ed[j];
                else:
                    if (myTree.Detector_pdg[j] == myTree.Prim_pdg and fNewTheta > Theta_min_cut):
                        edep0tot += myTree.Detector_Ed[j];
                        
            if (myTree.Detector_id[j] > Detector_Offset and myTree.Detector_id[j] <= NMaxPMT+Detector_Offset):
                if (Analyse_Secondaries == 1 and fNewTheta > Theta_min_cut):
                    edeptot[myTree.Detector_id[j]-1-Detector_Offset] += myTree.Detector_Ed[j]
                else:
                    if (myTree.Detector_pdg[j] == myTree.Prim_pdg and fNewTheta > Theta_min_cut):
                        edeptot[myTree.Detector_id[j]-1-Detector_Offset] += myTree.Detector_Ed[j]
                          
    if trigger:
        hFingerEd.Fill(edep0tot)
        hAnaBarEd.Fill(edeptot[6])
        
        hE1vsE2.Fill(edep0tot, edeptot[6])
        
        #fsum = 0
        #asum = 0
        #for i in range(NUMPADDLE+1):
        #    if (i<NUMPADDLE):
        #        asum += myTree.PMT_Nphotons[i]
        #    else:
        #        fsum += myTree.PMT_Nphotons[i]
        #        
        #print("Totals: ",fsum,asum)
        
        for iq in range(NUMPADDLE):
            #print(iq,myTree.PMT_Nphotons[iq])
            for jq in range(myTree.PMT_Nphotons[iq]):
                #if (myTree.PMT_KineticEnergy[iq*5000+jq]>0.1):
                #    print("Paddle: ",iq,iq*5000+jq,myTree.PMT_KineticEnergy[iq*5000+jq])
                if (myTree.PMT_KineticEnergy[iq*5000+jq] != 0):
                    hAnaBarPMTKEA1.Fill(1240.0/myTree.PMT_KineticEnergy[iq*5000+jq])
        
        for jq in range(myTree.PMT_Nphotons[14]):
            #if (myTree.PMT_KineticEnergy[14*5000+jq]>0.1):
            #    print("Finger: ",14*5000+jq,myTree.PMT_KineticEnergy[14*5000+jq])
            if (myTree.PMT_KineticEnergy[14*5000+jq] != 0):
                hFingerPMTKE.Fill(1240.0/myTree.PMT_KineticEnergy[14*5000+jq])
        
        #for i in range(20*5000):
        #    if (myTree.PMT_KineticEnergy[i] > 0.1):
        #        print (i,myTree.PMT_KineticEnergy[i])
        
        imult = 0
        for icount in range(NUMPADDLE):
            if(PMT_Nphotons_Noise[icount]>=8):
                imult += 1
        
        hAnaBarMult.Fill(imult)
        
    if trigger2:
        PMT_Nphotons_Total = 0
        for icount in range(14):
            PMT_Nphotons_Total += myTree.PMT_Nphotons[icount]
        
        hFinger_Edep_vs_Nphot.Fill(myTree.PMT_Nphotons[14],edep0tot)
        hNphot0_vs_Nphot1.Fill(PMT_Nphotons_Total,myTree.PMT_Nphotons[14])
        hAnaBar_Edep_vs_Nphot.Fill(myTree.PMT_Nphotons[6],edeptot[6])
            
        for i in range(NUMPADDLE):
            hAnaBarEdAll[i].Fill(edeptot[i])
        
    
    if trigger3:
        hAnaBarEdCut.Fill(edeptot[6])
        for i in range(NUMPADDLE):
            hAnaBarPMTNphotCut[i].Fill(PMT_Nphotons_Noise[i])
            hAnaBarEdAllCut[i].Fill(edeptot[i])
            if(anabar_hit_paddle[i] and edeptot[i]>=4.0):
                hAnaBarPMTNoiseCutNphot[i].Fill(PMT_Nphotons_Noise[i])
                
        hAnaBar_Edep_vs_NphotCut.Fill(myTree.PMT_Nphotons[6],edeptot[6])
        
    hPrimE.Fill(myTree.Prim_E)
    hPrimTh.Fill(fNewTheta)
    hPrimPh.Fill(fNewPhi)
    hPrimPdg.Fill(myTree.Prim_pdg)
    
    hPMTID.Fill(myTree.PMT_id)
    hDetectorNhits.Fill(myTree.Detector_Nhits)
    
    hPrimPx.Fill(fPx)
    hPrimPy.Fill(fPy)
    hPrimPz.Fill(fPz)
    

c1 = root.TCanvas("c1","c1",800,800)
c1.Divide(2,2,0.01,0.01,0)

c1.cd(1)
hFingerX.Draw()
c1.cd(2)
hFingerY.Draw()
c1.cd(3)
hFingerZ.Draw()
c1.cd(4)
hFingerT.Draw()


# In[13]:


c1.Draw()
c1.Print("plots/c1.pdf")


# In[14]:


c2 = root.TCanvas("c2","c2",800,800)
c2.Divide(2,2,0.01,0.01,0)

c2.cd(1)
hPrimE.Draw()
c2.cd(2)
hPrimTh.Draw()
c2.cd(3)
hPrimPh.Draw()
c2.cd(4)
hPrimPdg.Draw()


# In[15]:


c2.Draw()
c2.Print("plots/c2.pdf")


# In[16]:


c3 = root.TCanvas("c3","c3",800,800)
c3.Divide(2,2,0.01,0.01,0)

c3.cd(1)
hDetectorNhits.Draw()
c3.cd(2)
hDetectorPdg.Draw()
c3.cd(3)
hDetectorID.Draw()
c3.cd(4)
hPMTID.Draw()


# In[17]:


c3.Draw()
c3.Print("plots/c3.pdf")


# In[18]:


c4 = root.TCanvas("c4", "c4", 800,800)
c4.Divide(2,2, 0.01, 0.01, 0)
  
c4.cd(1)
hFingerEd.Draw()
c4.cd(2);
hFingerPMTNphot.Draw()
c4.cd(3)
hAnaBarEd.Draw()
c4.cd(4)
hAnaBarPMTNphot[0].Draw()


# In[19]:


c4.Draw()
c4.Print("plots/c4.pdf")


# In[20]:


c5 = root.TCanvas("c5", "c5", 800,800)
c5.Divide(2,2, 0.01, 0.01, 0)

c5.cd(1)
hAnaBarX.Draw()
c5.cd(2)
hAnaBarY.Draw()
c5.cd(3)
hAnaBarZ.Draw()
c5.cd(4)
hAnaBarT.Draw()


# In[21]:


c5.Draw()
c5.Print("plots/c5.pdf")


# In[22]:


c6 = root.TCanvas("c6", "c6", 800, 800)
c6.Divide(1,1, 0.01, 0.01, 0)

c6.cd(1)
hE1vsE2.Draw("COLZ")


# In[23]:


c6.Draw()
c6.Print("plots/c6.pdf")


# In[24]:


c7 = root.TCanvas("c7", "c7", 800, 800)
c7.Divide(2,2, 0.01, 0.01, 0)

c7PE_MeV = root.TCanvas("c7PE_MeV", "c7PE_MeV", 800,800)
c7Profile = root.TCanvas("c7Profile", "c7Profile", 800,800)

c7.cd(1)
hFinger_Edep_vs_Nphot.Draw("COLZ")
c7.cd(2)
hAnaBar_Edep_vs_Nphot.Draw("COLZ")
c7.cd(3)
hNphot0_vs_Nphot1.Draw("COLZ")

c7Profile.cd()
prof = hAnaBar_Edep_vs_Nphot.ProfileX()

prof.Fit("pol1")

c7PE_MeV.cd()
hAnaBar_Edep_vs_Nphot.Draw("COLZ")


# In[25]:


c7.Draw()
c7.Print("plots/c7.pdf")


# In[26]:


c7PE_MeV.Draw()
c7PE_MeV.Print("plots/c7PE_MeV.pdf")


# In[27]:


c7Profile.Draw()
c7Profile.Print("plots/c7Profile.pdf")


# In[28]:


c8 = root.TCanvas("c8", "c8", 800,800)
c8.Divide(2,2, 0.01, 0.01, 0)

cEd = root.TCanvas("cEd", "cEd", 800,800)
cEd.Divide(4,4)

cEdOne = root.TCanvas("cEdOne", "cEdOne", 800, 800)

c8.cd(1)
hAnaBarEdCut.Draw()
c8.cd(2)
hAnaBarPMTNphotCut[6].Draw()
c8.cd(3)
hAnaBar_Edep_vs_Nphot.Draw("COLZ")

means = []
meanErr = []

#function = root.TF1()

start = 5.5
gf = root.TF1("gf", "gaus", start, AnaBar_Edep_Max)

for i in range(NUMPADDLE):
    print("paddle ",i)
    if (i+1 == 3):
        cEdOne.cd()
    else:
        cEd.cd(i+1)
        
    hAnaBarEdAllCut[i].Draw()
    
    hAnaBarEdAllCut[i].Fit(gf, "R")
    function = hAnaBarEdAllCut[i].GetFunction("gf")
    function.SetLineColor(1)
    
    means.append(function.GetParameter(1))
    meanErr.append(function.GetParError(1))

for i in range(NUMPADDLE):
    print("Paddle " + str(i+1) + ": Mean peak Edep = " + str(means[i]) + " MeV")
    print("    \t Mean Edep error = " +str(meanErr[i]) + " MeV")
    
sumMeans = 0.0
sumMeanErrSqrs = 0.0

for i in range(NUMPADDLE):
    sumMeans += means[i]
    sumMeanErrSqrs += meanErr[i]*meanErr[i]

print("Sum of mean error squares = " + str(sumMeanErrSqrs))

meanMean = sumMeans/NUMPADDLE

sumErr = np.sqrt(sumMeanErrSqrs)
print("Error in sum of means = " + str(sumErr))
meanMeanErr = sumErr/NUMPADDLE

print("Mean peak Edep across all paddles: " + str(meanMean) + " MeV" )
print("Mean peak Edep uncertainty: " + str(meanMeanErr) + " Mev")


# In[29]:


c8.Draw()
c8.Print("plots/c8.pdf")


# In[30]:


#cEdOne.cd()
#start = 4.0
#gf = root.TF1("gf", "gaus", start, AnaBar_Edep_Max)
#hAnaBarEdAllCut[2].Fit(gf, "R")

cEdOne.Draw()
cEdOne.Print("plots/cEdOne.pdf")

#function = hAnaBarEdAllCut[2].GetFunction("gf")
#print(function.GetParameter(0))


# In[31]:


cEd.Draw()
cEd.Print("plots/cEd.pdf")


# In[32]:


c9 = root.TCanvas("c9", "c9", 800,800)
#c9.Divide(4,4, 0.01, 0.01, 0)

c9Single = root.TCanvas("c9Single", "c9Single", 800, 800)

print("Fittingi A1 ...\n")
fr = [float, float]
fp, fpe = [float, float, float, float], [float, float, float, float]
pllo = [0.05, 0.5, 1.0, 0.04]
plhi = [10.0, 50.0, 10000.0, 5.0]
sv = [1.8, 5.0, 1400.0, 3.0]
chisqr = float
ndf = int
SNRPeak, SNRFWHM = float, float


for i in range(NUMPADDLE):
    if(i == NUMPADDLE - 1):
        c9Single.cd()
        hAnaBarPMTNphot[i].Draw()
        hAnaBarPMTNoiseCutNphot[i].Draw("SAME")
        fr[0] = 0.7*hAnaBarPMTNphot[i].GetMean()
        fr[1] = 25.0*hAnaBarPMTNphot[i].GetMean()
    else:
        c9.cd()
        xl = 0.25*(i%4)
        xh = 0.25*(i%4)+0.25
        yl = 0.33*(i%3)
        yh = 0.33*(i%3)+0.33
        pad = root.TPad("pad","pad",xl,yl,xh,yh)
        pad.SetLogy(True)
        pad.Draw()
        pad.cd()
        hAnaBarPMTNphot[i].Draw()
        hAnaBarPMTNoiseCutNphot[i].Draw("SAME")
        
        #fr[0] = 0.7*hAnaBarPMTNphot[i].GetMean()
        #fr[1] = 25.0*hAnaBarPMTNphot[i].GetMean()
        #fitsnr = root.TF1(hAnaBarPMTNphot[i],fr,sv,pllo,plhi,fp,fpe,chisqr,ndf)
        #langaupro(fp,SNRPeak,SNRFWHM)


# In[33]:


c9.Draw()
c9.Print("plots/c9.pdf")


# In[34]:


c9Single.SetLogy()
c9Single.Draw()
c9Single.Print("plots/c9Single.pdf")


# In[35]:


c10 = root.TCanvas("c10", "c10", 800,800)
c10.Divide(1,2, 0.01, 0.01, 0)

c10.cd(1)
hFingerPMTKE.Draw()
c10.cd(2)
hAnaBarPMTKEA1.Draw()


# In[36]:


c10.Draw()
c10.Print("plots/c10.pdf")


# In[37]:


c11 = root.TCanvas("c11", "c11", 800,800)
c11.Divide(1,1, 0.01, 0.01, 0)

c11.cd(1)
hAnaBarMult.Draw()


# In[38]:


c11.Draw()
c11.Print("plots/c11.pdf")


# In[39]:


c12 = root.TCanvas("c12", "c12", 800,800)
c12.Divide(2,2, 0.01, 0.01, 0)

c12.cd(1)
hPrimPx.Draw()
c12.cd(2)
hPrimPy.Draw()
c12.cd(3)
hPrimPz.Draw()


# In[40]:


c12.Draw()
c12.Print("plots/c12.pdf")


# In[ ]:





# In[ ]:





# In[ ]:



