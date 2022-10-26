#include <iostream>


using namespace std;
using RNode = ROOT::RDF::RNode;

int Analyse_Secondaries = 1;
float Theta_min_cut = 2.524;
float ThetaVerticalCut = 3.02;

int MaxHits = 50000;
int MaxPMTNo = 20;
int MaxPMTHits = 5000;
float Finger_Edep_Max = 10.0;
float AnaBar_Edep_Max = 10.0;
float pedastel_sigma = 2.9;
int Detector_Offset = 0;

int Finger_NPhotons_Max = 250;
int AnaBar_NPhotons_Max = 200;

int NUMPADDLE = 14;

int NMaxPMT = 14;

bool getTrigger(int Detector_Nhits, int* Detector_id) {

    bool fhit = false;
    bool ahit = false;
    bool trigger = false;
    for (int j=0; j<Detector_Nhits; j++) {
        if (Detector_id[j] == Detector_Offset && !fhit) {
            fhit = true;
        }
        for (int ibar=1; ibar<15; ibar++){
            if (Detector_id[j+Detector_Offset] == ibar+Detector_Offset) {
                ahit = true;
            }
        }
        if (fhit && ahit) {
            trigger = true;
        }
    }
    return trigger;
}

bool getTrigger2(bool trigger, float fNewTheta) {

    bool trigger2 = false;
    if (fNewTheta > Theta_min_cut) {
        trigger2 = true;
    }
    return trigger2;
}

bool getTrigger3(bool trigger, float fNewTheta) {

    bool trigger3 = false;
    if (fNewTheta > ThetaVerticalCut) {
        trigger3 = true;
    }
    return trigger3;
}

float getMass(int Prim_pdg) {

    float fMass;
    if (Prim_pdg == 11) {
        fMass = 0.511;
    } else {
        if (Prim_pdg == 13) {
            fMass = 105.7;
        } else {
            if (Prim_pdg == 2212) {
                fMass = 938.28;
            } else {
                fMass = 939.65;
            }
        }
    }
    return fMass;
}

float getMomentum(float Prim_E, float fMass) {
    return sqrt(Prim_E*Prim_E - fMass*fMass);
}

float getPx(float fMomentum, float Prim_Th, float Prim_Ph) {
    return fMomentum*TMath::Sin(Prim_Th)*TMath::Cos(Prim_Ph);
}

float getPy(float fMomentum, float Prim_Th, float Prim_Ph) {
    return fMomentum*TMath::Sin(Prim_Th)*TMath::Sin(Prim_Ph);
}

float getPz(float fMomentum, float Prim_Th, float Prim_Ph) {
    return fMomentum*TMath::Cos(Prim_Th);
}

float getNewTheta(float fMomentum, float fPy) {
    return TMath::ACos(fPy/fMomentum);
}

float getNewPhi(float fMomentum, float fPx, float fPz) {
    float fNewPhi;
    if (fPx < 0.0) {
        fNewPhi = TMath::ATan(fPz/fPx) + TMath::Pi();
    } else {
        if (fPx > 0.0 && fPz < 0.0) {
            fNewPhi = TMath::ATan(fPz/fPx) + TMath::TwoPi();
        } else {
            fNewPhi = TMath::ATan(fPz/fPx);
        }
    }
    return fNewPhi;
}

int getAnaBarMult(bool trigger, int* PMT_Nphotons) {

    int imult = 0;
    float temp;
    TRandom3* fRand = new TRandom3(-1);
    
    for (int icount = 0;icount < NUMPADDLE;icount++){
        temp = PMT_Nphotons[icount]+fRand->Gaus(0.0,pedastel_sigma);
        if (temp>3.0*pedastel_sigma) {
            imult++;
        }
    }
    
    return imult;
}

std::vector<float> getFingerPMTKE(bool trigger, int* PMT_Nphotons, float* PMT_KineticEnergy) {
    
    std::vector<float> v;
    
    //for (Int_t jq=0; jq<PMT_Nphotons[14]; jq++){
    for (Int_t jq=0; jq<1; jq++){
        if(PMT_KineticEnergy[14*5000+jq] > 0.001) {
            v.push_back(1240.0/PMT_KineticEnergy[14*5000+jq]);
        };
    }
    
    return v;
}

std::vector<float> getAnaBarPMTKE(bool trigger, int* PMT_Nphotons, float* PMT_KineticEnergy) {
    
    std::vector<float> v;
    
    for (Int_t iq=0; iq<NUMPADDLE; iq++) {
        //for (Int_t jq=0; jq<PMT_Nphotons[iq]; jq++){
        for (Int_t jq=0; jq<1; jq++){
            if(PMT_KineticEnergy[iq*5000+jq] > 0.001) {
                v.push_back(1240.0/PMT_KineticEnergy[iq*5000+jq]);
            };
        }
    }
    
    return v;
}

std::vector<float> getFingerXVec(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg, float* Detector_x, int Prim_pdg) {

    std::vector<float> v;

    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] == Detector_Offset && Detector_pdg[j] == Prim_pdg) {
                v.push_back(Detector_x[j]);
            }
        }
    }
    return v;
}
std::vector<float> getFingerYVec(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg, float* Detector_y, int Prim_pdg) {
    
    std::vector<float> v;

    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] == Detector_Offset && Detector_pdg[j] == Prim_pdg) {
                v.push_back(Detector_y[j]);
            }
        }
    }
    return v;
}
std::vector<float> getFingerZVec(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg, float* Detector_z, int Prim_pdg) {
    
    std::vector<float> v;
    
    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] == Detector_Offset && Detector_pdg[j] == Prim_pdg) {
                v.push_back(Detector_z[j]);
            }
        }
    }
    return v;
}
std::vector<float> getFingerTVec(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg, float* Detector_t, int Prim_pdg) {
    
    std::vector<float> v;
    
    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] == Detector_Offset && Detector_pdg[j] == Prim_pdg) {
                v.push_back(Detector_t[j]);
            }
        }
    }
    return v;
}

std::vector<float> getAnaBarXVec(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg, float* Detector_x, int Prim_pdg) {
    
    std::vector<float> v;
    
    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] == 1 + Detector_Offset && Detector_pdg[j] == Prim_pdg) {
                v.push_back(Detector_x[j]);
            }
        }
    }
    return v;
}
std::vector<float> getAnaBarYVec(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg, float* Detector_y, int Prim_pdg) {
    
    std::vector<float> v;
    
    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] == 1 + Detector_Offset && Detector_pdg[j] == Prim_pdg) {
                v.push_back(Detector_y[j]);
            }
        }
    }
    return v;
}
std::vector<float> getAnaBarZVec(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg, float* Detector_z, int Prim_pdg) {
    
    std::vector<float> v;
    
    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] == 1 + Detector_Offset && Detector_pdg[j] == Prim_pdg) {
                v.push_back(Detector_z[j]);
            }
        }
    }
    return v;
}
std::vector<float> getAnaBarTVec(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg, float* Detector_t, int Prim_pdg) {
   
    std::vector<float> v;
    
    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] == 1 + Detector_Offset && Detector_pdg[j] == Prim_pdg) {
                v.push_back(Detector_t[j]);
            }
        }
    }
    return v;
}

std::vector<int> getFingerID(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg) {
    
    std::vector<int> v;
    
    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] == Detector_Offset) {
                v.push_back(Detector_id[j]);
            }
        }
    }
    return v;
}

std::vector<int> getFingerPDG(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg) {
    
    std::vector<int> v;
    
    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] == Detector_Offset) {
                v.push_back(Detector_pdg[j]);
            }
        }
    }
    return v;
}

std::vector<int> getAnaBarID(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg) {

    std::vector<int> v;
    
    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] == 1 + Detector_Offset) {
                v.push_back(Detector_id[j]);
            }
        }
    }
    return v;
}

std::vector<int> getAnaBarPDG(bool trigger, int Detector_Nhits, int* Detector_id, int* Detector_pdg) {
    
    std::vector<int> v;
   
    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] == 1 + Detector_Offset) {
                v.push_back(Detector_pdg[j]);
            }
        }
    }
    return v;
}

std::vector<float> getFingerPMTNPhotons(bool trigger, int* PMT_Nphotons) {
    
    std::vector<float> v;
    TRandom3* fRand = new TRandom3(-1);
    float pmt0tot = 0;
    

    if (trigger) {
        for (Int_t icount = 14;icount<15;icount++){
            pmt0tot += PMT_Nphotons[icount]+fRand->Gaus(0.0,pedastel_sigma);
        }
    }
    
    v.push_back(pmt0tot);
    
    return v;
}

std::vector<float> getAnaBarPMTNPhotons(bool trigger, int* PMT_Nphotons) {
    
    std::vector<float> v;
    TRandom3* fRand = new TRandom3(-1);
    float pmttot[NMaxPMT];

    if (trigger) {
        for (Int_t icount = 0;icount<14;icount++){
            pmttot[icount] = PMT_Nphotons[icount]+fRand->Gaus(0.0,pedastel_sigma);
            v.push_back(pmttot[icount]);
        }
    }
    
    return v;
}

std::vector<float> getAnaBarNPhotonsTotal(bool trigger, int* PMT_Nphotons) {
    
    std::vector<float> v;
    TRandom3* fRand = new TRandom3(-1);
    float pmttot = 0;

    if (trigger) {
        for (Int_t icount = 0;icount<14;icount++){
            pmttot = pmttot + PMT_Nphotons[icount]+fRand->Gaus(0.0,pedastel_sigma);
        }
    }
    
    v.push_back(pmttot);
    
    return v;
}

std::vector<float> getFingerEd(bool trigger, float fNewTheta, int Detector_Nhits, int Prim_pdg, int* Detector_id, int* Detector_pdg, float* Detector_Ed) {
    
    std::vector<float> v;
    float edep0tot = 0;
    
    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] == Detector_Offset) {
                if (Analyse_Secondaries == 1 && fNewTheta > Theta_min_cut) {
                    edep0tot += Detector_Ed[j];
                }
            }
        }
    }
    
    v.push_back(edep0tot);
    
    return v;
}

std::vector<float> getAnaBarEd(bool trigger, float fNewTheta, int Detector_Nhits, int Prim_pdg, int* Detector_id, int* Detector_pdg, float* Detector_Ed) {
    
    std::vector<float> v;
    float edeptot[NMaxPMT];
    
    for (int j=0; j<NMaxPMT; j++) {
        edeptot[j]=0.0;
    }
    
    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] > Detector_Offset && Detector_id[j] <= NMaxPMT+Detector_Offset) {
                if (Analyse_Secondaries == 1 && fNewTheta > Theta_min_cut) {
                    edeptot[Detector_id[j]-1-Detector_Offset] += Detector_Ed[j];
                } else {
                    if (Detector_pdg[j] == Prim_pdg && fNewTheta > Theta_min_cut) {
                        edeptot[Detector_id[j]-1-Detector_Offset] += Detector_Ed[j];
                    }
                }
            }
        }
    }
    
    for (int j=0; j<NMaxPMT; j++) {
        v.push_back(edeptot[j]);
    }
    
    return v;
}

std::vector<float> getAnaBarEdTotal(bool trigger, float fNewTheta, int Detector_Nhits, int Prim_pdg, int* Detector_id, int* Detector_pdg, float* Detector_Ed) {
    
    std::vector<float> v;
    float edeptotal = 0;
    float edeptot[NMaxPMT];
    
    for (int j=0; j<NMaxPMT; j++) {
        edeptot[j]=0.0;
    }
    
    for (int j=0; j < Detector_Nhits; j++) {
        if (trigger) {
            if (Detector_id[j] > Detector_Offset && Detector_id[j] <= NMaxPMT+Detector_Offset) {
                if (Analyse_Secondaries == 1 && fNewTheta > Theta_min_cut) {
                    edeptot[Detector_id[j]-1-Detector_Offset] += Detector_Ed[j];
                } else {
                    if (Detector_pdg[j] == Prim_pdg && fNewTheta > Theta_min_cut) {
                        edeptot[Detector_id[j]-1-Detector_Offset] += Detector_Ed[j];
                    }
                }
            }
        }
    }
    
    for (int j=0; j<NMaxPMT; j++) {
        edeptotal += edeptot[j];
    }
    
    v.push_back(edeptotal);
    
    return v;
}

RNode AnalyseSignalsRDataFrame() {

	auto fileName = "data/AnaBarMC_777777.root";
	auto treeName = "T";

	ROOT::RDataFrame d(treeName,fileName);

	//auto entries = d.Count();
	//cout << *entries << " entries in Tree with no filter" << endl;

	auto fdf = d.Define("trigger", "getTrigger(Detector_Nhits, &Detector_id[0])")
       			.Define("fMass", "getMass(Prim_pdg)")
       			.Define("fMomentum","getMomentum(Prim_E,fMass)")
       			.Define("fPx", "getPx(fMomentum,Prim_Th,Prim_Ph)")
       			.Define("fPy", "getPy(fMomentum,Prim_Th,Prim_Ph)")
       			.Define("fPz", "getPz(fMomentum,Prim_Th,Prim_Ph)")
       			.Define("fNewTheta", "getNewTheta(fMomentum,fPy)")
       			.Define("fNewPhi", "getNewPhi(fMomentum,fPx,fPz)")
       			.Define("trigger2", "getTrigger2(trigger,fNewTheta)")
       			.Define("trigger3", "getTrigger3(trigger,fNewTheta)")
       			.Define("fingerXVec","getFingerXVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_x[0],Prim_pdg)")
       			.Define("fingerYVec","getFingerYVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_y[0],Prim_pdg)")
       			.Define("fingerZVec","getFingerZVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_z[0],Prim_pdg)")
       			.Define("fingerTVec","getFingerTVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_t[0],Prim_pdg)")
       			.Define("anaBarXVec","getAnaBarXVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_x[0],Prim_pdg)")
       			.Define("anaBarYVec","getAnaBarYVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_y[0],Prim_pdg)")
       			.Define("anaBarZVec","getAnaBarZVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_z[0],Prim_pdg)")
       			.Define("anaBarTVec","getAnaBarTVec(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0],&Detector_t[0],Prim_pdg)")
       			.Define("fingerID","getFingerID(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])")
       			.Define("fingerPDG","getFingerPDG(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])")
       			.Define("anaBarID","getAnaBarID(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])")
       			.Define("anaBarPDG","getAnaBarPDG(trigger,Detector_Nhits,&Detector_id[0],&Detector_pdg[0])")
       			.Define("fingerPMTNPhotons","getFingerPMTNPhotons(trigger,&PMT_Nphotons[0])")
       			.Define("anaBarPMTNPhotons","getAnaBarPMTNPhotons(trigger,&PMT_Nphotons[0])")
       			.Define("anaBarNPhotonsTotal","getAnaBarNPhotonsTotal(trigger,&PMT_Nphotons[0])")
       			.Define("imult","getAnaBarMult(trigger,&PMT_Nphotons[0])")
       			.Define("fingerPMTKE","getFingerPMTKE(trigger,&PMT_Nphotons[0],&PMT_KineticEnergy[0])")
      			.Define("anaBarPMTKE","getAnaBarPMTKE(trigger,&PMT_Nphotons[0],&PMT_KineticEnergy[0])")
       			.Define("fingerEd","getFingerEd(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])")
       			.Define("anaBarEd","getAnaBarEd(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])")
       			.Define("anaBarEdTotal","getAnaBarEdTotal(trigger,fNewTheta,Detector_Nhits,Prim_pdg,&Detector_id[0],&Detector_pdg[0],&Detector_Ed[0])");

	//auto entries2 = fdf.Count();
	//cout << *entries2 << " entries in Expanded Dataframe with no filter" << endl;

	auto triggers = fdf.Filter("trigger==true").Count();
	cout << *triggers << " entries passed Main trigger" << endl;

	auto fdft = fdf.Filter("trigger==true");

	return fdft;
}

TCanvas* plotC1(){

  RNode fdft = AnalyseSignalsRDataFrame();

  auto hFingerX = fdft.Histo1D("fingerXVec");
  auto hFingerY = fdft.Histo1D("fingerYVec");
  auto hFingerZ = fdft.Histo1D("fingerZVec");
  auto hFingerT = fdft.Histo1D("fingerTVec");
  
  TCanvas *c1 = new TCanvas("c1", "c1", 100,100,500,270);
  c1->Divide(2,2, 0.01, 0.01, 0);

  c1->cd(1);
  hFingerX->Draw();
  c1->cd(2);
  hFingerY->Draw();
  c1->cd(3);
  hFingerZ->Draw();
  c1->cd(4);
  hFingerT->Draw();

  c1->DrawClone();
  c1->Print("plots/c1.pdf");

  return c1;

}
