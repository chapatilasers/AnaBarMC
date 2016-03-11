#include "TMath.h"

Double_t fitfunction(Double_t *x, Double_t *par) {
	Float_t xx = x[0];
	Float_t binsize = 1000.0/100.0;
	//Double_t f = binsize*par[0]/sqrt(2.0*3.14159265*pow(par[2],2))*exp(-0.5*(xx-par[1])*(xx-par[1])/(par[2]*par[2]));
	Double_t f = par[3]*xx*xx+par[4]*xx+par[5]+binsize*par[0]/sqrt(2.0*3.14159265*pow(par[2],2))*exp(-0.5*(xx-par[1])*(xx-par[1])/(par[2]*par[2]));
	return f;
}


void AnalyseSignals(Int_t Analysis_Run_Number = 1, Int_t Analyse_Secondaries = 1, Float_t Theta_min_cut = 0.0) {

  //-------------------------------------------------------------------
  //Set stuff up for reading
  //-------------------------------------------------------------------
  TString filename;
  filename.Form("data/G4_test%d.root",Analysis_Run_Number);
  TFile *f1 = new TFile(filename,"READ");
  TTree *tree1 = (TTree*)f1->Get("T");

  const int MaxHits = 10000;
  Float_t Prim_E, Prim_Th, Prim_Ph;
  Int_t Prim_pdg;
  Int_t Detector_Nhits;
  Int_t Detector_pdg[MaxHits];
  Int_t Detector_id[MaxHits];
  Float_t Detector_x[MaxHits], Detector_y[MaxHits], Detector_z[MaxHits], Detector_t[MaxHits];
  Float_t Detector_Ed[MaxHits];  
  Int_t PMT_id;
  Int_t PMT_Nphotons_Zero;
  Int_t PMT_Nphotons_One;

  tree1->SetBranchAddress("Prim_E", &Prim_E);
  tree1->SetBranchAddress("Prim_Th", &Prim_Th);
  tree1->SetBranchAddress("Prim_Ph", &Prim_Ph);
  tree1->SetBranchAddress("Prim_pdg", &Prim_pdg);
  tree1->SetBranchAddress("PMT_Nphotons_Zero", &PMT_Nphotons_Zero);
  tree1->SetBranchAddress("PMT_Nphotons_One", &PMT_Nphotons_One);
  tree1->SetBranchAddress("PMT_id", &PMT_id);
  tree1->SetBranchAddress("Detector_Nhits", &Detector_Nhits);
  tree1->SetBranchAddress("Detector_pdg", &Detector_pdg);
  tree1->SetBranchAddress("Detector_id", &Detector_id);
  tree1->SetBranchAddress("Detector_x", &Detector_x);
  tree1->SetBranchAddress("Detector_y", &Detector_y);
  tree1->SetBranchAddress("Detector_z", &Detector_z);
  tree1->SetBranchAddress("Detector_t", &Detector_t);
  tree1->SetBranchAddress("Detector_Ed", &Detector_Ed);

  //-------------------------------------------------------------------
  //Create histograms
  //-------------------------------------------------------------------
  
  TH1F *hPrimE = new TH1F("PrimE","Primary Energy", 100, 0, 25000);
  TH1F *hPrimTh = new TH1F("PrimTh","Primary Theta", 100, 0, TMath::Pi());
  TH1F *hPrimPh = new TH1F("PrimPh","Primary Phi", 100, -TMath::Pi(), TMath::Pi());
  TH1F *hPrimPdg = new TH1F("PrimPdg","Primary PDG ID", 20, 0, 20);

  TH1F *hDetectorNhits = new TH1F("DetectorNhits","Detector Number of Hits", 100, 0, 400);
  TH1F *hDetectorPdg = new TH1F("DetectorPdg","Detector PDG ID", 50, -20, 30);
  TH1F *hDetectorID = new TH1F("DetectorID","Detector ID Number", 5, 0, 5);
  TH1F *hPMTID = new TH1F("PMTID","PMT ID Number", 5, 0, 5);
  TH1F *hFingerPMTNphot = new TH1F("FingerPMTNphot","Finger PMT Number of Photons", 200, 0, 200);
  TH1F *hAnaBarPMTNphot = new TH1F("AnaBarPMTNphot","AnaBar PMT Number of Photons", 200, 0, 200);
  
  TH2F *hFinger_Edep_vs_Nphot = new TH2F("FingerEdepVsNphot", "Finger Edep vs. Number of Photons", 200, 0, 200, 100, 0.01, 10);
  TH2F *hAnaBar_Edep_vs_Nphot = new TH2F("AnaBarEdepVsNphot", "AnaBar Edep vs. Number of Photons", 200, 0, 200, 100, 0.01, 30);
  
  TH1F *hFingerX = new TH1F("FingerX","Finger X Position", 100, -120, 120);
  TH1F *hFingerY = new TH1F("FingerY","Finger Y Position", 100, -30, 30);
  TH1F *hFingerZ = new TH1F("FingerZ","Finger Z Position", 100, -30, 30);
  TH1F *hFingerT = new TH1F("FingerT","Finger Time", 100, 0, 5);
  TH1F *hAnaBarX = new TH1F("AnaBarX","AnaBar X Position", 100, -120, 120);
  TH1F *hAnaBarY = new TH1F("AnaBarY","AnaBar Y Position", 100, -30, 30);
  TH1F *hAnaBarZ = new TH1F("AnaBarZ","AnaBar Z Position", 100, -30, 30);
  TH1F *hAnaBarT = new TH1F("AnaBarT","AnaBar Time", 100, 0, 5);
  
  TH1F *hFingerEd = new TH1F("FingerEd","Finger Energy Deposited", 100, 0.01, 10);
  TH1F *hAnaBarEd = new TH1F("AnaBarEd","AnaBar Energy Deposited", 100, 0.01, 30);

  TH2F *hFinger_Edep_vs_Y = new TH2F("FingerEdepVsY", "Finger Edep vs. Y entrant", 100, -30, 30, 100, 0.01, 10);
  TH2F *hAnaBar_Edep_vs_Y = new TH2F("AnaBarEdepVsY", "AnaBar vs. Y entrant", 100, -30, 30, 100, 0.01, 30);
  TH2F *hE1vsE2 = new TH2F("E1vsE2", "AnaBar Edep vs. Finger Edep", 100, 0.01, 10, 100, 0.01, 30);

  TH2F *hyentran1_vs_xentran1 = new TH2F("Yentrance vs Xentrance", "Yentrance vs Xentrance", 100, -80, 80, 100, -30, 30);
  TH2F *hyexit1_vs_xexit1 = new TH2F("Yexit vs Xexit", "Yexit vs Xexit", 100, -80, 80, 100, -30, 30);
  TH2F *htracklength_vs_AnaBar_Edep = new TH2F("Tracklength vs AnaBar Edep", "Tracklength vs AnaBar Edep", 100, -2, 10, 100, 0.01 , 15);
  //-------------------------------------------------------------------
  //Limits and cuts
  //-------------------------------------------------------------------

  //-------------------------------------------------------------------
  //Event loop
  //-------------------------------------------------------------------
  
  Long64_t nentries = tree1->GetEntries();

  Long64_t counter = 0;
  
  for (Int_t i = 0; i < nentries; i++) {
    float edep0tot = 0.0;
    float edep1tot = 0.0;
    float yentrant0 = 0.0;
    float yentrant1 = 0.0;
    float xentrant0 = 0.0;
    float xentrant1 = 0.0;
    float zentrant0 = 0.0;	
    float zentrant1 = 0.0;	
    float yexit0 = 0.0;
    float yexit1 = 0.0;
    float xexit0 = 0.0;
    float xexit1 = 0.0;
    float zexit0 = 1000.0;
    float zexit1 = 1000.0;
    float tracklength = 0.0;
    tree1->GetEntry(i);

    hPrimE->Fill(Prim_E);
    hPrimTh->Fill(Prim_Th);
    hPrimPh->Fill(Prim_Ph);
    hPrimPdg->Fill(Prim_pdg);
    hPMTID->Fill(PMT_id);
    
    hDetectorNhits->Fill(Detector_Nhits);

    bool trigger = false;
    bool finger_hit = false;
    bool anabar_hit = false;
    int j_finger = 0;
    int j_anabar = 0;
    for (Int_t j=0; j < Detector_Nhits ; j++) {
	if (Detector_id[j] == 0 && !finger_hit) {
		finger_hit = true;
		j_finger = j;
		//cout<<"hit in finger";
	}
	if (Detector_id[j] == 1 && !anabar_hit) {
		anabar_hit = true;
		j_anabar = j;
	}
    }

    if (finger_hit && anabar_hit) trigger = true; 
    if (trigger) {
  	hAnaBarPMTNphot->Fill(PMT_Nphotons_Zero);
        hFingerPMTNphot->Fill(PMT_Nphotons_One);
    }

    for (Int_t j=0; j < Detector_Nhits ; j++) {

	if (trigger) {
		
		counter++;
		if (Detector_id[j] == 0 ) {
			if (j==j_finger) {
			   yentrant0 = Detector_y[j];
		 	   xentrant0 = Detector_x[j];
		 	   zentrant0 = Detector_z[j];
			}	

        		hDetectorPdg->Fill(Detector_pdg[j]);
        		hDetectorID->Fill(Detector_id[j]);

        		if (Detector_pdg[j] == 13) {
			  hFingerX->Fill(Detector_x[j]);
        		  hFingerY->Fill(Detector_y[j]);
        		  hFingerZ->Fill(Detector_z[j]);
    			  hFingerT->Fill(Detector_t[j]);

				if(zexit0 > Detector_z[j]) {
				    zexit0 = Detector_z[j];
				    yexit0 = Detector_y[j];
				    xexit0 = Detector_x[j];
				}
			
                        }
			if (Analyse_Secondaries == 1 && Prim_Th > Theta_min_cut) {
			  edep0tot += Detector_Ed[j];
			}else{ if (Detector_pdg[j] == 13 && Prim_Th > Theta_min_cut) {
				edep0tot += Detector_Ed[j];
			       }
			} 

		}
		if (Detector_id[j] == 1 ) {
			if (j==j_anabar) {
	                   yentrant1 = Detector_y[j];
                           xentrant1 = Detector_x[j];
                           zentrant1 = Detector_z[j];

			}

        		hDetectorPdg->Fill(Detector_pdg[j]);
        		hDetectorID->Fill(Detector_id[j]);

        		if (Detector_pdg[j] == 13) {
        		  hAnaBarX->Fill(Detector_x[j]);
        		  hAnaBarY->Fill(Detector_y[j]);
        		  hAnaBarZ->Fill(Detector_z[j]);
    			  hAnaBarT->Fill(Detector_t[j]);
                        
				if(zexit1 > Detector_z[j]) {
				    zexit1 = Detector_z[j];
				    yexit1 = Detector_y[j];
				    xexit1 = Detector_x[j];
				}
			}
			if (Analyse_Secondaries == 1 && Prim_Th > Theta_min_cut) {
			  edep1tot += Detector_Ed[j];
			}else{ if (Detector_pdg[j] == 13 && Prim_Th > Theta_min_cut) {
				edep1tot += Detector_Ed[j];
			       }
			}

		}
	}
    }
    //cout << "Energy deposited = " << edep0tot << endl;

	tracklength = sqrt((xexit1-xentrant1)*(xexit1-xentrant1)+(yexit1-yentrant1)*(yexit1-yentrant1)+(zexit1-zentrant1)*(zexit1-zentrant1));
 
    if (trigger) {
    	hFingerEd->Fill(edep0tot);
    	hAnaBarEd->Fill(edep1tot);
    	hFinger_Edep_vs_Y->Fill(yentrant0,edep0tot);
    	hAnaBar_Edep_vs_Y->Fill(yentrant1,edep1tot);
    	hFinger_Edep_vs_Nphot->Fill(PMT_Nphotons_One,edep0tot);
    	hAnaBar_Edep_vs_Nphot->Fill(PMT_Nphotons_Zero,edep1tot);
    	hE1vsE2->Fill(edep0tot,edep1tot);
   	hyentran1_vs_xentran1->Fill(xentrant1,yentrant1);
 	htracklength_vs_AnaBar_Edep->Fill(edep1tot,tracklength);
	hyexit1_vs_xexit1->Fill(xexit1,yexit1);
  
     }
  }
  
  //-------------------------------------------------------------------
  //Plotting and writing out
  //-------------------------------------------------------------------
  
  TCanvas *c2 = new TCanvas("c2", "c2", 100,100,600,400);
  c2->Divide(2,2, 0.01, 0.01, 0);
  TCanvas *c3 = new TCanvas("c3", "c3", 700,100,600,400);
  c3->Divide(2,2, 0.01, 0.01, 0);
  
  c2->cd(1);
  hPrimE->Draw();
  c2->cd(2);
  hPrimTh->Draw();
  c2->cd(3);
  hPrimPh->Draw();
  c2->cd(4);
  hPrimPdg->Draw();
  
  c3->cd(1);
  hDetectorNhits->Draw();
  c3->cd(2);
  gPad->SetLogy();
  hDetectorPdg->Draw();
  c3->cd(3);
  hDetectorID->Draw();
  c3->cd(4);
  hPMTID->Draw();

  TCanvas *c1 = new TCanvas("c1", "c1", 1300,100,600,400);
  c1->Divide(2,2, 0.01, 0.01, 0);
  TCanvas *c5 = new TCanvas("c5", "c5", 700,500,600,400);
  c5->Divide(2,2, 0.01, 0.01, 0);
  TCanvas *c4 = new TCanvas("c4", "c4", 100,500,600,400);
  c4->Divide(3,2, 0.01, 0.01, 0);
  TCanvas *c6 = new TCanvas("c6", "c6", 1300,500,600,400);
  c6->Divide(2,2, 0.01, 0.01, 0);
  TCanvas *c7 = new TCanvas("c7", "c7", 200,500,600,400);
  c7->Divide(2,1, 0.01, 0.01, 0);
 
  c1->cd(1);
  hFingerX->Draw();
  c1->cd(2);
  hFingerY->Draw();
  c1->cd(3);
  hFingerZ->Draw();
  c1->cd(4);
  hFingerT->Draw();
  
  c5->cd(1);
  hAnaBarX->Draw();
  c5->cd(2);
  hAnaBarY->Draw();
  c5->cd(3);
  hAnaBarZ->Draw();
  c5->cd(4);
  hAnaBarT->Draw();

  c4->cd(1);
  gPad->SetLogz();
  hFinger_Edep_vs_Y->Draw("COLZ");
  c4->cd(2);
  hFingerEd->Draw();
  c4->cd(3);
  hFingerPMTNphot->Draw();
  c4->cd(4);
  gPad->SetLogz();
  hAnaBar_Edep_vs_Y->Draw("COLZ");
  c4->cd(5);
  hAnaBarEd->Draw();
  c4->cd(6);
  hAnaBarPMTNphot->Draw();

  c7->cd(1);
  hFinger_Edep_vs_Nphot->Draw("COLZ");
  c7->cd(2);
  hAnaBar_Edep_vs_Nphot->Draw("COLZ");

  c6->cd(1);
  gPad->SetLogz();
  hE1vsE2->Draw("COLZ");
  c6->cd(2);
 gPad->SetLogz();
  hyentran1_vs_xentran1->Draw("COLZ");
  c6->cd(3);
 gPad->SetLogz();
  hyexit1_vs_xexit1->Draw("COLZ");
  c6->cd(4);
  htracklength_vs_AnaBar_Edep->Draw("COLZ");

}
