#include "TMath.h"

Double_t fitfunction(Double_t *x, Double_t *par) {
	Float_t xx = x[0];
	Float_t binsize = 1000.0/100.0;
	//Double_t f = binsize*par[0]/sqrt(2.0*3.14159265*pow(par[2],2))*exp(-0.5*(xx-par[1])*(xx-par[1])/(par[2]*par[2]));
	Double_t f = par[3]*xx*xx+par[4]*xx+par[5]+binsize*par[0]/sqrt(2.0*3.14159265*pow(par[2],2))*exp(-0.5*(xx-par[1])*(xx-par[1])/(par[2]*par[2]));
	return f;
}


void AnalyseSignals(Int_t Analysis_Run_Number = 1001, Int_t Analyse_Secondaries = 1, Float_t Theta_min_cut = 0.0, bool displayall = true) {

  //-------------------------------------------------------------------
  //Set stuff up for reading
  //-------------------------------------------------------------------
  TString filename;
  filename.Form("data/AnaBarMC_%d.root",Analysis_Run_Number);
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
  Float_t block_size = 4.0;
  Float_t pedthresh = 0.5;
  Int_t PMT_id;
  Int_t PMT_Nphotons;

  tree1->SetBranchAddress("Prim_E", &Prim_E);
  tree1->SetBranchAddress("Prim_Th", &Prim_Th);
  tree1->SetBranchAddress("Prim_Ph", &Prim_Ph);
  tree1->SetBranchAddress("Prim_pdg", &Prim_pdg);
  tree1->SetBranchAddress("PMT_Nphotons", &PMT_Nphotons);
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
  
  TH1F *hPrimE = new TH1F("PrimE","Primary Energy", 110, 0, 11000);
  TH1F *hPrimTh = new TH1F("PrimTh","Primary Theta", 100, 0, TMath::Pi());
  TH1F *hPrimPh = new TH1F("PrimPh","Primary Phi", 100, -TMath::Pi(), TMath::Pi());
  TH1F *hPrimPdg = new TH1F("PrimPdg","Primary PDG ID", 100, 0, 100);

  TH1F *hDetectorNhits = new TH1F("DetectorNhits","Detector Number of Hits", 400, 0, 400);
  TH1F *hDetectorPdg = new TH1F("DetectorPdg","Detector PDG ID", 200, -100, 100);
  TH1F *hDetectorID = new TH1F("DetectorID","Detector ID Number", 20, 0, 19);
  TH1F *hPMTID = new TH1F("PMTID","PMT ID Number", 5, 0, 5);
  TH1F *hAnaBarPMTNphot = new TH1F("AnaBarPMTNphot","AnaBar PMT Number of Photons", 299, 1, 300);
  
  TH2F *hAnaBar_Edep_vs_Nphot = new TH2F("AnaBarEdepVsNphot", "AnaBar Edep vs. Number of Photons", 299, 1, 299, 100, 0.01, 150);
  TH2F *hAnaBar_Nphot_vs_Eprimary = new TH2F("AnaBarNphotVsEprimary", "AnaBar Number of Photons vs. Eprimary", 110, 0, 11000, 299, 1, 300);
  
  TH1F *hAnaBarX = new TH1F("AnaBarX","AnaBar X Position", 100, -120, 120);
  TH1F *hAnaBarY = new TH1F("AnaBarY","AnaBar Y Position", 100, -30, 30);
  TH1F *hAnaBarZ = new TH1F("AnaBarZ","AnaBar Z Position", 100, -30, 30);
  TH1F *hAnaBarT = new TH1F("AnaBarT","AnaBar Time", 100, 0, 5);
  
  TH1F *hAnaBarEd = new TH1F("AnaBarEd","AnaBar Energy Central", 100, 0.01, 150);
  TH1F *hAnaBarEd_mu = new TH1F("AnaBarEd_mu","AnaBar Energy Middle Upper", 100, 0.01, 150);
  TH1F *hAnaBarEd_ml = new TH1F("AnaBarEd_ml","AnaBar Energy Middle Lower", 100, 0.01, 150);
  TH1F *hAnaBarEd_rm = new TH1F("AnaBarEd_rm","AnaBar Energy Right Middle", 100, 0.01, 150);
  TH1F *hAnaBarEd_lm = new TH1F("AnaBarEd_lm","AnaBar Energy Left Middle", 100, 0.01, 150);
  TH1F *hAnaBarEd_ru = new TH1F("AnaBarEd_ru","AnaBar Energy Right Upper", 100, 0.01, 150);
  TH1F *hAnaBarEd_lu = new TH1F("AnaBarEd_lu","AnaBar Energy Left Upper", 100, 0.01, 150);
  TH1F *hAnaBarEd_rl = new TH1F("AnaBarEd_rl","AnaBar Energy Right Lower", 100, 0.01, 150);
  TH1F *hAnaBarEd_ll = new TH1F("AnaBarEd_ll","AnaBar Energy Left Lower", 100, 0.01, 150);

  TH2F *hAnaBar_Edep_vs_Y = new TH2F("AnaBarEdepVsY", "AnaBar vs. Y entrant", 100, -30, 30, 100, 0.01, 150);

  TH2F *hyentran1_vs_xentran1 = new TH2F("Yentrance vs Xentrance", "Yentrance vs Xentrance", 100, -120, 120, 100, -30, 30);
  TH2F *hyexit1_vs_xexit1 = new TH2F("Yexit vs Xexit", "Yexit vs Xexit", 100, -120, 120, 100, -30, 30);
  TH2F *htracklength_vs_AnaBar_Edep = new TH2F("Tracklength vs AnaBar Edep", "Tracklength vs AnaBar Edep", 100, 0, 150, 100, 0.01 , 260);

  //
  // Histograms for position calculation
  //

  TH1F *xposition;
  TH1F *yposition;
  TH1F *eupdown;
  TH1F *eleftright;
  TH2F *xyposition;
  xposition = new TH1F("xposition","X Position",300,-6.0,6.0);
  yposition = new TH1F("yposition","Y Position",300,-6.0,6.0);
  eupdown = new TH1F("eupdown","E Up Down",300,-0.5,0.5);
  eleftright = new TH1F("eleftright","E Left Right",300,-0.5,0.5);
  xyposition = new TH2F("xyposition","XY Position",300,-6.0,6.0,300,-6.0,6.0);

  //-------------------------------------------------------------------
  //Limits and cuts
  //-------------------------------------------------------------------

  //-------------------------------------------------------------------
  //Event loop
  //-------------------------------------------------------------------
  
  Long64_t nentries = tree1->GetEntries();

  Long64_t counter = 0;
  float edep1tot[9];
  float n_good_edep = 0.0;
  
  for (Int_t i = 0; i < nentries; i++) {
    for (Int_t j = 0; j< 9; j++) {
	 edep1tot[j] = 0.0;
    }
    float yentrant1 = 0.0;
    float xentrant1 = 0.0;
    float zentrant1 = 0.0;	
    float yexit1 = 0.0;
    float xexit1 = -1000.0;
    float zexit1 = 0.0;
    float tracklength = 0.0;
    tree1->GetEntry(i);

    hPrimE->Fill(Prim_E);
    hPrimTh->Fill(Prim_Th);
    hPrimPh->Fill(Prim_Ph);
    hPrimPdg->Fill(Prim_pdg);
    hPMTID->Fill(PMT_id);
    
    hDetectorNhits->Fill(Detector_Nhits);

    bool trigger = false;
    bool anabar_hit = false;
    int j_anabar = 0;
    for (Int_t j=0; j < Detector_Nhits ; j++) {
	if (Detector_id[j] == 1 && !anabar_hit) {
		anabar_hit = true;
		j_anabar = j;
	}
    }

    if (anabar_hit) trigger = true; 
    if (trigger) {
  	if (PMT_Nphotons > 0) hAnaBarPMTNphot->Fill(PMT_Nphotons);
    }

    for (Int_t j=0; j < Detector_Nhits ; j++) {

	if (trigger) {
		counter++;
		if (Detector_id[j] >= 1 && Detector_id[j] <= 9 ) {
			if (j==j_anabar) {
	                   yentrant1 = Detector_y[j];
                           xentrant1 = Detector_x[j];
                           zentrant1 = Detector_z[j];

			}

        		hDetectorPdg->Fill(Detector_pdg[j]);
        		hDetectorID->Fill(Detector_id[j]);

        		if (Detector_pdg[j] == 22) {
        		  hAnaBarX->Fill(Detector_x[j]);
        		  hAnaBarY->Fill(Detector_y[j]);
        		  hAnaBarZ->Fill(Detector_z[j]);
    			  hAnaBarT->Fill(Detector_t[j]);
                        
				if(xexit1 < Detector_x[j]) {
				    zexit1 = Detector_z[j];
				    yexit1 = Detector_y[j];
				    xexit1 = Detector_x[j];
				}
			}
			if (Analyse_Secondaries == 1 && Prim_Th > Theta_min_cut) {
			  edep1tot[Detector_id[j]-1] += Detector_Ed[j];
			}else{ if (Detector_pdg[j] == 22 && Prim_Th > Theta_min_cut) {
				edep1tot[Detector_id[j]-1] += Detector_Ed[j];
			       }
			}

		}
	}
    }
    //cout << "Energy deposited = " << edep0tot << endl;

	tracklength = sqrt((xexit1-xentrant1)*(xexit1-xentrant1)+(yexit1-yentrant1)*(yexit1-yentrant1)+(zexit1-zentrant1)*(zexit1-zentrant1));
 
    if (trigger) {
        //std::cout << "Track Length = " << tracklength << std::endl;
        //std::cout << "Edep = " << edep1tot[0] << "  Number of photons = " << PMT_Nphotons << std::endl;
    	hAnaBarEd->Fill(edep1tot[0]);
        if (edep1tot[0] > pedthresh) n_good_edep=n_good_edep+1;
    	hAnaBarEd_rm->Fill(edep1tot[1]);
    	hAnaBarEd_lm->Fill(edep1tot[2]);
    	hAnaBarEd_mu->Fill(edep1tot[3]);
    	hAnaBarEd_ml->Fill(edep1tot[4]);
    	hAnaBarEd_ru->Fill(edep1tot[5]);
    	hAnaBarEd_lu->Fill(edep1tot[6]);
    	hAnaBarEd_rl->Fill(edep1tot[7]);
    	hAnaBarEd_ll->Fill(edep1tot[8]);
    	hAnaBar_Edep_vs_Y->Fill(yentrant1,edep1tot[0]);
    	if (PMT_Nphotons > 0) hAnaBar_Edep_vs_Nphot->Fill(PMT_Nphotons,edep1tot[0]);
    	if (PMT_Nphotons > 0) hAnaBar_Nphot_vs_Eprimary->Fill(Prim_E,PMT_Nphotons);
   	hyentran1_vs_xentran1->Fill(zentrant1,yentrant1);
 	htracklength_vs_AnaBar_Edep->Fill(edep1tot[0],tracklength);
	hyexit1_vs_xexit1->Fill(zexit1,yexit1);

	Float_t emid = edep1tot[0];
	Float_t eright = edep1tot[1];
	Float_t eleft = edep1tot[2];
 	Float_t eupper = edep1tot[3];
	Float_t elower = edep1tot[4];

	Float_t holdxpos, holdypos;

	Float_t sigmaxa,sigmaxb,xposcorr,xposoffset;
	Float_t sigmaya,sigmayb,yposcorr,yposoffset;

	xposcorr=1.0;
	xposoffset=0.0;
	yposcorr=1.0;
	yposoffset=0.0;
	sigmaxa=1.7;
	sigmaya=1.7;
	sigmaxb=1.0;
	sigmayb=1.0;

	//std::cout << "Eleft = " << eleft << "  Emid = " << emid << "  Eright = " << eright << std::endl;
        if(eleft < pedthresh && eright > pedthresh && emid > pedthresh) {
                eleftright->Fill((block_size/2.0-sigmaxa*log(0.5*(emid/eright)+1.0)));
                holdxpos=((block_size/2.0-sigmaxa*log(0.5*(emid/eright)+1.0))*xposcorr+xposoffset  );
        }else if (eleft > pedthresh && eright < pedthresh && emid > pedthresh) {
                eleftright->Fill((-1.0*block_size/2.0+sigmaxa*log(0.5*(emid/eleft)+1.0)));
                holdxpos=((-1.0*block_size/2.0+sigmaxa*log(0.5*(emid/eleft)+1.0))*xposcorr+xposoffset  );
        }else if (eleft > pedthresh && eright > pedthresh) {
                eleftright->Fill( (sigmaxb/2.0*log(eleft/eright))  );
                holdxpos=((sigmaxb/2.0*log(eleft/eright))*xposcorr+xposoffset  );
        }else{
                eleftright->Fill(0.0);
                holdxpos=0.0;
        }

        if(elower < pedthresh && eupper > pedthresh && emid > pedthresh) {
                eupdown->Fill((block_size/2.0-sigmaya*log(0.5*(emid/eupper)+1.0)));
                holdypos=((block_size/2.0-sigmaya*log(0.5*(emid/eupper)+1.0))*yposcorr+yposoffset  );
        }else if (elower > pedthresh && eupper < pedthresh && emid > pedthresh) {
                eupdown->Fill((-1.0*block_size/2.0+sigmaya*log(0.5*(emid/elower)+1.0)));
                holdypos=((-1.0*block_size/2.0+sigmaya*log(0.5*(emid/elower)+1.0))*yposcorr+yposoffset  );
        }else if (elower > pedthresh && eupper > pedthresh) {
                eupdown->Fill( (sigmayb/2.0*log(elower/eupper))  );
                holdypos=((sigmayb/2.0*log(elower/eupper))*yposcorr+yposoffset  );
        }else{
                eupdown->Fill(0.0);
                holdypos=0.0;
        }

        xposition->Fill(holdxpos);
        yposition->Fill(holdypos);
        xyposition->Fill(holdxpos,holdypos);

     }
  }
  
  //-------------------------------------------------------------------
  //Plotting and writing out
  //-------------------------------------------------------------------
  
  Double_t n_primaries = (Double_t) hPrimE->GetEntries();
  Double_t n_photons = (Double_t) hAnaBarPMTNphot->GetEntries();
  //Double_t eff = n_photons/n_primaries*100.0;
  Double_t eff = n_good_edep/n_primaries*100.0;
  Double_t deff = 1.0/sqrt(n_good_edep)*eff;
  std::cout << "Efficiency = " << eff << " +/- " << deff << " %" << std::endl;
  
  if (displayall) {
  TCanvas *c1 = new TCanvas("c1", "c1", 100,100,650,650);
  c1->Divide(3,3, 0.01, 0.01, 0);
  TCanvas *c2 = new TCanvas("c2", "c2", 150,150,650,450);
  c2->Divide(2,2, 0.01, 0.01, 0);
  TCanvas *c3 = new TCanvas("c3", "c3", 700,100,600,400);
  c3->Divide(2,2, 0.01, 0.01, 0);
  TCanvas *c4 = new TCanvas("c4", "c4", 100,500,600,400);
  c4->Divide(3,1, 0.01, 0.01, 0);
  TCanvas *c5 = new TCanvas("c5", "c5", 700,500,600,400);
  c5->Divide(2,2, 0.01, 0.01, 0);
  TCanvas *c6 = new TCanvas("c6", "c6", 1300,500,600,400);
  c6->Divide(2,2, 0.01, 0.01, 0);
  TCanvas *c7 = new TCanvas("c7", "c7", 200,500,600,400);
  c7->Divide(2,1, 0.01, 0.01, 0);
  TCanvas *c8 = new TCanvas("c8", "c8", 250,550,600,400);
  c8->Divide(2,3, 0.01, 0.01, 0);

  c8->cd(1);
  gPad->SetLogy();
  xposition->Draw();
  c8->cd(2);
  gPad->SetLogy();
  yposition->Draw();
  c8->cd(3);
  xyposition->Draw("COLZ");
  c8->cd(4);
  gPad->SetLogy();
  eupdown->Draw();
  c8->cd(5);
  gPad->SetLogy();
  eleftright->Draw();

  c1->cd(1);
  gPad->SetLogy();
  hAnaBarEd_lu->Draw();
  c1->cd(2);
  gPad->SetLogy();
  hAnaBarEd_mu->Draw();
  c1->cd(3);
  gPad->SetLogy();
  hAnaBarEd_ru->Draw();
  c1->cd(4);
  gPad->SetLogy();
  hAnaBarEd_lm->Draw();
  c1->cd(5);
  gPad->SetLogy();
  hAnaBarEd->Draw();
  c1->cd(6);
  gPad->SetLogy();
  hAnaBarEd_rm->Draw();
  c1->cd(7);
  gPad->SetLogy();
  hAnaBarEd_ll->Draw();
  c1->cd(8);
  gPad->SetLogy();
  hAnaBarEd_ml->Draw();
  c1->cd(9);
  gPad->SetLogy();
  hAnaBarEd_rl->Draw();
  
  c2->cd(1);
  hPrimE->Draw();
  c2->cd(2);
  hPrimTh->Draw();
  c2->cd(3);
  hPrimPh->Draw();
  c2->cd(4);
  hPrimPdg->Draw();
  
  c3->cd(1);
  gPad->SetLogy();
  hDetectorNhits->Draw();
  c3->cd(2);
  gPad->SetLogy();
  hDetectorPdg->Draw();
  c3->cd(3);
  hDetectorID->Draw();
  c3->cd(4);
  hPMTID->Draw();
  
  c5->cd(1);
  hAnaBarX->Draw();
  c5->cd(2);
  hAnaBarY->Draw();
  c5->cd(3);
  hAnaBarZ->Draw();
  c5->cd(4);
  hAnaBarT->Draw();


  c7->cd(1);
  hAnaBar_Edep_vs_Nphot->Draw("COLZ");
  c7->cd(2);
  hAnaBar_Nphot_vs_Eprimary->Draw("COLZ");

  c6->cd(1);
 gPad->SetLogz();
  hyentran1_vs_xentran1->Draw("COLZ");
  c6->cd(2);
 gPad->SetLogz();
  hyexit1_vs_xexit1->Draw("COLZ");
  c6->cd(3);
  htracklength_vs_AnaBar_Edep->Draw("COLZ");
  }

  
  c4->cd(1);
  gPad->SetLogz();
  hAnaBar_Edep_vs_Y->Draw("COLZ");
  c4->cd(2);
  hAnaBarEd->Draw();
  c4->cd(3);
  hAnaBarPMTNphot->Draw();



}
