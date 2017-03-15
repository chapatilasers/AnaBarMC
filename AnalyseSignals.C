#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"

Double_t langaufun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation), 
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      mpc = par[1] - mpshift * par[0]; 

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
}

TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
   // Once again, here are the Landau * Gaussian parameters:
   //   par[0]=Width (scale) parameter of Landau density
   //   par[1]=Most Probable (MP, location) parameter of Landau density
   //   par[2]=Total area (integral -inf to inf, normalization constant)
   //   par[3]=Width (sigma) of convoluted Gaussian function
   //
   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[4]  reasonable start values for the fit
   //   parlimitslo[4]  lower parameter limits
   //   parlimitshi[4]  upper parameter limits
   //   fitparams[4]    returns the final fit parameters
   //   fiterrors[4]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf

   Int_t i;
   Char_t FunName[100];

   sprintf(FunName,"Fitfcn_%s",his->GetName());

   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;

   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MP","Area","GSigma");
   
   for (i=0; i<4; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }

   his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot

   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<4; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf

   return (ffit);              // return fit function

}

Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {

   // Seaches for the location (x value) at the maximum of the 
   // Landau-Gaussian convolute and its full width at half-maximum.
   //
   // The search is probably not very efficient, but it's a first try.

   Double_t p,x,fy,fxr,fxl;
   Double_t step;
   Double_t l,lold;
   Int_t i = 0;
   Int_t MAXCALLS = 10000;


   // Search for maximum

   p = params[1] - 0.1 * params[0];
   step = 0.05 * params[0];
   lold = -2.0;
   l    = -1.0;


   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = langaufun(&x,params);
 
      if (l < lold)
         step = -step/10;
 
      p += step;
   }

   if (i == MAXCALLS)
      return (-1);

   maxx = x;

   fy = l/2;


   // Search for right x location of fy

   p = maxx + params[0];
   step = params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;


   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
 
      if (l > lold)
         step = -step/10;
 
      p += step;
   }

   if (i == MAXCALLS)
      return (-2);

   fxr = x;


   // Search for left x location of fy

   p = maxx - 0.5 * params[0];
   step = -params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;

   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
 
      if (l > lold)
         step = -step/10;
 
      p += step;
   }

   if (i == MAXCALLS)
      return (-3);


   fxl = x;

   FWHM = fxr - fxl;
   return (0);
}

void AnalyseSignals(Int_t Analysis_Run_Number = 1, Int_t Analyse_Secondaries = 1, Float_t Theta_min_cut = 0.0, bool displayall = false) {

  //-------------------------------------------------------------------
  //Set stuff up for reading
  //-------------------------------------------------------------------
  TString filename;
  filename.Form("data/AnaBarMC_%d.root",Analysis_Run_Number);
  TFile *f1 = new TFile(filename,"READ");
  TTree *tree1 = (TTree*)f1->Get("T");

  const int MaxHits = 10000;
  const int MaxPMTNo = 50;
  const int MaxPMTHits = 500;
  const Float_t Finger_Edep_Max = 10.0;
  const Float_t AnaBar_Edep_Max = 5.0;
  const Int_t Finger_NPhotons_Max = 1000;
  const Int_t AnaBar_NPhotons_Max = 30;
  Float_t Prim_E, Prim_Th, Prim_Ph;
  Int_t Prim_pdg;
  Int_t Detector_Nhits;
  Int_t Detector_pdg[MaxHits];
  Int_t Detector_id[MaxHits];
  Float_t Detector_x[MaxHits], Detector_y[MaxHits], Detector_z[MaxHits], Detector_t[MaxHits];
  Float_t Detector_Ed[MaxHits];  
  Int_t PMT_id;
  Int_t PMT_Nphotons[MaxPMTNo];
  Float_t PMT_KineticEnergy[MaxPMTNo][MaxPMTHits];

  tree1->SetBranchAddress("Prim_E", &Prim_E);
  tree1->SetBranchAddress("Prim_Th", &Prim_Th);
  tree1->SetBranchAddress("Prim_Ph", &Prim_Ph);
  tree1->SetBranchAddress("Prim_pdg", &Prim_pdg);
  tree1->SetBranchAddress("PMT_id", &PMT_id);
  tree1->SetBranchAddress("PMT_Nphotons", &PMT_Nphotons);
  tree1->SetBranchAddress("PMT_KineticEnergy", &PMT_KineticEnergy);
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
  TH1F *hFingerPMTNphot = new TH1F("FingerPMTNphot","Finger PMT Number of Photons", Finger_NPhotons_Max, 0, Finger_NPhotons_Max);
  TH1F *hAnaBarPMTNphotA1 = new TH1F("AnaBarPMTNphotA1","AnaBar PMT Number of Photons A1", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max);
  TH1F *hAnaBarPMTNphotA2 = new TH1F("AnaBarPMTNphotA2","AnaBar PMT Number of Photons A2", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max);
  TH1F *hAnaBarPMTNphotA3 = new TH1F("AnaBarPMTNphotA3","AnaBar PMT Number of Photons A3", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max);
  TH1F *hAnaBarPMTNphotA4 = new TH1F("AnaBarPMTNphotA4","AnaBar PMT Number of Photons A4", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max);
  TH1F *hAnaBarPMTNphotA5 = new TH1F("AnaBarPMTNphotA5","AnaBar PMT Number of Photons A5", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max);
  TH1F *hAnaBarPMTNphotA6 = new TH1F("AnaBarPMTNphotA6","AnaBar PMT Number of Photons A6", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max);
  TH1F *hAnaBarPMTNphotA7 = new TH1F("AnaBarPMTNphotA7","AnaBar PMT Number of Photons A7", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max);
  TH1F *hAnaBarPMTNphotA8 = new TH1F("AnaBarPMTNphotA8","AnaBar PMT Number of Photons A8", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max);
  TH1F *hAnaBarPMTNphotA9 = new TH1F("AnaBarPMTNphotA9","AnaBar PMT Number of Photons A9", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max);
  TH1F *hAnaBarPMTNphotA10 = new TH1F("AnaBarPMTNphotA10","AnaBar PMT Number of Photons A10", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max);
  TH1F *hAnaBarPMTNphotA11 = new TH1F("AnaBarPMTNphotA11","AnaBar PMT Number of Photons A11", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max);
  TH1F *hAnaBarPMTNphotA12 = new TH1F("AnaBarPMTNphotA12","AnaBar PMT Number of Photons A12", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max);
  TH1F *hAnaBarPMTNphotA13 = new TH1F("AnaBarPMTNphotA13","AnaBar PMT Number of Photons A13", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max);
  TH1F *hAnaBarPMTNphotA14 = new TH1F("AnaBarPMTNphotA14","AnaBar PMT Number of Photons A14", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max);
  TH1F *hFingerPMTKE = new TH1F("FingerPMTKE","Photon Wavelength Production Spectrum", 400, 300.0, 700.0);
  TH1F *hAnaBarPMTKEA1 = new TH1F("AnaBarPMTKEA1","Photon Wavelength in WLS at PMT", 400, 300.0, 700.0);
  
  TH2F *hFinger_Edep_vs_Nphot = new TH2F("FingerEdepVsNphot", "Finger Edep vs. Number of Photons", Finger_NPhotons_Max, 0, Finger_NPhotons_Max, 100, 0.01, Finger_Edep_Max);
  TH2F *hAnaBar_Edep_vs_Nphot = new TH2F("AnaBarEdepVsNphot", "AnaBar Edep vs. Number of Photons", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max, 100, 0.01, 5);
  TH2F *hNphot0_vs_Nphot1 = new TH2F("AnaBarVsFingerNphot", "AnaBar vs. Finger Number of Photons", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max, Finger_NPhotons_Max, 0, Finger_NPhotons_Max);
  
  TH1F *hFingerX = new TH1F("FingerX","Finger X Position", 100, -120, 120);
  TH1F *hFingerY = new TH1F("FingerY","Finger Y Position", 100, -30, 30);
  TH1F *hFingerZ = new TH1F("FingerZ","Finger Z Position", 100, -30, 30);
  TH1F *hFingerT = new TH1F("FingerT","Finger Time", 100, 0, .4);
  TH1F *hAnaBarX = new TH1F("AnaBarX","AnaBar X Position", 100, -120, 120);
  TH1F *hAnaBarY = new TH1F("AnaBarY","AnaBar Y Position", 100, -30, 30);
  TH1F *hAnaBarZ = new TH1F("AnaBarZ","AnaBar Z Position", 100, -30, 30);
  TH1F *hAnaBarT = new TH1F("AnaBarT","AnaBar Time", 100, 0, .4);
  
  TH1F *hFingerEd = new TH1F("FingerEd","Finger Energy Deposited", 100, 0.01, Finger_Edep_Max);
  TH1F *hAnaBarEd = new TH1F("AnaBarEd","AnaBar Energy Deposited", 100, 0.01, AnaBar_Edep_Max);

  TH2F *hFinger_Edep_vs_Y = new TH2F("FingerEdepVsY", "Finger Edep vs. Y entrant", 100, -30, 30, 100, 0.01, Finger_Edep_Max);
  TH2F *hAnaBar_Edep_vs_Y = new TH2F("AnaBarEdepVsY", "AnaBar Edep vs. Y entrant", 100, -30, 30, 100, 0.01, AnaBar_Edep_Max);
  TH2F *hE1vsE2 = new TH2F("E1vsE2", "AnaBar Edep vs. Finger Edep", 100, 0.01, Finger_Edep_Max, 100, 0.01, AnaBar_Edep_Max);

  TH2F *hyentran1_vs_xentran1 = new TH2F("Yentrance vs Xentrance", "Yentrance vs Xentrance", 100, -80, 80, 100, -30, 30);
  TH2F *hyexit1_vs_xexit1 = new TH2F("Yexit vs Xexit", "Yexit vs Xexit", 100, -80, 80, 100, -30, 30);
  TH2F *htracklength_vs_AnaBar_Edep = new TH2F("Tracklength vs AnaBar Edep", "Tracklength vs AnaBar Edep", 100, -2, AnaBar_Edep_Max, 100, -2 , 10);
  //-------------------------------------------------------------------
  //Limits and cuts
  //-------------------------------------------------------------------

  //-------------------------------------------------------------------
  //Event loop
  //-------------------------------------------------------------------
  
  Long64_t nentries = tree1->GetEntries();

  Long64_t counter = 0;

  const int NMaxPMT=14;
  float edeptot[NMaxPMT];
  
  for (Int_t i = 0; i < nentries; i++) {
    for (Int_t j=0; j<NMaxPMT; j++) {edeptot[j] = 0.0;};
    float edep0tot = 0.0;
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
    bool anabar_top_hit = false;
    bool anabar_bottom_hit = false;
    int j_finger = 0;
    int j_anabar = 0;
    for (Int_t j=0; j < Detector_Nhits ; j++) {
	if (Detector_id[j] == 0 && !finger_hit) {
		finger_hit = true;
		j_finger = j;
		//cout<<"hit in finger";
	}
	if (Detector_id[j] == 1 && !anabar_top_hit) {
		anabar_top_hit = true;
		j_anabar = j;
	}
	if (Detector_id[j] == 14 && !anabar_bottom_hit) {
		anabar_bottom_hit = true;
		j_anabar = j;
	}
    }

    if (finger_hit && anabar_top_hit && anabar_bottom_hit) trigger = true; 
    //if (finger_hit && anabar_top_hit) trigger = true; 
    //if (finger_hit) trigger = true; 
    
    if (trigger) {
        //for (Int_t j=0; j<15; j++) { std::cout << "j = " << j << " Nphotons = " << PMT_Nphotons[j] << std::endl; }
  	hAnaBarPMTNphotA1->Fill(PMT_Nphotons[0]);
  	hAnaBarPMTNphotA2->Fill(PMT_Nphotons[1]);
  	hAnaBarPMTNphotA3->Fill(PMT_Nphotons[2]);
  	hAnaBarPMTNphotA4->Fill(PMT_Nphotons[3]);
  	hAnaBarPMTNphotA5->Fill(PMT_Nphotons[4]);
  	hAnaBarPMTNphotA6->Fill(PMT_Nphotons[5]);
  	hAnaBarPMTNphotA7->Fill(PMT_Nphotons[6]);
  	hAnaBarPMTNphotA8->Fill(PMT_Nphotons[7]);
  	hAnaBarPMTNphotA9->Fill(PMT_Nphotons[8]);
  	hAnaBarPMTNphotA10->Fill(PMT_Nphotons[9]);
  	hAnaBarPMTNphotA11->Fill(PMT_Nphotons[10]);
  	hAnaBarPMTNphotA12->Fill(PMT_Nphotons[11]);
  	hAnaBarPMTNphotA13->Fill(PMT_Nphotons[12]);
  	hAnaBarPMTNphotA14->Fill(PMT_Nphotons[13]);
        hFingerPMTNphot->Fill(PMT_Nphotons[14]);
	for (Int_t jq=0; jq<PMT_Nphotons[14]; jq++){
		//std::cout << "Processing Finger hit = " << jq << std::endl; 
		hFingerPMTKE->Fill(1240.0/PMT_KineticEnergy[14][jq]);
	}
	for (Int_t iq=1; iq<14; iq++){
		for (Int_t jq=0; jq<PMT_Nphotons[iq]; jq++){
			//std::cout << "Processing Anabar pmt = " << iq << " hit = " << jq << " Energy = " << PMT_KineticEnergy[iq][jq] << std::endl;
			hAnaBarPMTKEA1->Fill(1240.0/PMT_KineticEnergy[iq][jq]);
		}
	}
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

		}

		if (Detector_id[j] > 0 && Detector_id[j] <= NMaxPMT) {
			if (Analyse_Secondaries == 1 && Prim_Th > Theta_min_cut) {
				edeptot[Detector_id[j]-1] += Detector_Ed[j];
			}else{ if (Detector_pdg[j] == 13 && Prim_Th > Theta_min_cut) {
					edeptot[Detector_id[j]-1] += Detector_Ed[j];
		     	       }
			}
		}
	}
    }
    //cout << "Energy deposited = " << edep0tot << endl;

	tracklength = sqrt((xexit1-xentrant1)*(xexit1-xentrant1)+(yexit1-yentrant1)*(yexit1-yentrant1)+(zexit1-zentrant1)*(zexit1-zentrant1));
 
    if (trigger) {
    	hFingerEd->Fill(edep0tot);
    	hAnaBarEd->Fill(edeptot[0]);
    	hFinger_Edep_vs_Y->Fill(yentrant0,edep0tot);
    	hAnaBar_Edep_vs_Y->Fill(yentrant1,edeptot[0]);
    	hFinger_Edep_vs_Nphot->Fill(PMT_Nphotons[14],edep0tot);
    	hAnaBar_Edep_vs_Nphot->Fill(PMT_Nphotons[0],edeptot[0]);
    	hNphot0_vs_Nphot1->Fill(PMT_Nphotons[0],PMT_Nphotons[14]);
    	hE1vsE2->Fill(edep0tot,edeptot[0]);
   	hyentran1_vs_xentran1->Fill(xentrant1,yentrant1);
 	htracklength_vs_AnaBar_Edep->Fill(edeptot[0],tracklength);
	hyexit1_vs_xexit1->Fill(xexit1,yexit1);
  
     }
  }
  
  //-------------------------------------------------------------------
  //Plotting and writing out
  //-------------------------------------------------------------------
  
  
  if (displayall) {
  TCanvas *c1 = new TCanvas("c1", "c1", 100,100,500,270);
  c1->Divide(2,2, 0.01, 0.01, 0);
  TCanvas *c2 = new TCanvas("c2", "c2", 600,100,500,270);
  c2->Divide(2,2, 0.01, 0.01, 0);
  TCanvas *c3 = new TCanvas("c3", "c3", 1100,100,500,270);
  c3->Divide(2,2, 0.01, 0.01, 0);
  TCanvas *c5 = new TCanvas("c5", "c5", 100,400,500,270);
  c5->Divide(2,2, 0.01, 0.01, 0);
  TCanvas *c6 = new TCanvas("c6", "c6", 1100,400,500,270);
  c6->Divide(2,2, 0.01, 0.01, 0);
  TCanvas *c7 = new TCanvas("c7", "c7", 100,700,500,270);
  c7->Divide(2,2, 0.01, 0.01, 0);
  TCanvas *c8 = new TCanvas("c8", "c8", 600,700,500,270);
  c8->Divide(2,2, 0.01, 0.01, 0);
  TCanvas *c9 = new TCanvas("c9", "c9", 1100,700,300,270);
  c9->Divide(4,4, 0.01, 0.01, 0);
  TCanvas *c10 = new TCanvas("c10", "c10", 1100,700,300,270);
  c10->Divide(1,2, 0.01, 0.01, 0);
  
  c1->cd(1);
  hFingerX->Draw();
  c1->cd(2);
  hFingerY->Draw();
  c1->cd(3);
  hFingerZ->Draw();
  c1->cd(4);
  hFingerT->Draw();
  
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
  
  c5->cd(1);
  hAnaBarX->Draw();
  c5->cd(2);
  hAnaBarY->Draw();
  c5->cd(3);
  hAnaBarZ->Draw();
  c5->cd(4);
  hAnaBarT->Draw();
  
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

  c7->cd(1);
  hFinger_Edep_vs_Nphot->Draw("COLZ");
  c7->cd(2);
  hAnaBar_Edep_vs_Nphot->Draw("COLZ");
  c7->cd(3);
  hNphot0_vs_Nphot1->Draw("COLZ");

  c8->cd(1);
  hAnaBarEd->Draw();
  c8->cd(2);
  hAnaBarPMTNphotA1->Draw();
  c8->cd(3);
  hAnaBar_Edep_vs_Nphot->Draw("COLZ");
  c8->cd(4);
  htracklength_vs_AnaBar_Edep->Draw("COLZ");
  
  printf("Fittingi A1 ...\n");
  // Setting fit range and start values
  Double_t fr[2];
  Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
  pllo[0]=0.05; pllo[1]=0.5; pllo[2]=1.0; pllo[3]=0.04;
  plhi[0]=10.0; plhi[1]=50.0; plhi[2]=10000.0; plhi[3]=5.0;
  sv[0]=1.8; sv[1]=5.0; sv[2]=1400.0; sv[3]=3.0;
  Double_t chisqr;
  Int_t    ndf;
  Double_t SNRPeak, SNRFWHM;
  
  c9->cd(1);
  hAnaBarPMTNphotA1->Draw();
  fr[0]=0.1*hAnaBarPMTNphotA1->GetMean();
  fr[1]=4.0*hAnaBarPMTNphotA1->GetMean();
  TF1 *fitsnr = langaufit(hAnaBarPMTNphotA1,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  langaupro(fp,SNRPeak,SNRFWHM);
  fitsnr->Draw("SAME");
  
  c9->cd(2);
  hAnaBarPMTNphotA2->Draw();
  fr[0]=0.1*hAnaBarPMTNphotA2->GetMean();
  fr[1]=4.0*hAnaBarPMTNphotA2->GetMean();
  TF1 *fitsnr = langaufit(hAnaBarPMTNphotA2,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  langaupro(fp,SNRPeak,SNRFWHM);
  fitsnr->Draw("SAME");
  
  c9->cd(3);
  hAnaBarPMTNphotA3->Draw();
  fr[0]=0.1*hAnaBarPMTNphotA3->GetMean();
  fr[1]=4.0*hAnaBarPMTNphotA3->GetMean();
  TF1 *fitsnr = langaufit(hAnaBarPMTNphotA3,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  langaupro(fp,SNRPeak,SNRFWHM);
  fitsnr->Draw("SAME");
  
  c9->cd(4);
  hAnaBarPMTNphotA4->Draw();
  fr[0]=0.1*hAnaBarPMTNphotA4->GetMean();
  fr[1]=4.0*hAnaBarPMTNphotA4->GetMean();
  TF1 *fitsnr = langaufit(hAnaBarPMTNphotA4,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  langaupro(fp,SNRPeak,SNRFWHM);
  fitsnr->Draw("SAME");
  
  c9->cd(5);
  hAnaBarPMTNphotA5->Draw();
  fr[0]=0.1*hAnaBarPMTNphotA5->GetMean();
  fr[1]=4.0*hAnaBarPMTNphotA5->GetMean();
  TF1 *fitsnr = langaufit(hAnaBarPMTNphotA5,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  langaupro(fp,SNRPeak,SNRFWHM);
  fitsnr->Draw("SAME");
  
  c9->cd(6);
  hAnaBarPMTNphotA6->Draw();
  fr[0]=0.1*hAnaBarPMTNphotA6->GetMean();
  fr[1]=4.0*hAnaBarPMTNphotA6->GetMean();
  TF1 *fitsnr = langaufit(hAnaBarPMTNphotA6,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  langaupro(fp,SNRPeak,SNRFWHM);
  fitsnr->Draw("SAME");
  
  c9->cd(7);
  hAnaBarPMTNphotA7->Draw();
  fr[0]=0.1*hAnaBarPMTNphotA7->GetMean();
  fr[1]=4.0*hAnaBarPMTNphotA7->GetMean();
  TF1 *fitsnr = langaufit(hAnaBarPMTNphotA7,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  langaupro(fp,SNRPeak,SNRFWHM);
  fitsnr->Draw("SAME");
  
  c9->cd(8);
  hAnaBarPMTNphotA8->Draw();
  fr[0]=0.1*hAnaBarPMTNphotA8->GetMean();
  fr[1]=4.0*hAnaBarPMTNphotA8->GetMean();
  TF1 *fitsnr = langaufit(hAnaBarPMTNphotA8,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  langaupro(fp,SNRPeak,SNRFWHM);
  fitsnr->Draw("SAME");
  
  c9->cd(9);
  hAnaBarPMTNphotA9->Draw();
  fr[0]=0.1*hAnaBarPMTNphotA9->GetMean();
  fr[1]=4.0*hAnaBarPMTNphotA9->GetMean();
  TF1 *fitsnr = langaufit(hAnaBarPMTNphotA9,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  langaupro(fp,SNRPeak,SNRFWHM);
  fitsnr->Draw("SAME");
  
  c9->cd(10);
  hAnaBarPMTNphotA10->Draw();
  fr[0]=0.1*hAnaBarPMTNphotA10->GetMean();
  fr[1]=4.0*hAnaBarPMTNphotA10->GetMean();
  TF1 *fitsnr = langaufit(hAnaBarPMTNphotA10,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  langaupro(fp,SNRPeak,SNRFWHM);
  fitsnr->Draw("SAME");
  
  c9->cd(11);
  hAnaBarPMTNphotA11->Draw();
  fr[0]=0.1*hAnaBarPMTNphotA11->GetMean();
  fr[1]=4.0*hAnaBarPMTNphotA11->GetMean();
  TF1 *fitsnr = langaufit(hAnaBarPMTNphotA11,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  langaupro(fp,SNRPeak,SNRFWHM);
  fitsnr->Draw("SAME");
  
  c9->cd(12);
  hAnaBarPMTNphotA12->Draw();
  fr[0]=0.1*hAnaBarPMTNphotA12->GetMean();
  fr[1]=4.0*hAnaBarPMTNphotA12->GetMean();
  TF1 *fitsnr = langaufit(hAnaBarPMTNphotA12,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  langaupro(fp,SNRPeak,SNRFWHM);
  fitsnr->Draw("SAME");
  
  c9->cd(13);
  hAnaBarPMTNphotA13->Draw();
  fr[0]=0.1*hAnaBarPMTNphotA13->GetMean();
  fr[1]=4.0*hAnaBarPMTNphotA13->GetMean();
  TF1 *fitsnr = langaufit(hAnaBarPMTNphotA13,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  langaupro(fp,SNRPeak,SNRFWHM);
  fitsnr->Draw("SAME");
  
  c9->cd(14);
  hAnaBarPMTNphotA14->Draw();
  fr[0]=0.1*hAnaBarPMTNphotA14->GetMean();
  fr[1]=4.0*hAnaBarPMTNphotA14->GetMean();
  TF1 *fitsnr = langaufit(hAnaBarPMTNphotA14,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  langaupro(fp,SNRPeak,SNRFWHM);
  fitsnr->Draw("SAME");
  

  c10->cd(1);
  hFingerPMTKE->Draw();
  c10->cd(2);
  hAnaBarPMTKEA1->Draw();
  
  }

  
  TCanvas *c4 = new TCanvas("c4", "c4", 600,400,500,270);
  c4->Divide(3,2, 0.01, 0.01, 0);
  
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
  hAnaBarPMTNphotA1->Draw();

}


void langaus() {
   // Fill Histogram
   Int_t data[100] = {0,0,0,0,0,0,2,6,11,18,18,55,90,141,255,323,454,563,681,
                    737,821,796,832,720,637,558,519,460,357,291,279,241,212,
                    153,164,139,106,95,91,76,80,80,59,58,51,30,49,23,35,28,23,
                    22,27,27,24,20,16,17,14,20,12,12,13,10,17,7,6,12,6,12,4,
                    9,9,10,3,4,5,2,4,1,5,5,1,7,1,6,3,3,3,4,5,4,4,2,2,7,2,4};
   TH1F *hSNR = new TH1F("snr","Signal-to-noise",400,0,400);

   for (Int_t i=0; i<100; i++) hSNR->Fill(i,data[i]);

}

