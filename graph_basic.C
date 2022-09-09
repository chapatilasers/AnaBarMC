//
// Fit to U.S. Standard Atmosphere Data
// From http://www.engineeringtoolbox.com/standard-atmosphere-d_604.html
//

#include <ostream>

Double_t fitfunction(Double_t *x, Double_t *par)
{
	Double_t f = par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0);
	//Double_t f = par[0]*exp(par[1]*x[0]*x[0]);
	return f;
}

void graph_basic() {
   //
   //Modified: Edward J. Brash
   //		
  
   //gStyle->SetOptFit(kFALSE);	
   gStyle->SetOptFit(1);	
   TCanvas *c1 = new TCanvas("c1","Absorption and Emission vs. Wavelength",200,10,700,500);
   FILE *f1 = fopen("emission_scint.dat","r");
   FILE *f2 = fopen("emission_wls.dat","r");
   FILE *f3 = fopen("absorption_wls.dat","r");

   c1->SetFillColor(42);
   c1->SetGrid();

   const Int_t n = 100;
   Double_t e1[n], l1[n], em1[n];
   Double_t e2[n], l2[n], em2[n];
   Double_t e3[n], l3[n], em3[n];
   Int_t i=0;

   while (!feof(f1)){
	  fscanf(f1,"%lf %lf \n",&e1[i],&em1[i]);
          l1[i]=1240.0/e1[i];
	  i++;
   }

   const Int_t n1=i;
   i=0;

   Double_t wls_emax = 0.0;
   while (!feof(f2)){
	  fscanf(f2,"%lf %lf \n",&e2[i],&em2[i]);
          l2[i]=1240.0/e2[i];
	  if (em2[i] > wls_emax) wls_emax = em2[i];
	  i++;
   }


   const Int_t n2=i;
   i=0;

   for (int ii=0; ii<n2; ii++) {em2[ii]=em2[ii]/wls_emax;}

   while (!feof(f3)){
	  fscanf(f3,"%lf %lf \n",&e3[i],&em3[i]);
          l3[i]=1240.0/e3[i];
	  em3[i]=em3[i]/1.4;
	  i++;
   }

   const Int_t n3=i;
   i=0;


//   gr1->SetLineWidth(4);
   gr1 = new TGraph(n1,l1,em1);
   gr1->SetMarkerColor(kBlue);
   gr1->SetLineColor(kBlue);
   gr1->SetMarkerStyle(21);
   gr1->SetTitle("Scintillator Emission");
   gr1->GetXaxis()->SetTitle("Wavelength (nm)");
   gr1->GetYaxis()->SetTitle("Emission (Relative)");
   gr1->Draw("ALP");
   
   gr2 = new TGraph(n2,l2,em2);
   gr2->SetMarkerColor(kRed);
   gr2->SetLineColor(kRed);
   gr2->SetMarkerStyle(21);
   gr2->SetTitle("WLS Emission");
   gr2->GetXaxis()->SetTitle("Wavelength (nm)");
   gr2->GetYaxis()->SetTitle("Emission (Relative)");
   gr2->Draw("LP");

   gr3 = new TGraph(n3,l3,em3);
   gr3->SetMarkerColor(kGreen);
   gr3->SetLineColor(kGreen);
   gr3->SetMarkerStyle(21);
   gr3->SetTitle("WLS Absorption");
   gr3->GetXaxis()->SetTitle("Wavelength (nm)");
   gr3->GetYaxis()->SetTitle("Absorption (Relative)");
   gr3->Draw("LP");

   // TCanvas::Update() draws the frame, after which one can change it
   c1->Update();
   c1->GetFrame()->SetFillColor(31);
   c1->GetFrame()->SetBorderSize(12);
   c1->Modified();
}
