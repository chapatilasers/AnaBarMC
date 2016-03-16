#include <iostream>

void graph_eff() {
   //
   //Modified: Edward J. Brash
   //		
  
   //gStyle->SetOptFit(kFALSE);	
   gStyle->SetOptFit(1);	
   TCanvas *c1 = new TCanvas("c1","AnaBar Neutron Efficiency",200,10,700,500);
   c1->SetFillColor(42);
   c1->SetGrid();
   
   FILE *f = fopen("eff.dat","r");

   const Int_t n = 10000;
   Double_t x[n], y[n], ey[n], ex[n];
   Int_t i=0;

   while (!feof(f)){
	  fscanf(f,"%lf %lf %lf\n",&x[i],&y[i],&ey[i]);
	  ex[i]=0.0;
	  i++;
   }
   
   const Int_t n = i;
   gr = new TGraphErrors(n,x,y,ex,ey);
//   gr->SetLineColor(2);
//   gr->SetLineWidth(4);
//   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->SetTitle("AnaBar Neutron Efficiency");
   gr->GetXaxis()->SetTitle("X Position");
   gr->GetYaxis()->SetTitle("Efficiency (%)");
//
   gr->Draw("AP");

   // TCanvas::Update() draws the frame, after which one can change it
   c1->Update();
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(12);
   c1->Modified();
}
