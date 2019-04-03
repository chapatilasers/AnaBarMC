

void PlotMeanNumPhotonRatios(){

  TCanvas *grc = new TCanvas("grc", "grc", 100,100,700,500);
  grc->SetGrid();
  grc->SetFillColor(18);


  TGraph *plot_R = new TGraph("meanNumPhotonRatios");
  plot_R->SetTitle("Ratio of Mean Num. Photons (With Mirror / Without Mirror)");
  plot_R->SetMinimum(1.0); // sets minimum y value of graph display
  plot_R->SetMaximum(1.6); // sets maximum y value of graph display
  plot_R->SetMarkerStyle(20);
  plot_R->SetMarkerColor(4);
  plot_R->SetMarkerSize(1.4);
  plot_R->Draw("ALP");
}
