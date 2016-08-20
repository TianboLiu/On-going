{
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(14);
  
  tole = sqrt(29.7);

  TCanvas *c1 = new TCanvas("c1","c1",1200,500);
  c1->Divide(3,1);

  c1->cd(1);
 
  gPad->SetRightMargin(0.01);
  gPad->SetLeftMargin(0.05);
  TH2F* h1 = new TH2F("h1","",100, 0.17, 0.63, 100, -3.5,-0.);
  h1->GetYaxis()->SetLabelSize(0);
  h1->GetYaxis()->SetNdivisions(0);
  h1->GetXaxis()->SetNdivisions(4);
  h1->GetXaxis()->SetLabelSize(0.08);
  h1->GetXaxis()->SetLabelOffset(-0.0);
  h1->Draw();

  TLatex* tu = new TLatex();
  tu->SetNDC();
  tu->SetTextSize(0.1);
  tu->DrawLatex(0.4,0.78,"#delta_{T} u^{[0,1]}");

  TGraphAsymmErrors *gre_u = new TGraphAsymmErrors(1);
  gre_u->SetMarkerStyle(8);
  gre_u->SetMarkerSize(2.0);
  gre_u->SetMarkerColor(1);
  gre_u->SetLineColor(1);
  gre_u->SetLineWidth(2);

  //#0: Kang 2015
  gre_u->SetPoint(0, 0.41258, -1);
  gre_u->SetPointError(0, 0.024363*tole, 0.024362*tole,0,0);

  gre_u->Draw("p");

  TGraphAsymmErrors *gre_u_clas = new TGraphAsymmErrors(1);
  gre_u_clas->SetMarkerStyle(8);
  gre_u_clas->SetMarkerSize(2.0);
  gre_u_clas->SetMarkerColor(4);
  gre_u_clas->SetLineColor(4);
  gre_u_clas->SetLineWidth(2);

  gre_u_clas->SetPoint(0, 0.41258, -2);
  gre_u_clas->SetPointError(0,0.0032968*tole,0.0032968*tole,0,0);

  gre_u_clas->Draw("p");


  TGraphAsymmErrors *gre_u_solid = new TGraphAsymmErrors(1);
  gre_u_solid->SetMarkerStyle(8);
  gre_u_solid->SetMarkerSize(2.0);
  gre_u_solid->SetMarkerColor(2);
  gre_u_solid->SetLineColor(2);
  gre_u_solid->SetLineWidth(2);
  
  gre_u_solid->SetPoint(0, 0.41258, -3);
  gre_u_solid->SetPointError(0,0.0023986*tole,0.0023986*tole,0,0);
  
  gre_u_solid->Draw("p");


  c1->cd(2);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.01);
  TH2F* h2 = new TH2F("h2","",100, -0.43, 0.03, 100, -3.5,-0.);
  h2->GetYaxis()->SetLabelSize(0);
  h2->GetYaxis()->SetNdivisions(0);
  h2->GetXaxis()->SetNdivisions(4);
  h2->GetXaxis()->SetLabelSize(0.08);
  h2->GetXaxis()->SetLabelOffset(-0.0);
  h2->Draw();

  TLatex* td = new TLatex();
  td->SetNDC();
  td->SetTextSize(0.1);
  td->DrawLatex(0.4,0.78,"#delta_{T} d^{[0,1]}");
 

  TGraphAsymmErrors *gre_d = new TGraphAsymmErrors(1);
  gre_d->SetMarkerStyle(8);
  gre_d->SetMarkerSize(2.0);
  gre_d->SetMarkerColor(1);
  gre_d->SetLineColor(1);
  gre_d->SetLineWidth(2);

  //#0: Kang 2015
  gre_d->SetPoint(0, -0.22857, -1);
  gre_d->SetPointError(0,0.017224*tole,0.017224*tole,0,0);
  
  gre_d->Draw("p");

  TGraphAsymmErrors *gre_d_clas = new TGraphAsymmErrors(1);
  gre_d_clas->SetMarkerStyle(8);
  gre_d_clas->SetMarkerSize(2.0);
  gre_d_clas->SetMarkerColor(4);
  gre_d_clas->SetLineColor(4);
  gre_d_clas->SetLineWidth(2);

  gre_d_clas->SetPoint(0,-0.22857, -2);
  gre_d_clas->SetPointError(0,0.011196*tole,0.011196*tole,0,0);

  gre_d_clas->Draw("p");

  TGraphAsymmErrors *gre_d_solid = new TGraphAsymmErrors(1);
  gre_d_solid->SetMarkerStyle(8);
  gre_d_solid->SetMarkerSize(2.0);
  gre_d_solid->SetMarkerColor(2);
  gre_d_solid->SetLineColor(2);
  gre_d_solid->SetLineWidth(2);

  gre_d_solid->SetPoint(0,-0.22857,-3);
  gre_d_solid->SetPointError(0,0.0013644*tole,0.0013644*tole,0,0);
 
  gre_d_solid->Draw("p");
 

  const Int_t ny = 3;
  char *label[ny] = {"KPSY2015", "CLAS12", "SoLID"};
 
  c1->cd(3);
  Float_t x, y;
  TLatex *t = new TLatex();
  t->SetNDC();

  int i=0;
  x = 0.03;
  double y_offset=0.23;
  t->SetTextSize(0.06);
  y = 0.65;
  t->SetTextColor(1);
  t->DrawLatex(x,y,label[0]);
  y -=y_offset;
  t->SetTextColor(4);
  t->DrawLatex(x,y,label[1]);
  y -=y_offset;
  t->SetTextColor(2);
  t->DrawLatex(x,y,label[2]);

  // x = 0.0;
  // double y_offset=0.106;
  // t->DrawLatex(x,y,label[6]);
  
  // t->SetTextFont(32);
  // t->SetTextColor(1);
  // t->DrawLatex(0.0,0.82,"Extraction from Experiments:");
  // t->SetTextColor(4);
  // t->DrawLatex(0.0,0.4,"Lattice QCD:");

  c1->Print("Tensor_KPSY_Full.eps");
  c1->Print("Tensor_KPSY_Full.pdf");
}
