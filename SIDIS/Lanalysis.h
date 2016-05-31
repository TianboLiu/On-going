#ifndef _LANALYSIS_H_
#define _LANALYSIS_H_

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

#include "TROOT.h"
#include "TApplication.h"
#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TEventList.h"

#include "wiser.h"
#include "dis.h"

#include "Lgenerator.h"
#include "Lstructure.h"

class Lanalysis{
 protected:
  double ST;
  double lumi;
  double timewindow;
  double Pn;
  double Pp;
 public:
  Lanalysis(double st);
  static int MakeBinInfoTree(const char bininfo[], const char bintree[], const double E0);
  int getsimulationinfo(const char datadir[], double * Nsim, double * genvol, double * days);
  int binanalysis(const char bintreep[], const char bintreem[], const char sidis[], const char datadir[], const int Ntree);
  int ThreetermMatrix(const double * hr, double * M3inv);
  double Threeterm1(const double * coef, const double * phih);
  double Threeterm2(const double * coef, const double * phi);
  double ThreetermDAUT2H(const double * Asym, const double * phi);
  double ThreetermDAUT2S(const double * Asym, const double * phi);
  int ThreetermStatistics(const double * hr, TH1D * h0, double * Estat);
  int EStatisticsUT3(TH1D * h0, double * Estat);
  int EStatisticsUT3(TH2D * hs0, double * Estat);
  int ETargetPolarization(const double * Asym, double * Etarpol);
  int EPolarizationDirection(const double * Asym, double * EtarPD);
  int EResolutionH(const double * hr, const double * Asym, TH2D * dH, double * EresH);
  int EResolutionS(const double * hr, const double * Asym, TH2D * dS, double * EresS);
  int ECoincidence(const double * Asym, TH1D * hsb, double * Ecoin); 
  int CalRMS(double * lab, double * rms, int calls);
  double GetVertexFactor(const double * lab, const double length);
  double RandomCoincidentSigma(const double * lab, int hadron);
};

Lanalysis * lan;
Lanalysis glan(0.7);

Lanalysis::Lanalysis(double st){
  ST = 0.7;
  lumi = 0.1e+10 * pow(0.197327,2);//Polarized Proton
  //0.0334~NH3, 0.01825~4He, 0.3705~total proton, 0.2703~total neutron
  timewindow = 6.0e-9;
  Pn = 0;
  Pp = 1;
}

int Lanalysis::MakeBinInfoTree(const char bininfo[], const char bintree[], const double E0 = 11.0){
  TString infofile = bininfo;
  ifstream infile(infofile);
  if (!infile.is_open()){
    std::cout << "Lanalysis::makebininfotree: No such bin info file!" << std::endl;
    return 1;
  }
  TString filename = bintree;
  TFile * fbin = new TFile(filename, "RECREATE");
  TTree * Fbin = new TTree("bin", "bin");
  Fbin->SetDirectory(fbin);
  double Ebeam = E0;
  double xl, xu;
  double Q2l, Q2u;
  double zl, zu;
  double Ptl, Ptu;
  Fbin->Branch("Ebeam", &Ebeam, "Ebeam/D");
  Fbin->Branch("xl", &xl, "xl/D");
  Fbin->Branch("xu", &xu, "xu/D");
  Fbin->Branch("Q2l", &Q2l, "Q2l/D");
  Fbin->Branch("Q2u", &Q2u, "Q2u/D");
  Fbin->Branch("zl", &zl, "zl/D");
  Fbin->Branch("zu", &zu, "zu/D");
  Fbin->Branch("Ptl", &Ptl, "Ptl/D");
  Fbin->Branch("Ptu", &Ptu, "Ptu/D");
  TString temp;
  infile >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp;
  while (infile >> Q2l >> Q2u >> zl >> zu >> Ptl >> Ptu >> xl >> xu){
    Fbin->Fill();
  }
  infile.close();
  fbin->Write();
  return 0;
}

int Lanalysis::getsimulationinfo(const char datadir[], double * Nsim, double * genvol, double * days){
  TString datadirectory = datadir;
  TString filename = datadirectory+"/out00/run.dat";
  double run[11];
  //lge->setrange(run, filename);
  Nsim[0] = run[0];
  genvol[0] = run[10];
  double E0 = run[9];
  if (E0 == 11.0) days[0] = 55.0;
  if (E0 == 8.8)  days[0] = 27.5;
  return 0;
}

int Lanalysis::binanalysis(const char bintreep[], const char bintreem[], const char sidis[], const char datadir[], const int Ntree = 5){
  TString bininfotreep = bintreep;
  TString bininfotreem = bintreem;
  TString datadirectory = datadir;
  TFile * fbinp = new TFile(bininfotreep,"r");
  if (!fbinp->IsOpen()){
    std::cout << "Lanalysis::binanalysis: Bin info file does not exist!" << std::endl;
    return 1;
  }
  TFile * fbinm = new TFile(bininfotreem,"r");
  if (!fbinm->IsOpen()){
    std::cout << "Lanalysis::binanalysis: Bin info file does not exist!" << std::endl;
    return 1;
  }
  TTree * Fbinp = (TTree *) fbinp->GetObjectChecked("bin", "TTree");
  TTree * Fbinm = (TTree *) fbinm->GetObjectChecked("bin", "TTree");
  int Nbin = Fbinp->GetEntries();
  double Ebeam;
  double xl, xu, xm, ym;
  double Q2l, Q2u, Q2m;
  double zl, zu, zm;
  double Ptl, Ptu, Ptm;
  double Nacc;
  Fbinp->SetBranchAddress("Ebeam", &Ebeam);
  Fbinp->SetBranchAddress("xl", &xl);
  Fbinp->SetBranchAddress("xu", &xu);
  Fbinp->SetBranchAddress("Q2l", &Q2l);
  Fbinp->SetBranchAddress("Q2u", &Q2u);
  Fbinp->SetBranchAddress("zl", &zl);
  Fbinp->SetBranchAddress("zu", &zu);
  Fbinp->SetBranchAddress("Ptl", &Ptl);
  Fbinp->SetBranchAddress("Ptu", &Ptu);
  Fbinm->SetBranchAddress("Ebeam", &Ebeam);
  Fbinm->SetBranchAddress("xl", &xl);
  Fbinm->SetBranchAddress("xu", &xu);
  Fbinm->SetBranchAddress("Q2l", &Q2l);
  Fbinm->SetBranchAddress("Q2u", &Q2u);
  Fbinm->SetBranchAddress("zl", &zl);
  Fbinm->SetBranchAddress("zu", &zu);
  Fbinm->SetBranchAddress("Ptl", &Ptl);
  Fbinm->SetBranchAddress("Ptu", &Ptu);
  double ha, hb;
  double fp, fNH3;
  double A1[2], A1p[2], A1NH3[2];
  double A2[2], A2p[2], A2NH3[2];
  double A3[2], A3p[2], A3NH3[2];  
  TString sidisbinfile = sidis;
  TFile * fs = new TFile(sidisbinfile,"RECREATE");
  TTree * Tp = new TTree("binplus","binplus");
  Tp->SetDirectory(fs);
  TTree * Tm = new TTree("binminus","binminus");
  Tm->SetDirectory(fs);
  Tp->Branch("Ebeam", &Ebeam, "Ebeam/D");
  Tp->Branch("xl", &xl, "xl/D");
  Tp->Branch("xu", &xu, "xu/D");
  Tp->Branch("xm", &xm, "xm/D");
  Tp->Branch("ym", &ym, "ym/D");
  Tp->Branch("Q2l", &Q2l, "Q2l/D");
  Tp->Branch("Q2u", &Q2u, "Q2u/D");
  Tp->Branch("Q2m", &Q2m, "Q2m/D");
  Tp->Branch("zl", &zl, "zl/D");
  Tp->Branch("zu", &zu, "zu/D");
  Tp->Branch("zm", &zm, "zm/D");
  Tp->Branch("Ptl", &Ptl, "Ptl/D");
  Tp->Branch("Ptu", &Ptu, "Ptu/D");
  Tp->Branch("Ptm", &Ptm, "Ptm/D");
  Tp->Branch("ha", &ha, "ha/D");
  Tp->Branch("hb", &hb, "hb/D");
  Tp->Branch("fp", &fp, "fp/D");
  Tp->Branch("fNH3", &fNH3, "fNH3/D");
  Tp->Branch("A1", &A1[0], "A1/D");
  Tp->Branch("A2", &A2[0], "A2/D");
  Tp->Branch("A3", &A3[0], "A3/D");
  Tp->Branch("A1p", &A1p[0], "A1p/D");
  Tp->Branch("A2p", &A2p[0], "A2p/D");
  Tp->Branch("A3p", &A3p[0], "A3p/D");
  Tp->Branch("A1NH3", &A1NH3[0], "A1NH3/D");
  Tp->Branch("A2NH3", &A2NH3[0], "A2NH3/D");
  Tp->Branch("A3NH3", &A3NH3[0], "A3NH3/D");
  Tp->Branch("Nacc", &Nacc, "Nacc/D");
  Tm->Branch("Ebeam", &Ebeam, "Ebeam/D");
  Tm->Branch("xl", &xl, "xl/D");
  Tm->Branch("xu", &xu, "xu/D");
  Tm->Branch("xm", &xm, "xm/D");
  Tm->Branch("ym", &ym, "ym/D");
  Tm->Branch("Q2l", &Q2l, "Q2l/D");
  Tm->Branch("Q2u", &Q2u, "Q2u/D");
  Tm->Branch("Q2m", &Q2m, "Q2m/D");
  Tm->Branch("zl", &zl, "zl/D");
  Tm->Branch("zu", &zu, "zu/D");
  Tm->Branch("zm", &zm, "zm/D");
  Tm->Branch("Ptl", &Ptl, "Ptl/D");
  Tm->Branch("Ptu", &Ptu, "Ptu/D");
  Tm->Branch("Ptm", &Ptm, "Ptm/D");
  Tm->Branch("ha", &ha, "ha/D");
  Tm->Branch("hb", &hb, "hb/D");
  Tm->Branch("fp", &fp, "fp/D");
  Tm->Branch("fNH3", &fNH3, "fNH3/D");
  Tm->Branch("A1", &A1[1], "A1/D");
  Tm->Branch("A2", &A2[1], "A2/D");
  Tm->Branch("A3", &A3[1], "A3/D");
  Tm->Branch("A1p", &A1p[1], "A1p/D");
  Tm->Branch("A2p", &A2p[1], "A2p/D");
  Tm->Branch("A3p", &A3p[1], "A3p/D");
  Tm->Branch("A1NH3", &A1NH3[1], "A1NH3/D");
  Tm->Branch("A2NH3", &A2NH3[1], "A2NH3/D");
  Tm->Branch("A3NH3", &A3NH3[1], "A3NH3/D");
  Tm->Branch("Nacc", &Nacc, "Nacc/D");
  double acc_ele, acc_pion[2];
  double sigma[2], sigmap[2], sigmaNH3[2];
  TString nq, nz;
  TString selectedfile;
  Long64_t Nevent;
  double AZ[4] = {3.705, 2.703, 1.0/3.705, 0};//total
  double AZNH3[4] = {3.333, 2.333, 0.3, 0};//NH3
  double kin[7], lab[7];
  const double lumi = glan.lumi;//Polarized Proton
  double Nsim, genvol, days;
  lan->getsimulationinfo(datadirectory, &Nsim, &genvol, &days);
  std::cout << "bin analysis: pi+" << std::endl;
  //for (int nb = 0; nb < 1; nb++){
  for (int nb = 0; nb < Nbin; nb++){
    Fbinp->GetEntry(nb);
    std::cout << "#" << nb << " in " << Nbin << std::endl;
    if (Q2l == 1.0 && Q2u == 2.0) nq = "1";
    else if (Q2l == 2.0 && Q2u == 3.0) nq = "2";
    else if (Q2l == 3.0 && Q2u == 4.0) nq = "3";
    else if (Q2l == 4.0 && Q2u == 5.0) nq = "4";
    else if (Q2l == 5.0 && Q2u == 6.0) nq = "5";
    else if (Q2l == 6.0 && Q2u == 8.0) nq = "6";
    else continue;
    if (zl == 0.3 && zu == 0.35) nz = "30";
    else if (zl == 0.35 && zu == 0.40) nz = "35";
    else if (zl == 0.40 && zu == 0.45) nz = "40";
    else if (zl == 0.45 && zu == 0.50) nz = "45";
    else if (zl == 0.50 && zu == 0.55) nz = "50";
    else if (zl == 0.55 && zu == 0.60) nz = "55";
    else if (zl == 0.60 && zu == 0.65) nz = "60";
    else if (zl == 0.65 && zu == 0.70) nz = "65";
    else continue;
    selectedfile = "sidis_select"+nq+nz+".root";
    TChain * Tdata = new TChain("T", "T");
    for (int it = 0; it < Ntree; it++){
      Tdata->Add(Form(datadirectory+"/out%.2d/Selected/"+selectedfile, it));
    }
    Nevent = Tdata->GetEntries();
    Tdata->SetBranchAddress("x", &kin[0]);
    Tdata->SetBranchAddress("y", &kin[1]);
    Tdata->SetBranchAddress("z", &kin[2]);
    Tdata->SetBranchAddress("Q2", &kin[3]);
    Tdata->SetBranchAddress("Pt", &kin[4]);
    Tdata->SetBranchAddress("phi_h", &kin[5]);
    Tdata->SetBranchAddress("phi_S", &kin[6]);
    Tdata->SetBranchAddress("Ebeam", &lab[0]);
    Tdata->SetBranchAddress("p_ele", &lab[1]);
    Tdata->SetBranchAddress("theta_ele", &lab[2]);
    Tdata->SetBranchAddress("phi_ele", &lab[3]);
    Tdata->SetBranchAddress("p_pion", &lab[4]);
    Tdata->SetBranchAddress("theta_pion", &lab[5]);
    Tdata->SetBranchAddress("phi_pion", &lab[6]);
    Tdata->SetBranchAddress("acc_ele", &acc_ele);
    Tdata->SetBranchAddress("acc_pion_p", &acc_pion[0]);
    Tdata->SetBranchAddress("acc_pion_m", &acc_pion[1]);
    TH1D * h0p = new TH1D("h0p", "h0p", 180, 0.0, M_PI);
    TH1D * hpp = new TH1D("hpp", "hpp", 180, 0.0, M_PI);
    TH1D * hnp = new TH1D("hnp", "hnp", 180, 0.0, M_PI);
    TH1D * hxp = new TH1D("hxp", "hxp", 180, 0.0, M_PI);
    TH1D * hyp = new TH1D("hyp", "hyp", 180, 0.0, M_PI);
    TH1D * hzp = new TH1D("hzp", "hzp", 180, 0.0, M_PI);
    TH1D * hQ2p = new TH1D("hQ2p", "hQ2p", 180, 0.0, M_PI);
    TH1D * hPtp = new TH1D("hPtp", "hPtp", 180, 0.0, M_PI);
    for (int ie = 0; ie < Nevent; ie++){
      Tdata->GetEntry(ie);
      if (kin[4] < Ptl || kin[4] > Ptu) continue;
      if (kin[0] < xl || kin[0] > xu) continue;
      lst->sigmaUUT(AZ, lab, sigma);
      lst->sigmaUUTp(lab, sigmap);
      lst->sigmaUUT(AZNH3, lab, sigmaNH3);
      h0p->Fill(std::abs(kin[5]), sigma[0]*acc_ele*acc_pion[0]);
      hpp->Fill(std::abs(kin[5]), sigmap[0]*acc_ele*acc_pion[0]);
      hnp->Fill(std::abs(kin[5]), sigmaNH3[0]*acc_ele*acc_pion[0]);
      hxp->Fill(std::abs(kin[5]), kin[0]*sigma[0]*acc_ele*acc_pion[0]);
      hyp->Fill(std::abs(kin[5]), kin[1]*sigma[0]*acc_ele*acc_pion[0]);
      hzp->Fill(std::abs(kin[5]), kin[2]*sigma[0]*acc_ele*acc_pion[0]);
      hQ2p->Fill(std::abs(kin[5]), kin[3]*sigma[0]*acc_ele*acc_pion[0]);
      hPtp->Fill(std::abs(kin[5]), kin[4]*sigma[0]*acc_ele*acc_pion[0]);
    }
    Nacc = h0p->Integral(1,180) * lumi * days * 24.0 * 3600.0 * genvol / Nsim / Ntree;
    if (Nacc > 1e+3){
      fp = hpp->Integral(1,180) / h0p->Integral(1,180);
      fNH3 = hnp->Integral(1,180) / h0p->Integral(1,180);
      xm = hxp->Integral(1,180) / h0p->Integral(1,180);
      ym = hyp->Integral(1,180) / h0p->Integral(1,180);
      zm = hzp->Integral(1,180) / h0p->Integral(1,180);
      Q2m = hQ2p->Integral(1,180) / h0p->Integral(1,180);
      Ptm = hPtp->Integral(1,180) / h0p->Integral(1,180);
      if (h0p->FindFirstBinAbove(1e-9) == 1) ha = 0.0;
      else ha = h0p->GetBinCenter(h0p->FindFirstBinAbove(1e-9));
      if (h0p->FindLastBinAbove(1e-9) == 180) hb = M_PI;
      else hb = h0p->GetBinCenter(h0p->FindLastBinAbove(1e-9));
      kin[0] = xm;
      kin[1] = ym;
      kin[2] = zm;
      kin[3] = Q2m;
      kin[4] = Ptm;
      lst->AsinHmSN(AZ, kin, A1); lst->AsinHmSp(kin, A1p); lst->AsinHmSN(AZNH3, kin, A1NH3);
      lst->AsinHpSN(AZ, kin, A2); lst->AsinHpSp(kin, A2p); lst->AsinHpSN(AZNH3, kin, A2NH3);
      lst->Asin3HmSN(AZ, kin, A3); lst->Asin3HmSp(kin, A3p); lst->Asin3HmSN(AZNH3, kin, A3NH3);
      Tp->Fill();
    }
    h0p->Delete();
    hpp->Delete();
    hnp->Delete();
    hxp->Delete();
    hyp->Delete();
    hzp->Delete();
    hQ2p->Delete();
    hPtp->Delete();
    Tdata->Delete();
  }
  Nbin = Fbinm->GetEntries();
  std::cout << "bin analysis: pi-" << std::endl;
  //for (int nb = 0; nb < 1; nb++){
  for (int nb = 0; nb < Nbin; nb++){
    Fbinm->GetEntry(nb);
    std::cout << "#" << nb << " in " << Nbin << std::endl;
    if (Q2l == 1.0 && Q2u == 2.0) nq = "1";
    else if (Q2l == 2.0 && Q2u == 3.0) nq = "2";
    else if (Q2l == 3.0 && Q2u == 4.0) nq = "3";
    else if (Q2l == 4.0 && Q2u == 5.0) nq = "4";
    else if (Q2l == 5.0 && Q2u == 6.0) nq = "5";
    else if (Q2l == 6.0 && Q2u == 8.0) nq = "6";
    else continue;
    if (zl == 0.3 && zu == 0.35) nz = "30";
    else if (zl == 0.35 && zu == 0.40) nz = "35";
    else if (zl == 0.40 && zu == 0.45) nz = "40";
    else if (zl == 0.45 && zu == 0.50) nz = "45";
    else if (zl == 0.50 && zu == 0.55) nz = "50";
    else if (zl == 0.55 && zu == 0.60) nz = "55";
    else if (zl == 0.60 && zu == 0.65) nz = "60";
    else if (zl == 0.65 && zu == 0.70) nz = "65";
    else continue;
    selectedfile = "sidis_select"+nq+nz+".root";
    TChain * Tdata = new TChain("T", "T");
    for (int it = 0; it < Ntree; it++){
      Tdata->Add(Form(datadirectory+"/out%.2d/Selected/"+selectedfile, it));
    }
    Nevent = Tdata->GetEntries();
    Tdata->SetBranchAddress("x", &kin[0]);
    Tdata->SetBranchAddress("y", &kin[1]);
    Tdata->SetBranchAddress("z", &kin[2]);
    Tdata->SetBranchAddress("Q2", &kin[3]);
    Tdata->SetBranchAddress("Pt", &kin[4]);
    Tdata->SetBranchAddress("phi_h", &kin[5]);
    Tdata->SetBranchAddress("phi_S", &kin[6]);
    Tdata->SetBranchAddress("Ebeam", &lab[0]);
    Tdata->SetBranchAddress("p_ele", &lab[1]);
    Tdata->SetBranchAddress("theta_ele", &lab[2]);
    Tdata->SetBranchAddress("phi_ele", &lab[3]);
    Tdata->SetBranchAddress("p_pion", &lab[4]);
    Tdata->SetBranchAddress("theta_pion", &lab[5]);
    Tdata->SetBranchAddress("phi_pion", &lab[6]);
    Tdata->SetBranchAddress("acc_ele", &acc_ele);
    Tdata->SetBranchAddress("acc_pion_p", &acc_pion[0]);
    Tdata->SetBranchAddress("acc_pion_m", &acc_pion[1]);
    TH1D * h0m = new TH1D("h0m", "h0m", 180, 0.0, M_PI);
    TH1D * hpm = new TH1D("hpm", "hpm", 180, 0.0, M_PI);
    TH1D * hnm = new TH1D("hnm", "hnm", 180, 0.0, M_PI);
    TH1D * hxm = new TH1D("hxm", "hxm", 180, 0.0, M_PI);
    TH1D * hym = new TH1D("hym", "hym", 180, 0.0, M_PI);
    TH1D * hzm = new TH1D("hzm", "hzm", 180, 0.0, M_PI);
    TH1D * hQ2m = new TH1D("hQ2m", "hQ2m", 180, 0.0, M_PI);
    TH1D * hPtm = new TH1D("hPtm", "hPtm", 180, 0.0, M_PI);
    for (int ie = 0; ie < Nevent; ie++){
      Tdata->GetEntry(ie);
      if (kin[4] < Ptl || kin[4] > Ptu) continue;
      if (kin[0] < xl || kin[0] > xu) continue;
      lst->sigmaUUT(AZ, lab, sigma);
      lst->sigmaUUTp(lab, sigmap);
      lst->sigmaUUT(AZNH3, lab, sigmaNH3);
      h0m->Fill(std::abs(kin[5]), sigma[1]*acc_ele*acc_pion[1]);
      hpm->Fill(std::abs(kin[5]), sigmap[1]*acc_ele*acc_pion[1]);
      hnm->Fill(std::abs(kin[5]), sigmaNH3[1]*acc_ele*acc_pion[1]);
      hxm->Fill(std::abs(kin[5]), kin[0]*sigma[1]*acc_ele*acc_pion[1]);
      hym->Fill(std::abs(kin[5]), kin[1]*sigma[1]*acc_ele*acc_pion[1]);
      hzm->Fill(std::abs(kin[5]), kin[2]*sigma[1]*acc_ele*acc_pion[1]);
      hQ2m->Fill(std::abs(kin[5]), kin[3]*sigma[1]*acc_ele*acc_pion[1]);
      hPtm->Fill(std::abs(kin[5]), kin[4]*sigma[1]*acc_ele*acc_pion[1]);
    }
    Nacc = h0m->Integral(1,180) * lumi * days * 24.0 * 3600.0 * genvol / Nsim / Ntree;
    if (Nacc > 1e+3){
      fp = hpm->Integral(1,180) / h0m->Integral(1,180);
      fNH3 = hnm->Integral(1,180) / h0m->Integral(1,180);
      xm = hxm->Integral(1,180) / h0m->Integral(1,180);
      ym = hym->Integral(1,180) / h0m->Integral(1,180);
      zm = hzm->Integral(1,180) / h0m->Integral(1,180);
      Q2m = hQ2m->Integral(1,180) / h0m->Integral(1,180);
      Ptm = hPtm->Integral(1,180) / h0m->Integral(1,180);
      if (h0m->FindFirstBinAbove(1e-9) == 1) ha = 0.0;
      else ha = h0m->GetBinCenter(h0m->FindFirstBinAbove(1e-9));
      if (h0m->FindLastBinAbove(1e-9) == 180) hb = M_PI;
      else hb = h0m->GetBinCenter(h0m->FindLastBinAbove(1e-9));
      kin[0] = xm;
      kin[1] = ym;
      kin[2] = zm;
      kin[3] = Q2m;
      kin[4] = Ptm;
      lst->AsinHmSN(AZ, kin, A1); lst->AsinHmSp(kin, A1p); lst->AsinHmSN(AZNH3, kin, A1NH3);
      lst->AsinHpSN(AZ, kin, A2); lst->AsinHpSp(kin, A2p); lst->AsinHpSN(AZNH3, kin, A2NH3);
      lst->Asin3HmSN(AZ, kin, A3); lst->Asin3HmSp(kin, A3p); lst->Asin3HmSN(AZNH3, kin, A3NH3);
      Tm->Fill();
    }
    h0m->Delete();
    hpm->Delete();
    hnm->Delete();
    hxm->Delete();
    hym->Delete();
    hzm->Delete();
    hQ2m->Delete();
    hPtm->Delete();
    Tdata->Delete();
  }
  fs->Write();
  return 0;
}

int Lanalysis::ThreetermMatrix(const double * hr, double * M3inv){
  double a = hr[0];
  double b = hr[1];
  const double pi = M_PI;
  TMatrixD M3(3,3);
  M3(0,0) = b - a; 
  M3(0,1) = 0.5 * (sin(2.0*a) - sin(2.0*b));
  M3(0,2) = 0.5 * (sin(2.0*b) - sin(2.0*a));
  M3(1,0) = 0.5 * (sin(2.0*a) - sin(2.0*b));
  M3(1,1) = b - a;
  M3(1,2) = 0.25 * (sin(4.0*a) - sin(4.0*b));
  M3(2,0) = 0.5 * (sin(2.0*b) - sin(2.0*a));
  M3(2,1) = 0.25 * (sin(4.0*a) - sin(4.0*b));
  M3(2,2) = b - a;
  M3.Invert();
  M3inv[0] = M3(0,0)/pi;
  M3inv[1] = M3(0,1)/pi;
  M3inv[2] = M3(0,2)/pi;
  M3inv[3] = M3(1,0)/pi;
  M3inv[4] = M3(1,1)/pi;
  M3inv[5] = M3(1,2)/pi;
  M3inv[6] = M3(2,0)/pi;
  M3inv[7] = M3(2,1)/pi;
  M3inv[8] = M3(2,2)/pi;
  return 0;
}

double Lanalysis::Threeterm1(const double * coef, const double * phi){
  double phih = phi[0];
  double a1 = coef[0];
  double a2 = coef[1];
  double a3 = coef[2];
  return 0.5 * M_PI * (a1*a1 + a2*a2 + a3*a3 - 2.0*a1*a2*cos(2.0*phih) - 2.0*a2*a3*cos(4.0*phih) + 2.0*a1*a3*cos(2.0*phih));
}

double Lanalysis::Threeterm2(const double * coef, const double * phi){
  double phih = phi[0];
  double phiS = phi[1];
  double a1 = coef[0];
  double a2 = coef[1];
  double a3 = coef[2];
  return pow(a1*sin(phih-phiS) + a2*sin(phih+phiS) + a3*sin(3.0*phih-phiS), 2);
}

int Lanalysis::ThreetermStatistics(const double * hr, TH1D * h0, double * Estat){
  //const double eps = 1e-5;
  double M3inv[9];
  lan->ThreetermMatrix(hr, M3inv);
  int Nb = h0->GetNbinsX();
  double average = h0->Integral(1, Nb) / Nb;
  double Nacc = average * Nb;
  double phih;
  double Es[3] = {0.0, 0.0, 0.0};
  double gh;
  for (int i = 1; i <= Nb; i++){
    phih = h0->GetBinCenter(i);
    gh = h0->GetBinContent(i) / average;
    Es[0] = Es[0] + lan->Threeterm1(&M3inv[0], &phih) * gh * 2.0 * M_PI * M_PI / Nacc;
    Es[1] = Es[1] + lan->Threeterm1(&M3inv[3], &phih) * gh * 2.0 * M_PI * M_PI / Nacc;
    Es[2] = Es[2] + lan->Threeterm1(&M3inv[6], &phih) * gh * 2.0 * M_PI * M_PI / Nacc;
  }
  double width = h0->GetBinWidth(1);
  Estat[0] = sqrt(2.0 * Es[0] * width) / glan.ST;
  Estat[1] = sqrt(2.0 * Es[1] * width) / glan.ST;
  Estat[2] = sqrt(2.0 * Es[2] * width) / glan.ST;
  return 0;
}

int Lanalysis::EStatisticsUT3(TH1D * h0, double * Estat){
  int Nb = h0->GetNbinsX();
  double average = h0->Integral(1, Nb) / Nb;
  double Nacc = average * Nb;
  double gh;
  double c0 = 0;
  double c2 = 0;
  double c4 = 0;
  double phih;
  for (int i = 1; i <= Nb; i++){
    phih = h0->GetBinCenter(i);
    gh = h0->GetBinContent(i) / average;
    c0 = c0 + M_PI * gh;
    c2 = c2 + M_PI * cos(2.0 * phih) * gh;
    c4 = c4 + M_PI * cos(4.0 * phih) * gh;
  }
  double width = h0->GetBinWidth(1);
  c0 = c0 * width;
  c2 = c2 * width;
  c4 = c4 * width;
  TMatrixD M3(3,3);
  M3(0, 0) = c0; M3(0, 1) = -c2; M3(0, 2) = c2;
  M3(1, 0) = -c2; M3(1, 1) = c0; M3(1, 2) = -c4;
  M3(2, 0) = c2; M3(2, 1) = -c4; M3(2, 2) = c0;
  M3.Invert();
  Estat[0] = sqrt(2.0*M_PI*M_PI/Nacc * (pow(M3(0,0),2) + pow(M3(0,1),2) + pow(M3(0,2),2)) * M_PI * M_PI) / glan.ST;
  Estat[1] = sqrt(2.0*M_PI*M_PI/Nacc * (pow(M3(1,0),2) + pow(M3(1,1),2) + pow(M3(1,2),2)) * M_PI * M_PI) / glan.ST;
  Estat[2] = sqrt(2.0*M_PI*M_PI/Nacc * (pow(M3(2,0),2) + pow(M3(2,1),2) + pow(M3(2,2),2)) * M_PI * M_PI) / glan.ST;
  return 0;
};

int Lanalysis::EStatisticsUT3(TH2D * hs0, double * Estat){
  int NX = hs0->GetNbinsX();
  int NY = hs0->GetNbinsY();
  double average = hs0->Integral(1, NX, 1, NY) / NX / NY;
  double Nacc = average * NX * NY;
  double ghs;
  double cc[6] = {0, 0, 0, 0, 0, 0};
  double phih, phiS;
  for (int i = 1; i < NX; i++){
    for (int j = 1; j < NY; j++){
      phih = hs0->GetXaxis()->GetBinCenter(i);
      phiS = hs0->GetYaxis()->GetBinCenter(j);
      ghs = hs0->GetBinContent(i, j) / average;
      cc[0] = cc[0] + sin(phih - phiS) * sin(phih - phiS) * ghs;
      cc[1] = cc[1] + sin(phih - phiS) * sin(phih + phiS) * ghs;
      cc[2] = cc[2] + sin(phih - phiS) * sin(3.0*phih - phiS) * ghs;
      cc[3] = cc[3] + sin(phih + phiS) * sin(phih + phiS) * ghs;
      cc[4] = cc[4] + sin(phih + phiS) * sin(3.0*phih - phiS) * ghs;
      cc[5] = cc[5] + sin(3.0*phih - phiS) * sin(3.0*phih - phiS) * ghs;
    }
  }
  double volume = 2.0 * M_PI * M_PI / NX / NY;
  for (int n = 0; n < 6; n++){
    cc[n] = cc[n] * volume;
  }
  TMatrixD M3(3,3);
  M3(0, 0) = cc[0]; M3(0, 1) = cc[1]; M3(0, 2) = cc[2];
  M3(1, 0) = cc[1]; M3(1, 1) = cc[3]; M3(1, 2) = cc[4];
  M3(2, 0) = cc[2]; M3(2, 1) = cc[4]; M3(2, 2) = cc[5];
  M3.Invert();
  Estat[0] = sqrt(2.0*M_PI*M_PI/Nacc * (pow(M3(0,0),2) + pow(M3(0,1),2) + pow(M3(0,2),2)) * M_PI * M_PI) / glan.ST;
  Estat[1] = sqrt(2.0*M_PI*M_PI/Nacc * (pow(M3(1,0),2) + pow(M3(1,1),2) + pow(M3(1,2),2)) * M_PI * M_PI) / glan.ST;
  Estat[2] = sqrt(2.0*M_PI*M_PI/Nacc * (pow(M3(2,0),2) + pow(M3(2,1),2) + pow(M3(2,2),2)) * M_PI * M_PI) / glan.ST;
  return 0;
}

int Lanalysis::ETargetPolarization(const double * Asym, double * Etarpol){
  for (int i = 0; i < 3; i++){
    Etarpol[i] = 0.03 * std::abs(Asym[i]);
  }
  return 0;
}

int Lanalysis::EPolarizationDirection(const double * Asym, double * EtarPD){
  const double dth = 0.2 * M_PI / 180.0;
  for (int i = 0; i < 3; i++){
    EtarPD[i] = std::abs(Asym[i]) * dth * dth / 2.0;
  }
  return 0;
} 

TFile * _file_e = new TFile("/var/phy/project/mepg/tl190/SIDIS/SIDIS_electron_resolution_2d.root","r");
TH2D * _theta_e = (TH2D *) _file_e->GetObjectChecked("theta_resolution","TH2D");
TH2D * _phi_e = (TH2D *) _file_e->GetObjectChecked("phi_resolution","TH2D");
TH2D * _p_e = (TH2D *) _file_e->GetObjectChecked("p_resolution","TH2D");
TH2D * _z_e = (TH2D *) _file_e->GetObjectChecked("vertexz_resolution","TH2D");
TFile * _file_pi = new TFile("/var/phy/project/mepg/tl190/SIDIS/SIDIS_pim_resolution_2d.root","r");
TH2D * _theta_pi = (TH2D *) _file_e->GetObjectChecked("theta_resolution","TH2D");
TH2D * _phi_pi = (TH2D *) _file_e->GetObjectChecked("phi_resolution","TH2D");
TH2D * _p_pi = (TH2D *) _file_e->GetObjectChecked("p_resolution","TH2D");
TH2D * _z_pi = (TH2D *) _file_e->GetObjectChecked("vertexz_resolution","TH2D");

int Lanalysis::CalRMS(double * lab, double * rms, int calls = 100){
  //lab: Ebeam, p_ele, theta_ele, phi_ele, p_pion, theta_pion, phi_pion
  const double radtodeg = 180.0 / M_PI;
  int bin_p_e = _theta_e->GetXaxis()->FindBin(lab[1]);
  int bin_theta_e = _theta_e->GetYaxis()->FindBin(lab[2]);
  int bin_p_pi = _theta_pi->GetXaxis()->FindBin(lab[4]);
  int bin_theta_pi = _theta_pi->GetYaxis()->FindBin(lab[5]);
  double dw[7];
  dw[0] = 0;
  dw[1] = _p_e->GetBinContent(bin_p_e, bin_theta_e) / 100.0;
  dw[2] = _theta_e->GetBinContent(bin_p_e, bin_theta_e) / 1000.0 * radtodeg;
  dw[3] = _phi_e->GetBinContent(bin_p_e, bin_theta_e) / 1000.0 * radtodeg;
  dw[4] = _p_pi->GetBinContent(bin_p_pi, bin_theta_pi) / 100.0;
  dw[5] = _theta_pi->GetBinContent(bin_p_pi, bin_theta_pi) / 1000.0 * radtodeg;
  dw[6] = _phi_pi->GetBinContent(bin_p_pi, bin_theta_pi) / 1000.0 * radtodeg;
  double phys[9];//x, y, z, Q2, Pt, phih, phiS, W, Wp
  lst->CalcVariables(lab, phys);
  double newlab[7], newphys[9], d2[7];
  double accum[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (int n = 0; n < calls; n++){
    for (int i = 0; i < 7; i++){
      newlab[i] = lab[i] + gRandom->Gaus(0.0, dw[i]);
    }
    lst->CalcVariables(newlab, newphys);
    if (newphys[5] - phys[5] > M_PI) newphys[5] = newphys[5] - 2.0*M_PI;
    if (newphys[5] - phys[5] < -M_PI) newphys[5] = newphys[5] + 2.0*M_PI;
    if (newphys[6] - phys[6] > M_PI) newphys[6] = newphys[6] - 2.0*M_PI;
    if (newphys[6] - phys[6] < -M_PI) newphys[6] = newphys[6] + 2.0*M_PI;
    for (int i = 0; i < 7; i++){
      d2[i] = pow(newphys[i] - phys[i], 2);
      accum[i] = accum[i] + d2[i];
    }
  }
  for (int i = 0; i < 7; i++){
    rms[i] = sqrt(accum[i]/calls);
  }
  rms[6] = sqrt(rms[6]*rms[6] + pow(0.2*M_PI/180.0,2));
  return 0;
}

double  Lanalysis::GetVertexFactor(const double * lab, const double length = 40.0){
  double ze = _z_e->GetBinContent(_z_e->GetXaxis()->FindBin(lab[1]),_z_e->GetYaxis()->FindBin(lab[2]));
  double zpi = _z_pi->GetBinContent(_z_pi->GetXaxis()->FindBin(lab[4]),_z_pi->GetYaxis()->FindBin(lab[5]));
  return length / (3.0 * (ze + zpi));
}
  

double Lanalysis::ThreetermDAUT2H(const double * Asym, const double * phi){
  double A1 = Asym[0];
  double A2 = Asym[1];
  double A3 = Asym[2];
  double phih = phi[0];
  double phiS = phi[1];
  return pow(A1*cos(phih-phiS) + A2*cos(phih+phiS) + 3.0*A3*cos(3.0*phih-phiS), 2);
}

double Lanalysis::ThreetermDAUT2S(const double * Asym, const double * phi){
  double A1 = Asym[0];
  double A2 = Asym[1];
  double A3 = Asym[2];
  double phih = phi[0];
  double phiS = phi[1];
  return pow(A1*cos(phih-phiS) - A2*cos(phih+phiS) + A3*cos(3.0*phih-phiS), 2);
}

int Lanalysis::EResolutionH(const double * hr, const double * Asym, TH2D * dH, double * EresH){
  const double eps = 1.0e-5;
  double M3inv[9];
  lan->ThreetermMatrix(hr, M3inv);
  double phi[2];
  int bina1 = dH->GetXaxis()->FindBin(-hr[1]+eps);
  int binb1 = dH->GetXaxis()->FindBin(-hr[0]-eps);
  int bina2 = dH->GetXaxis()->FindBin(hr[0]+eps);
  int binb2 = dH->GetXaxis()->FindBin(hr[1]-eps);
  int Ns = dH->GetYaxis()->FindBin(M_PI-eps);
  double wphih = dH->GetXaxis()->GetBinWidth(1);
  double wphiS = dH->GetYaxis()->GetBinWidth(1);
  double delta;
  double separate2[3];
  double dAUT2;
  double Es[3] = {0.0, 0.0, 0.0};
  for (int i = bina1; i <= binb1; i++){
    for (int j = 1; j <= Ns; j++){
      delta = dH->GetBinContent(i,j);
      phi[0] = dH->GetXaxis()->GetBinCenter(i);
      phi[1] = dH->GetYaxis()->GetBinCenter(j);
      dAUT2 = lan->ThreetermDAUT2H(Asym, phi);
      separate2[0] = lan->Threeterm2(&M3inv[0], phi);
      separate2[1] = lan->Threeterm2(&M3inv[3], phi);
      separate2[2] = lan->Threeterm2(&M3inv[6], phi);
      Es[0] = Es[0] + delta*delta*dAUT2*separate2[0];
      Es[1] = Es[1] + delta*delta*dAUT2*separate2[1];
      Es[2] = Es[2] + delta*delta*dAUT2*separate2[2];
    }
  }
  for (int i = bina2; i <= binb2; i++){
    for (int j = 1; j <= Ns; j++){
      delta = dH->GetBinContent(i,j);
      phi[0] = dH->GetXaxis()->GetBinCenter(i);
      phi[1] = dH->GetYaxis()->GetBinCenter(j);
      dAUT2 = lan->ThreetermDAUT2H(Asym, phi);
      separate2[0] = lan->Threeterm2(&M3inv[0], phi);
      separate2[1] = lan->Threeterm2(&M3inv[3], phi);
      separate2[2] = lan->Threeterm2(&M3inv[6], phi);
      Es[0] = Es[0] + delta*delta*dAUT2*separate2[0];
      Es[1] = Es[1] + delta*delta*dAUT2*separate2[1];
      Es[2] = Es[2] + delta*delta*dAUT2*separate2[2];
    }
  }
  EresH[0] = sqrt(Es[0] * wphih * wphiS);
  EresH[1] = sqrt(Es[1] * wphih * wphiS);
  EresH[2] = sqrt(Es[2] * wphih * wphiS);
  return 0;
}

int Lanalysis::EResolutionS(const double * hr, const double * Asym, TH2D * dS, double * EresS){
  const double eps = 1.0e-5;
  double M3inv[9];
  lan->ThreetermMatrix(hr, M3inv);
  double phi[2];
  int bina1 = dS->GetXaxis()->FindBin(-hr[1]+eps);
  int binb1 = dS->GetXaxis()->FindBin(-hr[0]-eps);
  int bina2 = dS->GetXaxis()->FindBin(hr[0]+eps);
  int binb2 = dS->GetXaxis()->FindBin(hr[1]-eps);
  int Ns = dS->GetYaxis()->FindBin(M_PI-eps);
  double wphih = dS->GetXaxis()->GetBinWidth(1);
  double wphiS = dS->GetYaxis()->GetBinWidth(1);
  double delta;
  double separate2[3];
  double dAUT2;
  double Es[3] = {0.0, 0.0, 0.0};
  for (int i = bina1; i <= binb1; i++){
    for (int j = 1; j <= Ns; j++){
      delta = dS->GetBinContent(i,j);
      phi[0] = dS->GetXaxis()->GetBinCenter(i);
      phi[1] = dS->GetYaxis()->GetBinCenter(j);
      dAUT2 = lan->ThreetermDAUT2S(Asym, phi);
      separate2[0] = lan->Threeterm2(&M3inv[0], phi);
      separate2[1] = lan->Threeterm2(&M3inv[3], phi);
      separate2[2] = lan->Threeterm2(&M3inv[6], phi);
      Es[0] = Es[0] + delta*delta*dAUT2*separate2[0];
      Es[1] = Es[1] + delta*delta*dAUT2*separate2[1];
      Es[2] = Es[2] + delta*delta*dAUT2*separate2[2];
    }
  }
  for (int i = bina2; i <= binb2; i++){
    for (int j = 1; j <= Ns; j++){
      delta = dS->GetBinContent(i,j);
      phi[0] = dS->GetXaxis()->GetBinCenter(i);
      phi[1] = dS->GetYaxis()->GetBinCenter(j);
      dAUT2 = lan->ThreetermDAUT2S(Asym, phi);
      separate2[0] = lan->Threeterm2(&M3inv[0], phi);
      separate2[1] = lan->Threeterm2(&M3inv[3], phi);
      separate2[2] = lan->Threeterm2(&M3inv[6], phi);
      Es[0] = Es[0] + delta*delta*dAUT2*separate2[0];
      Es[1] = Es[1] + delta*delta*dAUT2*separate2[1];
      Es[2] = Es[2] + delta*delta*dAUT2*separate2[2];
    }
  }
  EresS[0] = sqrt(Es[0] * wphih * wphiS);
  EresS[1] = sqrt(Es[1] * wphih * wphiS);
  EresS[2] = sqrt(Es[2] * wphih * wphiS);
  return 0;
}

double Lanalysis::RandomCoincidentSigma(const double * lab, int hadron){
  const double degtorad = M_PI / 180.0;
  double radlen;
  if (lab[0] == 11.0) radlen = 0.04896;
  if (lab[0] == 8.8) radlen = 0.04793;
  double sigma_pi = wiser_sigma(lab[0], lab[4], lab[5] * degtorad, radlen, hadron) / 3.89379e+5;
  double sigma_e = calcombine(lab[0], lab[1], lab[2] * degtorad, 3.705, 2.705);
  //double vertexfactor = lan->GetVertexFactor(lab);
  return sigma_e * sigma_pi * glan.lumi * glan.timewindow / 1.0;
}

int Lanalysis::ECoincidence(const double * Asym, TH1D * hsb, double * Ecoin){
  double sbratio = hsb->GetMean();
  double relative = 0.2 / sbratio;
  for (int i = 0; i < 3; i++){
    Ecoin[i] = std::abs(Asym[i] * relative);
  }
  return 0;
}



#endif
