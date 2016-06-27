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
#include "TF1.h"
#include "TF2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TEventList.h"

#include "wiser.h"

#include "Lgenerator.h"
#include "Lstructure.h"

class Lanalysis{
 protected:
  TString _datadir;
  TString _binfile1;
  TString _binfile2;
  double _eff;
  double _ST;
  double _lumi;
  double _days;
  double _simdensity;
  int _Ntree;
 private:
  int _res_ctrl;
  TFile * _file_e;
  TH2D * _theta_e;
  TH2D * _phi_e;
  TH2D * _p_e;
  TH2D * _z_e;
  TFile * _file_pi;
  TH2D * _theta_pi;
  TH2D * _phi_pi;
  TH2D * _p_pi;
  TH2D * _z_pi;
 public:
  Lanalysis(TString datadir);
  Lanalysis(TString datadir, TString binfile1, TString binfile2);
  static int MakeBinInfoTree(const char bininfo[], const char bintree[], const double E0);
  static int FitCos1(TH1D * h0, const double * range, double * cc);
  static int FitCosSin1(TH1D * h0, const double * range, double * cc);
  int SetSimInfo(double lumi, double days, double ST, int Ntree);
  int BinAnalysisNeutron(const char savefile[]);
  int BinAnalysisProton(const char savefile[]);
  int BinResolutionNeutron(const char bintree[], const char savefile[]);
  int BinResolutionProton(const char bintree[], const char savefile[]);
  int BinAcceptanceNeutron(const char bintree[], const char savefile[]);
  int BinAcceptanceProton(const char bintree[], const char savefile[]);
  int ECoincidenceNeutron(const char bintree[], const char rmstree[], const char savefile[]);
  int ECoincidenceProton(const char bintree[], const char rmstree[], const char savefile[]);
  int EResolutionNeutron(const char bintree[], const char rmstree[], const char acctree[], const char savefile[]);
  int EResolutionProton(const char bintree[], const char rmstree[], const char acctree[], const char savefile[]);
  int ENuclearPDF(const char bintree[], const char savefile[]);
  int ENuclearNeutron(const char bintree[], const char savefile[]);
  int EDilutionNeutron(const char bintree[], const char savefile[]);
  int EDilutionProton(const char bintree[], const char savefile[]);
  int ERadiativeNeutron(const char bintree[], const char savefile[]);
  int ERadiativeProton(const char bintree[], const char savefile[]);
  int EExclusiveNeutron(const char bintree[], const char savefile[]);
  int EExclusiveProton(const char bintree[], const char savefile[]);
  int EDiffractiveNeutron(const char bintree[], const char savefile[]);
  int EDiffractiveProton(const char bintree[], const char savefile[]);
  int ETotalNeutron(const char dir[], const char savefile[]);
  int ETotalProton(const char dir[], const char savefile[]);
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
  int GetResolutionFiles();
  int EResolutionH(const double * hr, const double * Asym, TH2D * dH, double * EresH);
  int EResolutionS(const double * hr, const double * Asym, TH2D * dS, double * EresS);
  int CalRMS(double * lab, double * rms, int calls);
  double GetVertexFactor(const double * lab, const double length);
  int RandomCoincidenceSigmaN(const double * AZ, const double * lab, double * xs);
};

Lanalysis::Lanalysis(TString datadir){
  _datadir = datadir;
  _res_ctrl = 0;
}

Lanalysis::Lanalysis(TString datadir, TString binfile1, TString binfile2){
  _datadir = datadir;
  _binfile1 = binfile1;
  _binfile2 = binfile2;
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

int Lanalysis::FitCos1(TH1D * h0, const double * range, double * cc){
  TF1 * fc = new TF1("fc", "[0]*(1.0+[1]*cos(x)+[2]*cos(2.0*x)+[3]*cos(3.0*x)+[4]*cos(4.0*x)+[5]*cos(5.0*x))", 0.0, M_PI);
  fc->SetParameters(h0->Integral(1, -1), 0.0, 0.0, 0.0, 0.0, 0.0);
  h0->Fit("fc", "WIQO", "", range[0], range[1]);
  TF1 * rs = h0->GetFunction("fc");
  cc[0] = rs->GetParameter(0);
  cc[1] = rs->GetParameter(1);
  cc[2] = rs->GetParameter(2);
  cc[3] = rs->GetParameter(3);
  cc[4] = rs->GetParameter(4);
  cc[5] = rs->GetParameter(5);
  return 0;
}

int Lanalysis::FitCosSin1(TH1D * h0, const double * range, double * cc){
  TF1 * fc = new TF1("fc", "[0]*(1.0+[1]*cos(x)+[2]*cos(2.0*x)+[3]*cos(3.0*x)+[4]*cos(4.0*x)+[5]*cos(5.0*x)+[6]*sin(x)+[7]*sin(2.0*x)+[8]*sin(3.0*x)+[9]*sin(4.0*x)+[10]*sin(5.0*x))", -M_PI, M_PI);
  fc->SetParameters(h0->Integral(1, -1), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  h0->Fit("fc", "WIQO", "", range[0], range[1]);
  TF1 * rs = h0->GetFunction("fc");
  cc[0] = rs->GetParameter(0);
  cc[1] = rs->GetParameter(1);
  cc[2] = rs->GetParameter(2);
  cc[3] = rs->GetParameter(3);
  cc[4] = rs->GetParameter(4);
  cc[5] = rs->GetParameter(5);
  cc[6] = rs->GetParameter(6);
  cc[7] = rs->GetParameter(7);
  cc[8] = rs->GetParameter(8);
  cc[9] = rs->GetParameter(9);
  cc[10] = rs->GetParameter(10);
  return 0;
}
  
int Lanalysis::SetSimInfo(double lumi, double days, double ST, int Ntree){
  _eff = 0.85;
  _lumi = lumi;
  _days = days;
  _ST = ST;
  _Ntree = Ntree;
  TString filename = _datadir+"/out00/run.dat";
  ifstream infile(filename);
  double run[11];
  if (!infile.is_open()){
    std::cout << "Run file does not exist!" << std::endl;
    return 1;
  }
  TString temp;
  int i = 0;
  while (infile >> temp >> run[i]){i++;}
  infile.close();
  double genvol = (run[2] - run[1]) * (cos(run[3] * M_PI / 180.0) - cos(run[4] * M_PI / 180.0)) * (run[6] - run[5]) * (cos(run[7] * M_PI / 180.0) - cos(run[8] * M_PI / 180.0)) * (2.0 * M_PI) * (2.0 * M_PI);
  _simdensity = run[0] / genvol * Ntree;
  return 0;
}  

int Lanalysis::BinAnalysisNeutron(const char savefile[]){
  TFile * fbinp = new TFile(_binfile1,"r");
  if (!fbinp->IsOpen()){
    std::cout << "Lanalysis::binanalysis: Bin info file does not exist!" << std::endl;
    return 1;
  }
  TFile * fbinm = new TFile(_binfile2,"r");
  if (!fbinm->IsOpen()){
    std::cout << "Lanalysis::binanalysis: Bin info file does not exist!" << std::endl;
    return 1;
  }
  TTree * Fbinp = (TTree *) fbinp->GetObjectChecked("bin", "TTree");
  TTree * Fbinm = (TTree *) fbinm->GetObjectChecked("bin", "TTree");
  int Nbinp = Fbinp->GetEntries();
  int Nbinm = Fbinm->GetEntries();
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
  double BinNumber, Nsim;
  double ha, hb, fn;
  double A1[2], A1p[2], A1n[2];
  double A2[2], A2p[2], A2n[2];
  double A3[2], A3p[2], A3n[2];
  double Estat[3];
  TFile * fs = new TFile(savefile,"RECREATE");
  TTree * Tp = new TTree("binplus","binplus");
  Tp->SetDirectory(fs);
  TTree * Tm = new TTree("binminus","binminus");
  Tm->SetDirectory(fs);
  Tp->Branch("BinNumber", &BinNumber, "BinNumber/D");
  Tp->Branch("Nsim", &Nsim, "Nsim/D");
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
  Tp->Branch("fn", &fn, "fn/D");
  Tp->Branch("A1", &A1[0], "A1/D");
  Tp->Branch("A2", &A2[0], "A2/D");
  Tp->Branch("A3", &A3[0], "A3/D");
  Tp->Branch("A1p", &A1p[0], "A1p/D");
  Tp->Branch("A2p", &A2p[0], "A2p/D");
  Tp->Branch("A3p", &A3p[0], "A3p/D");
  Tp->Branch("A1n", &A1n[0], "A1n/D");
  Tp->Branch("A2n", &A2n[0], "A2n/D");
  Tp->Branch("A3n", &A3n[0], "A3n/D");
  Tp->Branch("Nacc", &Nacc, "Nacc/D");
  Tp->Branch("E1stat", &Estat[0], "E1stat/D");
  Tp->Branch("E2stat", &Estat[1], "E2stat/D");
  Tp->Branch("E3stat", &Estat[2], "E3stat/D");
  Tm->Branch("BinNumber", &BinNumber, "BinNumber/D");
  Tm->Branch("Nsim", &Nsim, "Nsim/D");
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
  Tm->Branch("fn", &fn, "fn/D");
  Tm->Branch("A1", &A1[1], "A1/D");
  Tm->Branch("A2", &A2[1], "A2/D");
  Tm->Branch("A3", &A3[1], "A3/D");
  Tm->Branch("A1p", &A1p[1], "A1p/D");
  Tm->Branch("A2p", &A2p[1], "A2p/D");
  Tm->Branch("A3p", &A3p[1], "A3p/D");
  Tm->Branch("A1n", &A1n[1], "A1n/D");
  Tm->Branch("A2n", &A2n[1], "A2n/D");
  Tm->Branch("A3n", &A3n[1], "A3n/D");
  Tm->Branch("Nacc", &Nacc, "Nacc/D");
  Tm->Branch("E1stat", &Estat[0], "E1stat/D");
  Tm->Branch("E2stat", &Estat[1], "E2stat/D");
  Tm->Branch("E3stat", &Estat[2], "E3stat/D");
  double acc_ele, acc_pion[2];
  double sigmaN[2], sigman[2];
  double AZ[4] = {2, 1, -0.028, 0.86};//total
  TString nq, nz;
  TString selectedfile;
  Long64_t Nevent;
  double kin[7], lab[7];
  std::cout << "bin analysis: pi+" << std::endl;
  //for (int nb = 684; nb < 686; nb++){
  for (int nb = 0; nb < Nbinp; nb++){
    Fbinp->GetEntry(nb);
    BinNumber = nb;
    std::cout << "#" << nb << " in " << Nbinp << std::endl;
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
    for (int it = 0; it < _Ntree; it++){
      Tdata->Add(Form(_datadir+"/out%.2d/Selected/"+selectedfile, it));
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
    TH1D * hsim = new TH1D("hsim", "hsim", 1, 0.0, M_PI);
    TH1D * h0p = new TH1D("h0p", "h0p", 180, 0.0, M_PI);
    TH1D * hnp = new TH1D("hnp", "hnp", 1, 0.0, M_PI);
    TH1D * hxp = new TH1D("hxp", "hxp", 1, 0.0, M_PI);
    TH1D * hyp = new TH1D("hyp", "hyp", 1, 0.0, M_PI);
    TH1D * hzp = new TH1D("hzp", "hzp", 1, 0.0, M_PI);
    TH1D * hQ2p = new TH1D("hQ2p", "hQ2p", 1, 0.0, M_PI);
    TH1D * hPtp = new TH1D("hPtp", "hPtp", 1, 0.0, M_PI);
    for (double ie = 0; ie < Nevent; ie = ie + 1){
      Tdata->GetEntry(ie);
      if (kin[4] < Ptl || kin[4] > Ptu) continue;
      if (kin[0] < xl || kin[0] > xu) continue;
      Lstructure::sigmaUUT(AZ, lab, sigmaN);
      Lstructure::sigmaUUTn(lab, sigman);
      hsim->Fill(std::abs(kin[5]),1);
      h0p->Fill(std::abs(kin[5]), sigmaN[0]*acc_ele*acc_pion[0]);
      hnp->Fill(std::abs(kin[5]), sigman[0]*acc_ele*acc_pion[0]);
      hxp->Fill(std::abs(kin[5]), kin[0]*sigmaN[0]*acc_ele*acc_pion[0]);
      hyp->Fill(std::abs(kin[5]), kin[1]*sigmaN[0]*acc_ele*acc_pion[0]);
      hzp->Fill(std::abs(kin[5]), kin[2]*sigmaN[0]*acc_ele*acc_pion[0]);
      hQ2p->Fill(std::abs(kin[5]), kin[3]*sigmaN[0]*acc_ele*acc_pion[0]);
      hPtp->Fill(std::abs(kin[5]), kin[4]*sigmaN[0]*acc_ele*acc_pion[0]);
    }
    Nsim = hsim->Integral(1, -1);
    Nacc = h0p->Integral(1, -1) * _eff * _lumi * _days * 24.0 * 3600.0 / _simdensity;
    std::cout << Nsim << "  " << Nacc << std::endl;
    if (Nacc > 1.0e+5 && Nsim > 10){
      fn = hnp->Integral(1,-1) / h0p->Integral(1,-1);
      xm = hxp->Integral(1,-1) / h0p->Integral(1,-1);
      ym = hyp->Integral(1,-1) / h0p->Integral(1,-1);
      zm = hzp->Integral(1,-1) / h0p->Integral(1,-1);
      Q2m = hQ2p->Integral(1,-1) / h0p->Integral(1,-1);
      Ptm = hPtp->Integral(1,-1) / h0p->Integral(1,-1);
      if (h0p->FindFirstBinAbove(1e-9) == 1) ha = 0.0;
      else ha = h0p->GetBinCenter(h0p->FindFirstBinAbove(1e-9));
      if (h0p->FindLastBinAbove(1e-9) == 180) hb = M_PI;
      else hb = h0p->GetBinCenter(h0p->FindLastBinAbove(1e-9));
      kin[0] = xm;
      kin[1] = ym;
      kin[2] = zm;
      kin[3] = Q2m;
      kin[4] = Ptm;
      Lstructure::AsinHmSN(AZ, kin, A1); 
      Lstructure::AsinHmSp(kin, A1p);
      Lstructure::AsinHmSn(kin, A1n);
      Lstructure::AsinHpSN(AZ, kin, A2); 
      Lstructure::AsinHpSp(kin, A2p);
      Lstructure::AsinHpSn(kin, A2n);
      Lstructure::Asin3HmSN(AZ, kin, A3); 
      Lstructure::Asin3HmSp(kin, A3p);
      Lstructure::Asin3HmSn(kin, A3n);
      EStatisticsUT3(h0p, Estat);
      Tp->Fill();
    }
    hsim->Delete();
    h0p->Delete();
    hnp->Delete();
    hxp->Delete();
    hyp->Delete();
    hzp->Delete();
    hQ2p->Delete();
    hPtp->Delete();
    Tdata->Delete();
  }
  std::cout << "bin analysis: pi-" << std::endl;
  //for (int nb = 0; nb < 1; nb++){
  for (int nb = 0; nb < Nbinm; nb++){
    Fbinm->GetEntry(nb);
    BinNumber = nb;
    std::cout << "#" << nb << " in " << Nbinm << std::endl;
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
    for (int it = 0; it < _Ntree; it++){
      Tdata->Add(Form(_datadir+"/out%.2d/Selected/"+selectedfile, it));
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
    TH1D * hsim = new TH1D("hsim", "hsim", 1, 0.0, M_PI);
    TH1D * h0m = new TH1D("h0m", "h0m", 180, 0.0, M_PI);
    TH1D * hnm = new TH1D("hnm", "hnm", 1, 0.0, M_PI);
    TH1D * hxm = new TH1D("hxm", "hxm", 1, 0.0, M_PI);
    TH1D * hym = new TH1D("hym", "hym", 1, 0.0, M_PI);
    TH1D * hzm = new TH1D("hzm", "hzm", 1, 0.0, M_PI);
    TH1D * hQ2m = new TH1D("hQ2m", "hQ2m", 1, 0.0, M_PI);
    TH1D * hPtm = new TH1D("hPtm", "hPtm", 1, 0.0, M_PI);
    for (double ie = 0; ie < Nevent; ie = ie + 1){
      Tdata->GetEntry(ie);
      if (kin[4] < Ptl || kin[4] > Ptu) continue;
      if (kin[0] < xl || kin[0] > xu) continue;
      Lstructure::sigmaUUT(AZ, lab, sigmaN);
      Lstructure::sigmaUUTn(lab, sigman);
      hsim->Fill(std::abs(kin[5]), 1);
      h0m->Fill(std::abs(kin[5]), sigmaN[1]*acc_ele*acc_pion[1]);
      hnm->Fill(std::abs(kin[5]), sigman[1]*acc_ele*acc_pion[1]);
      hxm->Fill(std::abs(kin[5]), kin[0]*sigmaN[1]*acc_ele*acc_pion[1]);
      hym->Fill(std::abs(kin[5]), kin[1]*sigmaN[1]*acc_ele*acc_pion[1]);
      hzm->Fill(std::abs(kin[5]), kin[2]*sigmaN[1]*acc_ele*acc_pion[1]);
      hQ2m->Fill(std::abs(kin[5]), kin[3]*sigmaN[1]*acc_ele*acc_pion[1]);
      hPtm->Fill(std::abs(kin[5]), kin[4]*sigmaN[1]*acc_ele*acc_pion[1]);
    }
    Nsim = hsim->Integral(1, -1);
    Nacc = h0m->Integral(1, -1) *_eff * _lumi * _days * 24.0 * 3600.0 / _simdensity;
    std::cout << Nsim << "  " << Nacc << std::endl;
    if (Nacc > 1.0e+5 && Nsim > 10){
      fn = hnm->Integral(1, -1) / h0m->Integral(1, -1);
      xm = hxm->Integral(1, -1) / h0m->Integral(1, -1);
      ym = hym->Integral(1, -1) / h0m->Integral(1, -1);
      zm = hzm->Integral(1, -1) / h0m->Integral(1, -1);
      Q2m = hQ2m->Integral(1, -1) / h0m->Integral(1, -1);
      Ptm = hPtm->Integral(1, -1) / h0m->Integral(1, -1);
      if (h0m->FindFirstBinAbove(1e-9) == 1) ha = 0.0;
      else ha = h0m->GetBinCenter(h0m->FindFirstBinAbove(1e-9));
      if (h0m->FindLastBinAbove(1e-9) == 180) hb = M_PI;
      else hb = h0m->GetBinCenter(h0m->FindLastBinAbove(1e-9));
      kin[0] = xm;
      kin[1] = ym;
      kin[2] = zm;
      kin[3] = Q2m;
      kin[4] = Ptm;
      Lstructure::AsinHmSN(AZ, kin, A1); 
      Lstructure::AsinHmSp(kin, A1p); 
      Lstructure::AsinHmSn(kin, A1n);
      Lstructure::AsinHpSN(AZ, kin, A2); 
      Lstructure::AsinHpSp(kin, A2p); 
      Lstructure::AsinHpSn(kin, A2n);
      Lstructure::Asin3HmSN(AZ, kin, A3); 
      Lstructure::Asin3HmSp(kin, A3p); 
      Lstructure::Asin3HmSn(kin, A3n);
      EStatisticsUT3(h0m, Estat);
      Tm->Fill();
    }
    hsim->Delete();
    h0m->Delete();
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

int Lanalysis::BinAnalysisProton(const char savefile[]){
  TFile * fbinp = new TFile(_binfile1,"r");
  if (!fbinp->IsOpen()){
    std::cout << "Lanalysis::binanalysis: Bin info file does not exist!" << std::endl;
    return 1;
  }
  TFile * fbinm = new TFile(_binfile2,"r");
  if (!fbinm->IsOpen()){
    std::cout << "Lanalysis::binanalysis: Bin info file does not exist!" << std::endl;
    return 1;
  }
  TTree * Fbinp = (TTree *) fbinp->GetObjectChecked("bin", "TTree");
  TTree * Fbinm = (TTree *) fbinm->GetObjectChecked("bin", "TTree");
  int Nbinp = Fbinp->GetEntries();
  int Nbinm = Fbinm->GetEntries();
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
  double BinNumber, Nsim;
  double fp, fNH3;
  double A1[2], A1p[2], A1NH3[2];
  double A2[2], A2p[2], A2NH3[2];
  double A3[2], A3p[2], A3NH3[2];
  double Estat[3];
  TFile * fs = new TFile(savefile,"RECREATE");
  TTree * Tp = new TTree("binplus","binplus");
  Tp->SetDirectory(fs);
  TTree * Tm = new TTree("binminus","binminus");
  Tm->SetDirectory(fs);
  Tp->Branch("BinNumber", &BinNumber, "BinNumber/D");
  Tp->Branch("Nsim", &Nsim, "Nsim/D");
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
  Tp->Branch("E1stat", &Estat[0], "E1stat/D");
  Tp->Branch("E2stat", &Estat[1], "E2stat/D");
  Tp->Branch("E3stat", &Estat[2], "E3stat/D");
  Tm->Branch("BinNumber", &BinNumber, "BinNumber/D");
  Tm->Branch("Nsim", &Nsim, "Nsim/D");
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
  Tm->Branch("E1stat", &Estat[0], "E1stat/D");
  Tm->Branch("E2stat", &Estat[1], "E2stat/D");
  Tm->Branch("E3stat", &Estat[2], "E3stat/D");
  double acc_ele, acc_pion[2];
  double sigmap[2], sigmaNH3[2], sigmaAll[2];
  double AZNH3[4] = {3.34, 0.334*7.0, -0.3, 0};
  double AZ[4] = {0.334*10.0+0.593*2.0, 0.334*7.0+0.593*2.0, 0.22, 0};
  TString nq, nz;
  TString selectedfile;
  Long64_t Nevent;
  double kin[7], lab[7];
  std::cout << "bin analysis: pi+" << std::endl;
  //for (int nb = 112; nb < 113; nb++){
  for (int nb = 0; nb < Nbinp; nb++){
    Fbinp->GetEntry(nb);
    BinNumber = nb;
    std::cout << "#" << nb << " in " << Nbinp << std::endl;
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
    for (int it = 0; it < _Ntree; it++){
      Tdata->Add(Form(_datadir+"/out%.2d/Selected/"+selectedfile, it));
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
    TH2D * hs0 = new TH2D("hs0", "hs0", 60, -M_PI, M_PI, 12, 0.0, M_PI);
    TH1D * hsim = new TH1D("hsim", "hsim", 1, -M_PI, M_PI);
    TH1D * h0p = new TH1D("h0p", "h0p", 1, -M_PI, M_PI);
    TH1D * hpp = new TH1D("hnp", "hnp", 1, -M_PI, M_PI);
    TH1D * hNp = new TH1D("hNp", "hNp", 1, -M_PI, M_PI);
    TH1D * hxp = new TH1D("hxp", "hxp", 1, -M_PI, M_PI);
    TH1D * hyp = new TH1D("hyp", "hyp", 1, -M_PI, M_PI);
    TH1D * hzp = new TH1D("hzp", "hzp", 1, -M_PI, M_PI);
    TH1D * hQ2p = new TH1D("hQ2p", "hQ2p", 1, -M_PI, M_PI);
    TH1D * hPtp = new TH1D("hPtp", "hPtp", 1, -M_PI, M_PI);
    for (int ie = 0; ie < Nevent; ie++){
      Tdata->GetEntry(ie);
      if (kin[4] < Ptl || kin[4] > Ptu) continue;
      if (kin[0] < xl || kin[0] > xu) continue;
      Lstructure::sigmaUUT(AZ, lab, sigmaAll);
      Lstructure::sigmaUUT(AZNH3, lab, sigmaNH3);
      Lstructure::sigmaUUTp(lab, sigmap);
      hsim->Fill(kin[5],1);
      hs0->Fill(kin[5], kin[6], sigmaAll[0]*acc_ele*acc_pion[0]);
      h0p->Fill(kin[5], sigmaAll[0]*acc_ele*acc_pion[0]);
      hNp->Fill(kin[5], sigmaNH3[0]*acc_ele*acc_pion[0]);
      hpp->Fill(kin[5], sigmap[0]*acc_ele*acc_pion[0]);
      hxp->Fill(kin[5], kin[0]*sigmaAll[0]*acc_ele*acc_pion[0]);
      hyp->Fill(kin[5], kin[1]*sigmaAll[0]*acc_ele*acc_pion[0]);
      hzp->Fill(kin[5], kin[2]*sigmaAll[0]*acc_ele*acc_pion[0]);
      hQ2p->Fill(kin[5], kin[3]*sigmaAll[0]*acc_ele*acc_pion[0]);
      hPtp->Fill(kin[5], kin[4]*sigmaAll[0]*acc_ele*acc_pion[0]);
    }
    Nsim = hsim->Integral(1, -1);
    Nacc = h0p->Integral(1, -1) * _eff * _lumi * _days * 24.0 * 3600.0 / _simdensity;
    std::cout << Nsim << "  " << Nacc << std::endl;
    if (Nacc > 1.0e+4){
      fp = hpp->Integral(1,-1) / h0p->Integral(1,-1);
      fNH3 = hNp->Integral(1, -1) / h0p->Integral(1, -1);
      xm = hxp->Integral(1,-1) / h0p->Integral(1,-1);
      ym = hyp->Integral(1,-1) / h0p->Integral(1,-1);
      zm = hzp->Integral(1,-1) / h0p->Integral(1,-1);
      Q2m = hQ2p->Integral(1,-1) / h0p->Integral(1,-1);
      Ptm = hPtp->Integral(1,-1) / h0p->Integral(1,-1);
      kin[0] = xm;
      kin[1] = ym;
      kin[2] = zm;
      kin[3] = Q2m;
      kin[4] = Ptm;
      Lstructure::AsinHmSN(AZ, kin, A1); 
      Lstructure::AsinHmSp(kin, A1p);
      Lstructure::AsinHmSN(AZNH3, kin, A1NH3);
      Lstructure::AsinHpSN(AZ, kin, A2); 
      Lstructure::AsinHpSp(kin, A2p);
      Lstructure::AsinHpSN(AZNH3, kin, A2NH3);
      Lstructure::Asin3HmSN(AZ, kin, A3); 
      Lstructure::Asin3HmSp(kin, A3p);
      Lstructure::Asin3HmSN(AZNH3, kin, A3NH3);
      EStatisticsUT3(hs0, Estat);
      Tp->Fill();
    }
    hs0->Delete();
    hsim->Delete();
    h0p->Delete();
    hpp->Delete();
    hNp->Delete();
    hxp->Delete();
    hyp->Delete();
    hzp->Delete();
    hQ2p->Delete();
    hPtp->Delete();
    Tdata->Delete();
  }
  std::cout << "bin analysis: pi-" << std::endl;
  //for (int nb = 0; nb < 1; nb++){
  for (int nb = 0; nb < Nbinm; nb++){
    Fbinm->GetEntry(nb);
    BinNumber = nb;
    std::cout << "#" << nb << " in " << Nbinm << std::endl;
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
    for (int it = 0; it < _Ntree; it++){
      Tdata->Add(Form(_datadir+"/out%.2d/Selected/"+selectedfile, it));
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
    TH2D * hs0 = new TH2D("hs0", "hs0", 60, -M_PI, M_PI, 12, 0.0, M_PI);
    TH1D * hsim = new TH1D("hsim", "hsim", 1, -M_PI, M_PI);
    TH1D * h0m = new TH1D("h0m", "h0m", 1, -M_PI, M_PI);
    TH1D * hpm = new TH1D("hnm", "hnm", 1, -M_PI, M_PI);
    TH1D * hNm = new TH1D("hNm", "hNm", 1, -M_PI, M_PI);
    TH1D * hxm = new TH1D("hxm", "hxm", 1, -M_PI, M_PI);
    TH1D * hym = new TH1D("hym", "hym", 1, -M_PI, M_PI);
    TH1D * hzm = new TH1D("hzm", "hzm", 1, -M_PI, M_PI);
    TH1D * hQ2m = new TH1D("hQ2m", "hQ2m", 1, -M_PI, M_PI);
    TH1D * hPtm = new TH1D("hPtm", "hPtm", 1, -M_PI, M_PI);
    for (int ie = 0; ie < Nevent; ie++){
      Tdata->GetEntry(ie);
      if (kin[4] < Ptl || kin[4] > Ptu) continue;
      if (kin[0] < xl || kin[0] > xu) continue;
      Lstructure::sigmaUUT(AZ, lab, sigmaAll);
      Lstructure::sigmaUUT(AZNH3, lab, sigmaNH3);
      Lstructure::sigmaUUTp(lab, sigmap);
      hsim->Fill(kin[5], 1);
      hs0->Fill(kin[5], kin[6], sigmaAll[1]*acc_ele*acc_pion[1]);
      h0m->Fill(kin[5], sigmaAll[1]*acc_ele*acc_pion[1]);
      hpm->Fill(kin[5], sigmap[1]*acc_ele*acc_pion[1]);
      hNm->Fill(kin[5], sigmaNH3[1]*acc_ele*acc_pion[1]);
      hxm->Fill(kin[5], kin[0]*sigmaAll[1]*acc_ele*acc_pion[1]);
      hym->Fill(kin[5], kin[1]*sigmaAll[1]*acc_ele*acc_pion[1]);
      hzm->Fill(kin[5], kin[2]*sigmaAll[1]*acc_ele*acc_pion[1]);
      hQ2m->Fill(kin[5], kin[3]*sigmaAll[1]*acc_ele*acc_pion[1]);
      hPtm->Fill(kin[5], kin[4]*sigmaAll[1]*acc_ele*acc_pion[1]);
    }
    Nsim = hsim->Integral(1, -1);
    Nacc = h0m->Integral(1, -1) *_eff * _lumi * _days * 24.0 * 3600.0 / _simdensity;
    std::cout << Nsim << "  " << Nacc << std::endl;
    if (Nacc > 1.0e+4){
      fp = hpm->Integral(1, -1) / h0m->Integral(1, -1);
      fNH3 = hNm->Integral(1, -1) / h0m->Integral(1, -1);
      xm = hxm->Integral(1, -1) / h0m->Integral(1, -1);
      ym = hym->Integral(1, -1) / h0m->Integral(1, -1);
      zm = hzm->Integral(1, -1) / h0m->Integral(1, -1);
      Q2m = hQ2m->Integral(1, -1) / h0m->Integral(1, -1);
      Ptm = hPtm->Integral(1, -1) / h0m->Integral(1, -1);
      kin[0] = xm;
      kin[1] = ym;
      kin[2] = zm;
      kin[3] = Q2m;
      kin[4] = Ptm;
      Lstructure::AsinHmSN(AZ, kin, A1); 
      Lstructure::AsinHmSp(kin, A1p); 
      Lstructure::AsinHmSN(AZNH3, kin, A1NH3);
      Lstructure::AsinHpSN(AZ, kin, A2); 
      Lstructure::AsinHpSp(kin, A2p); 
      Lstructure::AsinHpSN(AZNH3, kin, A2NH3);
      Lstructure::Asin3HmSN(AZ, kin, A3); 
      Lstructure::Asin3HmSp(kin, A3p); 
      Lstructure::Asin3HmSN(AZNH3, kin, A3NH3);
      EStatisticsUT3(hs0, Estat);
      Tm->Fill();
    }
    hs0->Delete();
    hsim->Delete();
    h0m->Delete();
    hpm->Delete();
    hNm->Delete();
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

int Lanalysis::BinResolutionNeutron(const char bintree[], const char savefile[]){
  TFile * fb = new TFile(bintree,"r");
  TTree * Tp = (TTree *) fb->Get("binplus");
  TTree * Tm = (TTree *) fb->Get("binminus");
  double BinNumber;
  double Q2m, zm;
  double xl, xu, Ptl, Ptu;
  double Nbinp = Tp->GetEntries();
  double Nbinm = Tm->GetEntries();
  Tp->SetBranchAddress("BinNumber", &BinNumber);
  Tp->SetBranchAddress("zm", &zm);
  Tp->SetBranchAddress("Q2m", &Q2m);
  Tp->SetBranchAddress("xl", &xl);
  Tp->SetBranchAddress("xu", &xu);
  Tp->SetBranchAddress("Ptl", &Ptl);
  Tp->SetBranchAddress("Ptu", &Ptu);
  Tm->SetBranchAddress("BinNumber", &BinNumber);
  Tm->SetBranchAddress("zm", &zm);
  Tm->SetBranchAddress("Q2m", &Q2m);
  Tm->SetBranchAddress("xl", &xl);
  Tm->SetBranchAddress("xu", &xu);
  Tm->SetBranchAddress("Ptl", &Ptl);
  Tm->SetBranchAddress("Ptu", &Ptu);
  double rms[7];
  double dz_e, dz_pi;
  TFile * fs = new TFile(savefile, "RECREATE");
  TTree * Rp = new TTree("rmsplus", "rmsplus");
  Rp->SetDirectory(fs);
  TTree * Rm = new TTree("rmsminus", "rmsminus");
  Rm->SetDirectory(fs);
  Rp->Branch("BinNumber", &BinNumber, "BinNumber/D");
  Rp->Branch("dx", &rms[0], "dx/D");
  Rp->Branch("dy", &rms[1], "dy/D");
  Rp->Branch("dz", &rms[2], "dz/D");
  Rp->Branch("dQ2", &rms[3], "dQ2/D");
  Rp->Branch("dPt", &rms[4], "dPt/D");
  Rp->Branch("dphih", &rms[5], "dphih/D");
  Rp->Branch("dphiS", &rms[6], "dphiS/D");
  Rp->Branch("dz_e", &dz_e, "dz_e/D");
  Rp->Branch("dz_pi", &dz_pi, "dz_pi/D");
  Rm->Branch("BinNumber", &BinNumber, "BinNumber/D");
  Rm->Branch("dx", &rms[0], "dx/D");
  Rm->Branch("dy", &rms[1], "dy/D");
  Rm->Branch("dz", &rms[2], "dz/D");
  Rm->Branch("dQ2", &rms[3], "dQ2/D");
  Rm->Branch("dPt", &rms[4], "dPt/D");
  Rm->Branch("dphih", &rms[5], "dphih/D");
  Rm->Branch("dphiS", &rms[6], "dphiS/D");
  Rm->Branch("dz_e", &dz_e, "dz_e/D");
  Rm->Branch("dz_pi", &dz_pi, "dz_pi/D");
  double acc_ele, acc_pion[2];
  double sigma[2];
  double AZ[4] = {2, 1, -0.028, 0.86};//total
  TString nq, nz;
  TString selectedfile;
  Long64_t Nevent;
  double kin[7], lab[7];
  std::cout << "bin resolution: pi+" << std::endl;
  //for (int nb = 680; nb < 681; nb++){
  for (int nb = 0; nb < Nbinp; nb++){
    Tp->GetEntry(nb);
    std::cout << "#" << nb << " in " << Nbinp << std::endl;
    if (Q2m >= 1.0 && Q2m < 2.0) nq = "1";
    else if (Q2m >= 2.0 && Q2m < 3.0) nq = "2";
    else if (Q2m >= 3.0 && Q2m < 4.0) nq = "3";
    else if (Q2m >= 4.0 && Q2m < 5.0) nq = "4";
    else if (Q2m >= 5.0 && Q2m < 6.0) nq = "5";
    else if (Q2m >= 6.0) nq = "6";
    else continue;
    if (zm >= 0.3 && zm < 0.35) nz = "30";
    else if (zm >= 0.35 && zm < 0.40) nz = "35";
    else if (zm >= 0.40 && zm < 0.45) nz = "40";
    else if (zm >= 0.45 && zm < 0.50) nz = "45";
    else if (zm >= 0.50 && zm < 0.55) nz = "50";
    else if (zm >= 0.55 && zm < 0.60) nz = "55";
    else if (zm >= 0.60 && zm < 0.65) nz = "60";
    else if (zm >= 0.65) nz = "65";
    else continue;
    selectedfile = "sidis_select"+nq+nz+".root";
    TChain * Tdata = new TChain("T", "T");
    for (int it = 0; it < _Ntree; it++){
      Tdata->Add(Form(_datadir+"/out%.2d/Selected/"+selectedfile, it));
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
    TH1D * h0 = new TH1D("h0", "h0", 1, 0.0, 2.0);
    TH1D * hx = new TH1D("hx", "hx", 1, 0.0, 2.0);
    TH1D * hy = new TH1D("hy", "hy", 1, 0.0, 2.0);
    TH1D * hz = new TH1D("hz", "hz", 1, 0.0, 2.0);
    TH1D * hQ2 = new TH1D("hQ2", "hQ2", 1, 0.0, 2.0);
    TH1D * hPt = new TH1D("hPt", "hPt", 1, 0.0, 2.0);
    TH1D * hphih = new TH1D("hphih", "hphih", 1, 0.0, 2.0);
    TH1D * hphiS = new TH1D("hphiS", "hphiS", 1, 0.0, 2.0);
    TH1D * hze = new TH1D("hze", "hze", 1.0, 0.0, 2.0);
    TH1D * hzpi = new TH1D("hzpi", "hzpi", 1.0, 0.0, 2.0);
    for (int ie = 0; ie < Nevent; ie++){
      Tdata->GetEntry(ie);
      if (kin[4] < Ptl || kin[4] > Ptu) continue;
      if (kin[0] < xl || kin[0] > xu) continue;
      Lstructure::sigmaUUT(AZ, lab, sigma);
      CalRMS(lab, rms, 100);
      dz_e = _z_e->GetBinContent(_z_e->GetXaxis()->FindBin(lab[1]),_z_e->GetYaxis()->FindBin(lab[2]));
      dz_pi = _z_pi->GetBinContent(_z_pi->GetXaxis()->FindBin(lab[4]),_z_pi->GetYaxis()->FindBin(lab[5]));
      h0->Fill(1.0, sigma[0]*acc_ele*acc_pion[0]);
      hx->Fill(1.0, rms[0]*sigma[0]*acc_ele*acc_pion[0]);
      hy->Fill(1.0, rms[1]*sigma[0]*acc_ele*acc_pion[0]);
      hz->Fill(1.0, rms[2]*sigma[0]*acc_ele*acc_pion[0]);
      hQ2->Fill(1.0, rms[3]*sigma[0]*acc_ele*acc_pion[0]);
      hPt->Fill(1.0, rms[4]*sigma[0]*acc_ele*acc_pion[0]);
      hphih->Fill(1.0, rms[5]*sigma[0]*acc_ele*acc_pion[0]);
      hphiS->Fill(1.0, rms[6]*sigma[0]*acc_ele*acc_pion[0]);
      hze->Fill(1.0, dz_e*sigma[0]*acc_ele*acc_pion[0]);
      hzpi->Fill(1.0, dz_pi*sigma[0]*acc_ele*acc_pion[0]);
    }
    rms[0] = hx->Integral(1, -1) / h0->Integral(1, -1);
    rms[1] = hy->Integral(1, -1) / h0->Integral(1, -1);
    rms[2] = hz->Integral(1, -1) / h0->Integral(1, -1);
    rms[3] = hQ2->Integral(1, -1) / h0->Integral(1, -1);
    rms[4] = hPt->Integral(1, -1) / h0->Integral(1, -1);
    rms[5] = hphih->Integral(1, -1) / h0->Integral(1, -1);
    rms[6] = hphiS->Integral(1, -1) / h0->Integral(1, -1);
    dz_e = hze->Integral(1, -1) / h0->Integral(1, -1);
    dz_pi = hzpi->Integral(1, -1) / h0->Integral(1, -1);
    Rp->Fill();
    h0->Delete();
    hx->Delete();
    hy->Delete();
    hz->Delete();
    hQ2->Delete();
    hPt->Delete();
    hphih->Delete();
    hphiS->Delete();
    hze->Delete();
    hzpi->Delete();
    Tdata->Delete();
  }
  std::cout << "bin analysis: pi-" << std::endl;
  //for (int nb = 0; nb < 1; nb++){
  for (int nb = 0; nb < Nbinm; nb++){
    Tm->GetEntry(nb);
    std::cout << "#" << nb << " in " << Nbinm << std::endl;
    if (Q2m >= 1.0 && Q2m < 2.0) nq = "1";
    else if (Q2m >= 2.0 && Q2m < 3.0) nq = "2";
    else if (Q2m >= 3.0 && Q2m < 4.0) nq = "3";
    else if (Q2m >= 4.0 && Q2m < 5.0) nq = "4";
    else if (Q2m >= 5.0 && Q2m < 6.0) nq = "5";
    else if (Q2m >= 6.0) nq = "6";
    else continue;
    if (zm >= 0.3 && zm < 0.35) nz = "30";
    else if (zm >= 0.35 && zm < 0.40) nz = "35";
    else if (zm >= 0.40 && zm < 0.45) nz = "40";
    else if (zm >= 0.45 && zm < 0.50) nz = "45";
    else if (zm >= 0.50 && zm < 0.55) nz = "50";
    else if (zm >= 0.55 && zm < 0.60) nz = "55";
    else if (zm >= 0.60 && zm < 0.65) nz = "60";
    else if (zm >= 0.65) nz = "65";
    else continue;
    selectedfile = "sidis_select"+nq+nz+".root";
    TChain * Tdata = new TChain("T", "T");
    for (int it = 0; it < _Ntree; it++){
      Tdata->Add(Form(_datadir+"/out%.2d/Selected/"+selectedfile, it));
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
    TH1D * h0 = new TH1D("h0", "h0", 1, 0.0, 2.0);
    TH1D * hx = new TH1D("hx", "hx", 1, 0.0, 2.0);
    TH1D * hy = new TH1D("hy", "hy", 1, 0.0, 2.0);
    TH1D * hz = new TH1D("hz", "hz", 1, 0.0, 2.0);
    TH1D * hQ2 = new TH1D("hQ2", "hQ2", 1, 0.0, 2.0);
    TH1D * hPt = new TH1D("hPt", "hPt", 1, 0.0, 2.0);
    TH1D * hphih = new TH1D("hphih", "hphih", 1, 0.0, 2.0);
    TH1D * hphiS = new TH1D("hphiS", "hphiS", 1, 0.0, 2.0);
    TH1D * hze = new TH1D("hze", "hze", 1.0, 0.0, 2.0);
    TH1D * hzpi = new TH1D("hzpi", "hzpi", 1.0, 0.0, 2.0);
    for (int ie = 0; ie < Nevent; ie++){
      Tdata->GetEntry(ie);
      if (kin[4] < Ptl || kin[4] > Ptu) continue;
      if (kin[0] < xl || kin[0] > xu) continue;
      Lstructure::sigmaUUT(AZ, lab, sigma);
      CalRMS(lab, rms, 100);
      dz_e = _z_e->GetBinContent(_z_e->GetXaxis()->FindBin(lab[1]),_z_e->GetYaxis()->FindBin(lab[2]));
      dz_pi = _z_pi->GetBinContent(_z_pi->GetXaxis()->FindBin(lab[4]),_z_pi->GetYaxis()->FindBin(lab[5]));
      h0->Fill(1.0, sigma[1]*acc_ele*acc_pion[1]);
      hx->Fill(1.0, rms[0]*sigma[1]*acc_ele*acc_pion[1]);
      hy->Fill(1.0, rms[1]*sigma[1]*acc_ele*acc_pion[1]);
      hz->Fill(1.0, rms[2]*sigma[1]*acc_ele*acc_pion[1]);
      hQ2->Fill(1.0, rms[3]*sigma[1]*acc_ele*acc_pion[1]);
      hPt->Fill(1.0, rms[4]*sigma[1]*acc_ele*acc_pion[1]);
      hphih->Fill(1.0, rms[5]*sigma[1]*acc_ele*acc_pion[1]);
      hphiS->Fill(1.0, rms[6]*sigma[1]*acc_ele*acc_pion[1]);
      hze->Fill(1.0, dz_e*sigma[1]*acc_ele*acc_pion[1]);
      hzpi->Fill(1.0, dz_pi*sigma[1]*acc_ele*acc_pion[1]);
    }
    rms[0] = hx->Integral(1, -1) / h0->Integral(1, -1);
    rms[1] = hy->Integral(1, -1) / h0->Integral(1, -1);
    rms[2] = hz->Integral(1, -1) / h0->Integral(1, -1);
    rms[3] = hQ2->Integral(1, -1) / h0->Integral(1, -1);
    rms[4] = hPt->Integral(1, -1) / h0->Integral(1, -1);
    rms[5] = hphih->Integral(1, -1) / h0->Integral(1, -1);
    rms[6] = hphiS->Integral(1, -1) / h0->Integral(1, -1);
    dz_e = hze->Integral(1, -1) / h0->Integral(1, -1);
    dz_pi = hzpi->Integral(1, -1) / h0->Integral(1, -1);
    Rm->Fill();
    h0->Delete();
    hx->Delete();
    hy->Delete();
    hz->Delete();
    hQ2->Delete();
    hPt->Delete();
    hphih->Delete();
    hphiS->Delete();
    hze->Delete();
    hzpi->Delete();
    Tdata->Delete();
  }
  fs->Write();
  return 0;
}

int Lanalysis::BinAcceptanceNeutron(const char bintree[], const char savefile[]){
  TFile * fb = new TFile(bintree, "r");
  TTree * Tp = (TTree *) fb->Get("binplus");
  double Nbinp = Tp->GetEntries();
  TTree * Tm = (TTree *) fb->Get("binminus");
  double Nbinm = Tm->GetEntries();
  double BinNumber;
  double zm, Q2m, xl, xu, Ptl, Ptu;
  Tp->SetBranchAddress("BinNumber", &BinNumber);
  Tp->SetBranchAddress("zm", &zm);
  Tp->SetBranchAddress("Q2m", &Q2m);
  Tp->SetBranchAddress("xl", &xl);
  Tp->SetBranchAddress("xu", &xu);
  Tp->SetBranchAddress("Ptl", &Ptl);
  Tp->SetBranchAddress("Ptu", &Ptu);
  Tm->SetBranchAddress("BinNumber", &BinNumber);
  Tm->SetBranchAddress("zm", &zm);
  Tm->SetBranchAddress("Q2m", &Q2m);
  Tm->SetBranchAddress("xl", &xl);
  Tm->SetBranchAddress("xu", &xu);
  Tm->SetBranchAddress("Ptl", &Ptl);
  Tm->SetBranchAddress("Ptu", &Ptu);
  TFile * fs = new TFile(savefile, "RECREATE");
  TTree * Ap = new TTree("accplus", "accplus");
  Ap->SetDirectory(fs);
  TTree * Am = new TTree("accminus", "accminus");
  Am->SetDirectory(fs);
  double cc[6], range[2];
  Ap->Branch("BinNumber", &BinNumber, "BinNumber/D");
  Ap->Branch("ha", &range[0], "ha/D");
  Ap->Branch("hb", &range[1], "hb/D");
  Ap->Branch("c0", &cc[0], "c0/D");
  Ap->Branch("c1", &cc[1], "c1/D");
  Ap->Branch("c2", &cc[2], "c2/D");
  Ap->Branch("c3", &cc[3], "c3/D");
  Ap->Branch("c4", &cc[4], "c4/D");
  Ap->Branch("c5", &cc[5], "c5/D");
  Am->Branch("BinNumber", &BinNumber, "BinNumber/D");
  Am->Branch("ha", &range[0], "ha/D");
  Am->Branch("hb", &range[1], "hb/D");
  Am->Branch("c0", &cc[0], "c0/D");
  Am->Branch("c1", &cc[1], "c1/D");
  Am->Branch("c2", &cc[2], "c2/D");
  Am->Branch("c3", &cc[3], "c3/D");
  Am->Branch("c4", &cc[4], "c4/D");
  Am->Branch("c5", &cc[5], "c5/D");
  double acc_ele, acc_pion[2];
  double sigma[2];
  double AZ[4] = {2, 1, -0.028, 0.86};//total
  TString nq, nz;
  TString selectedfile;
  Long64_t Nevent;
  double lab[7];
  double x, Pt, phih;
  std::cout << "Bin acceptance: pi+" << std::endl;
  //for (int nb = 12; nb < 13; nb++){
  for (int nb = 0; nb < Nbinp; nb++){
    Tp->GetEntry(nb);
    std::cout << "#" << nb << " in " << Nbinp << std::endl;
    if (Q2m >= 1.0 && Q2m < 2.0) nq = "1";
    else if (Q2m >= 2.0 && Q2m < 3.0) nq = "2";
    else if (Q2m >= 3.0 && Q2m < 4.0) nq = "3";
    else if (Q2m >= 4.0 && Q2m < 5.0) nq = "4";
    else if (Q2m >= 5.0 && Q2m < 6.0) nq = "5";
    else if (Q2m >= 6.0) nq = "6";
    else continue;
    if (zm >= 0.3 && zm < 0.35) nz = "30";
    else if (zm >= 0.35 && zm < 0.40) nz = "35";
    else if (zm >= 0.40 && zm < 0.45) nz = "40";
    else if (zm >= 0.45 && zm < 0.50) nz = "45";
    else if (zm >= 0.50 && zm < 0.55) nz = "50";
    else if (zm >= 0.55 && zm < 0.60) nz = "55";
    else if (zm >= 0.60 && zm < 0.65) nz = "60";
    else if (zm >= 0.65) nz = "65";
    else continue;
    selectedfile = "sidis_select"+nq+nz+".root";
    TChain * Tdata = new TChain("T", "T");
    for (int it = 0; it < _Ntree; it++){
      Tdata->Add(Form(_datadir+"/out%.2d/Selected/"+selectedfile, it));
    }
    Nevent = Tdata->GetEntries();
    Tdata->SetBranchAddress("x", &x);
    Tdata->SetBranchAddress("Pt", &Pt);
    Tdata->SetBranchAddress("Ebeam", &lab[0]);
    Tdata->SetBranchAddress("p_ele", &lab[1]);
    Tdata->SetBranchAddress("theta_ele", &lab[2]);
    Tdata->SetBranchAddress("phi_ele", &lab[3]);
    Tdata->SetBranchAddress("p_pion", &lab[4]);
    Tdata->SetBranchAddress("theta_pion", &lab[5]);
    Tdata->SetBranchAddress("phi_pion", &lab[6]);
    Tdata->SetBranchAddress("phi_h", &phih);
    Tdata->SetBranchAddress("acc_ele", &acc_ele);
    Tdata->SetBranchAddress("acc_pion_p", &acc_pion[0]);
    Tdata->SetBranchAddress("acc_pion_m", &acc_pion[1]);
    TH1D * h0 = new TH1D("h0", "h0", 360, 0.0, M_PI);
    for (double ie = 0; ie < Nevent; ie = ie + 1){
      Tdata->GetEntry(ie);
      if (Pt < Ptl || Pt > Ptu) continue;
      if (x < xl || x > xu) continue;
      Lstructure::sigmaUUT(AZ, lab, sigma);
      h0->Fill(std::abs(phih), sigma[0]*acc_ele*acc_pion[0]);
    }
    range[0] = h0->GetBinCenter(h0->FindFirstBinAbove());
    range[1] = h0->GetBinCenter(h0->FindLastBinAbove());
    FitCos1(h0, range, cc);
    Ap->Fill();
    h0->Delete();
    Tdata->Delete();
  }
  std::cout << "Bin acceptance: pi-" << std::endl;
  //for (int nb = 12; nb < 13; nb++){
  for (int nb = 0; nb < Nbinm; nb++){
    Tp->GetEntry(nb);
    std::cout << "#" << nb << " in " << Nbinm << std::endl;
    if (Q2m >= 1.0 && Q2m < 2.0) nq = "1";
    else if (Q2m >= 2.0 && Q2m < 3.0) nq = "2";
    else if (Q2m >= 3.0 && Q2m < 4.0) nq = "3";
    else if (Q2m >= 4.0 && Q2m < 5.0) nq = "4";
    else if (Q2m >= 5.0 && Q2m < 6.0) nq = "5";
    else if (Q2m >= 6.0) nq = "6";
    else continue;
    if (zm >= 0.3 && zm < 0.35) nz = "30";
    else if (zm >= 0.35 && zm < 0.40) nz = "35";
    else if (zm >= 0.40 && zm < 0.45) nz = "40";
    else if (zm >= 0.45 && zm < 0.50) nz = "45";
    else if (zm >= 0.50 && zm < 0.55) nz = "50";
    else if (zm >= 0.55 && zm < 0.60) nz = "55";
    else if (zm >= 0.60 && zm < 0.65) nz = "60";
    else if (zm >= 0.65) nz = "65";
    else continue;
    selectedfile = "sidis_select"+nq+nz+".root";
    TChain * Tdata = new TChain("T", "T");
    for (int it = 0; it < _Ntree; it++){
      Tdata->Add(Form(_datadir+"/out%.2d/Selected/"+selectedfile, it));
    }
    Nevent = Tdata->GetEntries();
    Tdata->SetBranchAddress("x", &x);
    Tdata->SetBranchAddress("Pt", &Pt);
    Tdata->SetBranchAddress("Ebeam", &lab[0]);
    Tdata->SetBranchAddress("p_ele", &lab[1]);
    Tdata->SetBranchAddress("theta_ele", &lab[2]);
    Tdata->SetBranchAddress("phi_ele", &lab[3]);
    Tdata->SetBranchAddress("p_pion", &lab[4]);
    Tdata->SetBranchAddress("theta_pion", &lab[5]);
    Tdata->SetBranchAddress("phi_pion", &lab[6]);
    Tdata->SetBranchAddress("phi_h", &phih);
    Tdata->SetBranchAddress("acc_ele", &acc_ele);
    Tdata->SetBranchAddress("acc_pion_p", &acc_pion[0]);
    Tdata->SetBranchAddress("acc_pion_m", &acc_pion[1]);
    TH1D * h0 = new TH1D("h0", "h0", 360, 0.0, M_PI);
    for (int ie = 0; ie < Nevent; ie++){
      Tdata->GetEntry(ie);
      if (Pt < Ptl || Pt > Ptu) continue;
      if (x < xl || x > xu) continue;
      Lstructure::sigmaUUT(AZ, lab, sigma);
      h0->Fill(std::abs(phih), sigma[1]*acc_ele*acc_pion[1]);
    }
    range[0] = h0->GetBinCenter(h0->FindFirstBinAbove());
    range[1] = h0->GetBinCenter(h0->FindLastBinAbove());
    FitCos1(h0, range, cc);
    Am->Fill();
    h0->Delete();
    Tdata->Delete();
  }
  fs->Write();
  return 0;
}

int Lanalysis::BinAcceptanceProton(const char bintree[], const char savefile[]){
  TFile * fb = new TFile(bintree, "r");
  TTree * Tp = (TTree *) fb->Get("binplus");
  double Nbinp = Tp->GetEntries();
  TTree * Tm = (TTree *) fb->Get("binminus");
  double Nbinm = Tm->GetEntries();
  double BinNumber;
  double zm, Q2m, xl, xu, Ptl, Ptu;
  Tp->SetBranchAddress("BinNumber", &BinNumber);
  Tp->SetBranchAddress("zm", &zm);
  Tp->SetBranchAddress("Q2m", &Q2m);
  Tp->SetBranchAddress("xl", &xl);
  Tp->SetBranchAddress("xu", &xu);
  Tp->SetBranchAddress("Ptl", &Ptl);
  Tp->SetBranchAddress("Ptu", &Ptu);
  Tm->SetBranchAddress("BinNumber", &BinNumber);
  Tm->SetBranchAddress("zm", &zm);
  Tm->SetBranchAddress("Q2m", &Q2m);
  Tm->SetBranchAddress("xl", &xl);
  Tm->SetBranchAddress("xu", &xu);
  Tm->SetBranchAddress("Ptl", &Ptl);
  Tm->SetBranchAddress("Ptu", &Ptu);
  TFile * fs = new TFile(savefile, "RECREATE");
  TTree * Ap = new TTree("accplus", "accplus");
  Ap->SetDirectory(fs);
  TTree * Am = new TTree("accminus", "accminus");
  Am->SetDirectory(fs);
  double cc[11], range[2];
  Ap->Branch("BinNumber", &BinNumber, "BinNumber/D");
  Ap->Branch("ha", &range[0], "ha/D");
  Ap->Branch("hb", &range[1], "hb/D");
  Ap->Branch("c0", &cc[0], "c0/D");
  Ap->Branch("c1", &cc[1], "c1/D");
  Ap->Branch("c2", &cc[2], "c2/D");
  Ap->Branch("c3", &cc[3], "c3/D");
  Ap->Branch("c4", &cc[4], "c4/D");
  Ap->Branch("c5", &cc[5], "c5/D");
  Ap->Branch("c6", &cc[6], "c6/D");
  Ap->Branch("c7", &cc[7], "c7/D");
  Ap->Branch("c8", &cc[8], "c8/D");
  Ap->Branch("c9", &cc[9], "c9/D");
  Ap->Branch("c10", &cc[10], "c10/D");
  Am->Branch("BinNumber", &BinNumber, "BinNumber/D");
  Am->Branch("ha", &range[0], "ha/D");
  Am->Branch("hb", &range[1], "hb/D");
  Am->Branch("c0", &cc[0], "c0/D");
  Am->Branch("c1", &cc[1], "c1/D");
  Am->Branch("c2", &cc[2], "c2/D");
  Am->Branch("c3", &cc[3], "c3/D");
  Am->Branch("c4", &cc[4], "c4/D");
  Am->Branch("c5", &cc[5], "c5/D");
  Am->Branch("c6", &cc[6], "c6/D");
  Am->Branch("c7", &cc[7], "c7/D");
  Am->Branch("c8", &cc[8], "c8/D");
  Am->Branch("c9", &cc[9], "c9/D");
  Am->Branch("c10", &cc[10], "c10/D");
  double acc_ele, acc_pion[2];
  double sigma[2];
  double AZ[4] = {0.334*10.0+0.593*2.0, 0.334*7.0+0.593*2.0, 0.22, 0};//total
  TString nq, nz;
  TString selectedfile;
  Long64_t Nevent;
  double lab[7];
  double x, Pt, phih;
  std::cout << "Bin acceptance: pi+" << std::endl;
  //for (int nb = 12; nb < 13; nb++){
  for (int nb = 0; nb < Nbinp; nb++){
    Tp->GetEntry(nb);
    std::cout << "#" << nb << " in " << Nbinp << std::endl;
    if (Q2m >= 1.0 && Q2m < 2.0) nq = "1";
    else if (Q2m >= 2.0 && Q2m < 3.0) nq = "2";
    else if (Q2m >= 3.0 && Q2m < 4.0) nq = "3";
    else if (Q2m >= 4.0 && Q2m < 5.0) nq = "4";
    else if (Q2m >= 5.0 && Q2m < 6.0) nq = "5";
    else if (Q2m >= 6.0) nq = "6";
    else continue;
    if (zm >= 0.3 && zm < 0.35) nz = "30";
    else if (zm >= 0.35 && zm < 0.40) nz = "35";
    else if (zm >= 0.40 && zm < 0.45) nz = "40";
    else if (zm >= 0.45 && zm < 0.50) nz = "45";
    else if (zm >= 0.50 && zm < 0.55) nz = "50";
    else if (zm >= 0.55 && zm < 0.60) nz = "55";
    else if (zm >= 0.60 && zm < 0.65) nz = "60";
    else if (zm >= 0.65) nz = "65";
    else continue;
    selectedfile = "sidis_select"+nq+nz+".root";
    TChain * Tdata = new TChain("T", "T");
    for (int it = 0; it < _Ntree; it++){
      Tdata->Add(Form(_datadir+"/out%.2d/Selected/"+selectedfile, it));
    }
    Nevent = Tdata->GetEntries();
    Tdata->SetBranchAddress("x", &x);
    Tdata->SetBranchAddress("Pt", &Pt);
    Tdata->SetBranchAddress("Ebeam", &lab[0]);
    Tdata->SetBranchAddress("p_ele", &lab[1]);
    Tdata->SetBranchAddress("theta_ele", &lab[2]);
    Tdata->SetBranchAddress("phi_ele", &lab[3]);
    Tdata->SetBranchAddress("p_pion", &lab[4]);
    Tdata->SetBranchAddress("theta_pion", &lab[5]);
    Tdata->SetBranchAddress("phi_pion", &lab[6]);
    Tdata->SetBranchAddress("phi_h", &phih);
    Tdata->SetBranchAddress("acc_ele", &acc_ele);
    Tdata->SetBranchAddress("acc_pion_p", &acc_pion[0]);
    Tdata->SetBranchAddress("acc_pion_m", &acc_pion[1]);
    TH1D * h0 = new TH1D("h0", "h0", 360, -M_PI, M_PI);
    for (double ie = 0; ie < Nevent; ie = ie + 1){
      Tdata->GetEntry(ie);
      if (Pt < Ptl || Pt > Ptu) continue;
      if (x < xl || x > xu) continue;
      Lstructure::sigmaUUT(AZ, lab, sigma);
      h0->Fill(phih, sigma[0]*acc_ele*acc_pion[0]);
    }
    range[0] = h0->GetBinCenter(h0->FindFirstBinAbove());
    range[1] = h0->GetBinCenter(h0->FindLastBinAbove());
    FitCosSin1(h0, range, cc);
    Ap->Fill();
    h0->Delete();
    Tdata->Delete();
  }
  std::cout << "Bin acceptance: pi-" << std::endl;
  //for (int nb = 12; nb < 13; nb++){
  for (int nb = 0; nb < Nbinm; nb++){
    Tp->GetEntry(nb);
    std::cout << "#" << nb << " in " << Nbinm << std::endl;
    if (Q2m >= 1.0 && Q2m < 2.0) nq = "1";
    else if (Q2m >= 2.0 && Q2m < 3.0) nq = "2";
    else if (Q2m >= 3.0 && Q2m < 4.0) nq = "3";
    else if (Q2m >= 4.0 && Q2m < 5.0) nq = "4";
    else if (Q2m >= 5.0 && Q2m < 6.0) nq = "5";
    else if (Q2m >= 6.0) nq = "6";
    else continue;
    if (zm >= 0.3 && zm < 0.35) nz = "30";
    else if (zm >= 0.35 && zm < 0.40) nz = "35";
    else if (zm >= 0.40 && zm < 0.45) nz = "40";
    else if (zm >= 0.45 && zm < 0.50) nz = "45";
    else if (zm >= 0.50 && zm < 0.55) nz = "50";
    else if (zm >= 0.55 && zm < 0.60) nz = "55";
    else if (zm >= 0.60 && zm < 0.65) nz = "60";
    else if (zm >= 0.65) nz = "65";
    else continue;
    selectedfile = "sidis_select"+nq+nz+".root";
    TChain * Tdata = new TChain("T", "T");
    for (int it = 0; it < _Ntree; it++){
      Tdata->Add(Form(_datadir+"/out%.2d/Selected/"+selectedfile, it));
    }
    Nevent = Tdata->GetEntries();
    Tdata->SetBranchAddress("x", &x);
    Tdata->SetBranchAddress("Pt", &Pt);
    Tdata->SetBranchAddress("Ebeam", &lab[0]);
    Tdata->SetBranchAddress("p_ele", &lab[1]);
    Tdata->SetBranchAddress("theta_ele", &lab[2]);
    Tdata->SetBranchAddress("phi_ele", &lab[3]);
    Tdata->SetBranchAddress("p_pion", &lab[4]);
    Tdata->SetBranchAddress("theta_pion", &lab[5]);
    Tdata->SetBranchAddress("phi_pion", &lab[6]);
    Tdata->SetBranchAddress("phi_h", &phih);
    Tdata->SetBranchAddress("acc_ele", &acc_ele);
    Tdata->SetBranchAddress("acc_pion_p", &acc_pion[0]);
    Tdata->SetBranchAddress("acc_pion_m", &acc_pion[1]);
    TH1D * h0 = new TH1D("h0", "h0", 360, -M_PI, M_PI);
    for (int ie = 0; ie < Nevent; ie++){
      Tdata->GetEntry(ie);
      if (Pt < Ptl || Pt > Ptu) continue;
      if (x < xl || x > xu) continue;
      Lstructure::sigmaUUT(AZ, lab, sigma);
      h0->Fill(phih, sigma[1]*acc_ele*acc_pion[1]);
    }
    range[0] = h0->GetBinCenter(h0->FindFirstBinAbove());
    range[1] = h0->GetBinCenter(h0->FindLastBinAbove());
    FitCosSin1(h0, range, cc);
    Am->Fill();
    h0->Delete();
    Tdata->Delete();
  }
  fs->Write();
  return 0;
}

int Lanalysis::ECoincidenceNeutron(const char bintree[], const char rmstree[], const char savefile[]){
  double dt = 6.0e-9;
  double fv = 1.0;
  double scale_e = 1.0;
  double scale_pi = 1.0;
  TFile * fb = new TFile(bintree, "r");
  TTree * Tp = (TTree *) fb->Get("binplus");
  double Nbinp = Tp->GetEntries();
  TTree * Tm = (TTree *) fb->Get("binminus");
  double Nbinm = Tm->GetEntries();
  double BinNumber;
  double zm, Q2m, xl, xu, Ptl, Ptu;
  double Nacc;
  Tp->SetBranchAddress("BinNumber", &BinNumber);
  Tp->SetBranchAddress("zm", &zm);
  Tp->SetBranchAddress("Q2m", &Q2m);
  Tp->SetBranchAddress("xl", &xl);
  Tp->SetBranchAddress("xu", &xu);
  Tp->SetBranchAddress("Ptl", &Ptl);
  Tp->SetBranchAddress("Ptu", &Ptu);
  Tp->SetBranchAddress("Nacc", &Nacc);
  Tm->SetBranchAddress("BinNumber", &BinNumber);
  Tm->SetBranchAddress("zm", &zm);
  Tm->SetBranchAddress("Q2m", &Q2m);
  Tm->SetBranchAddress("xl", &xl);
  Tm->SetBranchAddress("xu", &xu);
  Tm->SetBranchAddress("Ptl", &Ptl);
  Tm->SetBranchAddress("Ptu", &Ptu);
  Tm->SetBranchAddress("Nacc", &Nacc);
  TFile * fr = new TFile(rmstree, "r");
  TTree * Rp = (TTree *) fr->Get("rmsplus");
  TTree * Rm = (TTree *) fr->Get("rmsminus");
  double dz_e, dz_pi;
  double safe = 1.5;
  Rp->SetBranchAddress("dz_e", &dz_e);
  Rp->SetBranchAddress("dz_pi", &dz_pi);
  Rm->SetBranchAddress("dz_e", &dz_e);
  Rm->SetBranchAddress("dz_pi", &dz_pi);
  TFile * fs = new TFile(savefile, "RECREATE");
  TTree * Cp = new TTree("coinplus", "coinplus");
  Cp->SetDirectory(fs);
  TTree * Cm = new TTree("coinminus", "coinminus");
  Cm->SetDirectory(fs);
  double Ncoin, SB, ErrRel, Lv;
  Cp->Branch("BinNumber", &BinNumber, "BinNumber/D");
  Cp->Branch("Nacc", &Nacc, "Nacc/D");
  Cp->Branch("Ncoin", &Ncoin, "Ncoin/D");
  Cp->Branch("SB", &SB, "SB/D");
  Cp->Branch("Lv", &Lv, "Lv/D");
  Cp->Branch("fv", &fv, "fv/D");
  Cp->Branch("ErrRel", &ErrRel, "ErrRel/D");
  Cm->Branch("BinNumber", &BinNumber, "BinNumber/D");
  Cm->Branch("Nacc", &Nacc, "Nacc/D");
  Cm->Branch("Ncoin", &Ncoin, "Ncoin/D");
  Cm->Branch("SB", &SB, "SB/D");
  Cm->Branch("Lv", &Lv, "Lv/D");
  Cm->Branch("fv", &fv, "fv/D");
  Cm->Branch("ErrRel", &ErrRel, "ErrRel/D");
  double acc_ele, acc_pion[2];
  double sigma[2];
  double AZ[4] = {2, 1, -0.028, 0.86};//total
  TString nq, nz;
  TString selectedfile;
  Long64_t Nevent;
  double lab[7];
  double x, Pt;
  std::cout << "Random coincidence: pi+" << std::endl;
  //for (int nb = 12; nb < 13; nb++){
  for (int nb = 0; nb < Nbinp; nb++){
    Tp->GetEntry(nb);
    Rp->GetEntry(nb);
    std::cout << "#" << nb << " in " << Nbinp << std::endl;
    if (Q2m >= 1.0 && Q2m < 2.0) nq = "1";
    else if (Q2m >= 2.0 && Q2m < 3.0) nq = "2";
    else if (Q2m >= 3.0 && Q2m < 4.0) nq = "3";
    else if (Q2m >= 4.0 && Q2m < 5.0) nq = "4";
    else if (Q2m >= 5.0 && Q2m < 6.0) nq = "5";
    else if (Q2m >= 6.0) nq = "6";
    else continue;
    if (zm >= 0.3 && zm < 0.35) nz = "30";
    else if (zm >= 0.35 && zm < 0.40) nz = "35";
    else if (zm >= 0.40 && zm < 0.45) nz = "40";
    else if (zm >= 0.45 && zm < 0.50) nz = "45";
    else if (zm >= 0.50 && zm < 0.55) nz = "50";
    else if (zm >= 0.55 && zm < 0.60) nz = "55";
    else if (zm >= 0.60 && zm < 0.65) nz = "60";
    else if (zm >= 0.65) nz = "65";
    else continue;
    selectedfile = "sidis_select"+nq+nz+".root";
    TChain * Tdata = new TChain("T", "T");
    for (int it = 0; it < _Ntree; it++){
      Tdata->Add(Form(_datadir+"/out%.2d/Selected/"+selectedfile, it));
    }
    Nevent = Tdata->GetEntries();
    Tdata->SetBranchAddress("x", &x);
    Tdata->SetBranchAddress("Pt", &Pt);
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
    TH1D * h0 = new TH1D("h0", "h0", 1, 0.0, 2.0);
    for (int ie = 0; ie < Nevent; ie++){
      Tdata->GetEntry(ie);
      if (Pt < Ptl || Pt > Ptu) continue;
      if (x < xl || x > xu) continue;
      RandomCoincidenceSigmaN(AZ, lab, sigma);
      h0->Fill(1.0, sigma[0]*acc_ele*acc_pion[0]*dt);
    }
    Lv = 3.0 * (dz_e + dz_pi) * safe;
    fv = 40.0 / Lv;
    Ncoin = (h0->Integral(1, -1)) * _eff * _lumi * _days * 24.0 * 3600.0 / _simdensity;
    SB = Nacc / (Ncoin / fv);
    ErrRel = sqrt(scale_e * scale_pi * Ncoin / fv / fv) / (Nacc + Ncoin / fv); 
    Cp->Fill();
    h0->Delete();
    Tdata->Delete();
  }
  std::cout << "Random coincidence: pi-" << std::endl;
  //for (int nb = 0; nb < 1; nb++){
  for (int nb = 0; nb < Nbinm; nb++){
    Tm->GetEntry(nb);
    Rm->GetEntry(nb);
    std::cout << "#" << nb << " in " << Nbinm << std::endl;
    if (Q2m >= 1.0 && Q2m < 2.0) nq = "1";
    else if (Q2m >= 2.0 && Q2m < 3.0) nq = "2";
    else if (Q2m >= 3.0 && Q2m < 4.0) nq = "3";
    else if (Q2m >= 4.0 && Q2m < 5.0) nq = "4";
    else if (Q2m >= 5.0 && Q2m < 6.0) nq = "5";
    else if (Q2m >= 6.0) nq = "6";
    else continue;
    if (zm >= 0.3 && zm < 0.35) nz = "30";
    else if (zm >= 0.35 && zm < 0.40) nz = "35";
    else if (zm >= 0.40 && zm < 0.45) nz = "40";
    else if (zm >= 0.45 && zm < 0.50) nz = "45";
    else if (zm >= 0.50 && zm < 0.55) nz = "50";
    else if (zm >= 0.55 && zm < 0.60) nz = "55";
    else if (zm >= 0.60 && zm < 0.65) nz = "60";
    else if (zm >= 0.65) nz = "65";
    else continue;
    selectedfile = "sidis_select"+nq+nz+".root";
    TChain * Tdata = new TChain("T", "T");
    for (int it = 0; it < _Ntree; it++){
      Tdata->Add(Form(_datadir+"/out%.2d/Selected/"+selectedfile, it));
    }
    Nevent = Tdata->GetEntries();
    Tdata->SetBranchAddress("x", &x);
    Tdata->SetBranchAddress("Pt", &Pt);
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
    TH1D * h0 = new TH1D("h0", "h0", 1, 0.0, 2.0);
    for (int ie = 0; ie < Nevent; ie++){
      Tdata->GetEntry(ie);
      if (Pt < Ptl || Pt > Ptu) continue;
      if (x < xl || x > xu) continue;
      RandomCoincidenceSigmaN(AZ, lab, sigma);
      h0->Fill(1.0, sigma[1]*acc_ele*acc_pion[1]*dt);
    }
    Lv = 3.0 * (dz_e + dz_pi) * safe;
    fv = 40.0 / Lv;
    Ncoin = (h0->Integral(1, -1)) * _eff * _lumi * _days * 24.0 * 3600.0 / _simdensity;
    SB = Nacc / (Ncoin / fv);
    ErrRel = sqrt(scale_e * scale_pi * Ncoin / fv / fv) / (Nacc + Ncoin / fv); 
    Cm->Fill();
    h0->Delete();
    Tdata->Delete();
  }
  fs->Write();
  return 0;
}

int Lanalysis::ECoincidenceProton(const char bintree[], const char rmstree[], const char savefile[]){
  double dt = 6.0e-9;
  double fv = 1.0;
  double scale_e = 1.0;
  double scale_pi = 1.0;
  TFile * fb = new TFile(bintree, "r");
  TTree * Tp = (TTree *) fb->Get("binplus");
  double Nbinp = Tp->GetEntries();
  TTree * Tm = (TTree *) fb->Get("binminus");
  double Nbinm = Tm->GetEntries();
  double BinNumber;
  double zm, Q2m, xl, xu, Ptl, Ptu;
  double Nacc;
  Tp->SetBranchAddress("BinNumber", &BinNumber);
  Tp->SetBranchAddress("zm", &zm);
  Tp->SetBranchAddress("Q2m", &Q2m);
  Tp->SetBranchAddress("xl", &xl);
  Tp->SetBranchAddress("xu", &xu);
  Tp->SetBranchAddress("Ptl", &Ptl);
  Tp->SetBranchAddress("Ptu", &Ptu);
  Tp->SetBranchAddress("Nacc", &Nacc);
  Tm->SetBranchAddress("BinNumber", &BinNumber);
  Tm->SetBranchAddress("zm", &zm);
  Tm->SetBranchAddress("Q2m", &Q2m);
  Tm->SetBranchAddress("xl", &xl);
  Tm->SetBranchAddress("xu", &xu);
  Tm->SetBranchAddress("Ptl", &Ptl);
  Tm->SetBranchAddress("Ptu", &Ptu);
  Tm->SetBranchAddress("Nacc", &Nacc);
  TFile * fs = new TFile(savefile, "RECREATE");
  TTree * Cp = new TTree("coinplus", "coinplus");
  Cp->SetDirectory(fs);
  TTree * Cm = new TTree("coinminus", "coinminus");
  Cm->SetDirectory(fs);
  double Ncoin, SB, ErrRel, Lv = 6.0;
  Cp->Branch("BinNumber", &BinNumber, "BinNumber/D");
  Cp->Branch("Nacc", &Nacc, "Nacc/D");
  Cp->Branch("Ncoin", &Ncoin, "Ncoin/D");
  Cp->Branch("SB", &SB, "SB/D");
  Cp->Branch("Lv", &Lv, "Lv/D");
  Cp->Branch("fv", &fv, "fv/D");
  Cp->Branch("ErrRel", &ErrRel, "ErrRel/D");
  Cm->Branch("BinNumber", &BinNumber, "BinNumber/D");
  Cm->Branch("Nacc", &Nacc, "Nacc/D");
  Cm->Branch("Ncoin", &Ncoin, "Ncoin/D");
  Cm->Branch("SB", &SB, "SB/D");
  Cm->Branch("Lv", &Lv, "Lv/D");
  Cm->Branch("fv", &fv, "fv/D");
  Cm->Branch("ErrRel", &ErrRel, "ErrRel/D");
  double acc_ele, acc_pion[2];
  double sigma[2];
  double AZ[4] = {0.334*10.0+0.593*2.0, 0.334*7.0+0.593*2.0, 0.22, 0};//total
  TString nq, nz;
  TString selectedfile;
  Long64_t Nevent;
  double lab[7];
  double x, Pt;
  std::cout << "Random coincidence: pi+" << std::endl;
  //for (int nb = 12; nb < 13; nb++){
  for (int nb = 0; nb < Nbinp; nb++){
    Tp->GetEntry(nb);
    std::cout << "#" << nb << " in " << Nbinp << std::endl;
    if (Q2m >= 1.0 && Q2m < 2.0) nq = "1";
    else if (Q2m >= 2.0 && Q2m < 3.0) nq = "2";
    else if (Q2m >= 3.0 && Q2m < 4.0) nq = "3";
    else if (Q2m >= 4.0 && Q2m < 5.0) nq = "4";
    else if (Q2m >= 5.0 && Q2m < 6.0) nq = "5";
    else if (Q2m >= 6.0) nq = "6";
    else continue;
    if (zm >= 0.3 && zm < 0.35) nz = "30";
    else if (zm >= 0.35 && zm < 0.40) nz = "35";
    else if (zm >= 0.40 && zm < 0.45) nz = "40";
    else if (zm >= 0.45 && zm < 0.50) nz = "45";
    else if (zm >= 0.50 && zm < 0.55) nz = "50";
    else if (zm >= 0.55 && zm < 0.60) nz = "55";
    else if (zm >= 0.60 && zm < 0.65) nz = "60";
    else if (zm >= 0.65) nz = "65";
    else continue;
    selectedfile = "sidis_select"+nq+nz+".root";
    TChain * Tdata = new TChain("T", "T");
    for (int it = 0; it < _Ntree; it++){
      Tdata->Add(Form(_datadir+"/out%.2d/Selected/"+selectedfile, it));
    }
    Nevent = Tdata->GetEntries();
    Tdata->SetBranchAddress("x", &x);
    Tdata->SetBranchAddress("Pt", &Pt);
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
    TH1D * h0 = new TH1D("h0", "h0", 1, 0.0, 2.0);
    for (int ie = 0; ie < Nevent; ie++){
      Tdata->GetEntry(ie);
      if (Pt < Ptl || Pt > Ptu) continue;
      if (x < xl || x > xu) continue;
      RandomCoincidenceSigmaN(AZ, lab, sigma);
      h0->Fill(1.0, sigma[0]*acc_ele*acc_pion[0]*dt);
    }
    fv = 1.0;
    Ncoin = (h0->Integral(1, -1)) * _eff * _lumi * _days * 24.0 * 3600.0 / _simdensity;
    SB = Nacc / (Ncoin / fv);
    ErrRel = sqrt(scale_e * scale_pi * Ncoin / fv / fv) / (Nacc + Ncoin / fv); 
    Cp->Fill();
    h0->Delete();
    Tdata->Delete();
  }
  std::cout << "Random coincidence: pi-" << std::endl;
  //for (int nb = 0; nb < 1; nb++){
  for (int nb = 0; nb < Nbinm; nb++){
    Tm->GetEntry(nb);
    std::cout << "#" << nb << " in " << Nbinm << std::endl;
    if (Q2m >= 1.0 && Q2m < 2.0) nq = "1";
    else if (Q2m >= 2.0 && Q2m < 3.0) nq = "2";
    else if (Q2m >= 3.0 && Q2m < 4.0) nq = "3";
    else if (Q2m >= 4.0 && Q2m < 5.0) nq = "4";
    else if (Q2m >= 5.0 && Q2m < 6.0) nq = "5";
    else if (Q2m >= 6.0) nq = "6";
    else continue;
    if (zm >= 0.3 && zm < 0.35) nz = "30";
    else if (zm >= 0.35 && zm < 0.40) nz = "35";
    else if (zm >= 0.40 && zm < 0.45) nz = "40";
    else if (zm >= 0.45 && zm < 0.50) nz = "45";
    else if (zm >= 0.50 && zm < 0.55) nz = "50";
    else if (zm >= 0.55 && zm < 0.60) nz = "55";
    else if (zm >= 0.60 && zm < 0.65) nz = "60";
    else if (zm >= 0.65) nz = "65";
    else continue;
    selectedfile = "sidis_select"+nq+nz+".root";
    TChain * Tdata = new TChain("T", "T");
    for (int it = 0; it < _Ntree; it++){
      Tdata->Add(Form(_datadir+"/out%.2d/Selected/"+selectedfile, it));
    }
    Nevent = Tdata->GetEntries();
    Tdata->SetBranchAddress("x", &x);
    Tdata->SetBranchAddress("Pt", &Pt);
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
    TH1D * h0 = new TH1D("h0", "h0", 1, 0.0, 2.0);
    for (int ie = 0; ie < Nevent; ie++){
      Tdata->GetEntry(ie);
      if (Pt < Ptl || Pt > Ptu) continue;
      if (x < xl || x > xu) continue;
      RandomCoincidenceSigmaN(AZ, lab, sigma);
      h0->Fill(1.0, sigma[1]*acc_ele*acc_pion[1]*dt);
    }
    fv = 1.0;
    Ncoin = (h0->Integral(1, -1)) * _eff * _lumi * _days * 24.0 * 3600.0 / _simdensity;
    SB = Nacc / (Ncoin / fv);
    ErrRel = sqrt(scale_e * scale_pi * Ncoin / fv / fv) / (Nacc + Ncoin / fv); 
    Cm->Fill();
    h0->Delete();
    Tdata->Delete();
  }
  fs->Write();
  return 0;
}

int Lanalysis::EResolutionNeutron(const char bintree[], const char rmstree[], const char acctree[], const char savefile[]){
  TFile * fb = new TFile(bintree, "r");
  TTree * Tp = (TTree *) fb->Get("binplus");
  double Nbinp = Tp->GetEntries();
  TTree * Tm = (TTree *) fb->Get("binminus");
  double Nbinm = Tm->GetEntries();
  double BinNumber;
  double A1, A2, A3;
  Tp->SetBranchAddress("BinNumber", &BinNumber);
  Tp->SetBranchAddress("A1", &A1);
  Tp->SetBranchAddress("A2", &A2);
  Tp->SetBranchAddress("A3", &A3);
  Tp->SetBranchAddress("BinNumber", &BinNumber);
  Tp->SetBranchAddress("A1", &A1);
  Tp->SetBranchAddress("A2", &A2);
  Tp->SetBranchAddress("A3", &A3);
  TFile * fr = new TFile(rmstree, "r");
  TTree * Rp = (TTree *) fr->Get("rmsplus");
  TTree * Rm = (TTree *) fr->Get("rmsminus");
  double dphih, dphiS;
  Rp->SetBranchAddress("dphih", &dphih);
  Rp->SetBranchAddress("dphiS", &dphiS);
  Rm->SetBranchAddress("dphih", &dphih);
  Rm->SetBranchAddress("dphiS", &dphiS);
  TFile * fa = new TFile(acctree, "r");
  TTree * Ap = (TTree *) fa->Get("accplus");
  TTree * Am = (TTree *) fa->Get("accminus");
  double ha, hb;
  double cc[6];
  Ap->SetBranchAddress("ha", &ha);
  Ap->SetBranchAddress("hb", &hb);
  Ap->SetBranchAddress("c0", &cc[0]);
  Ap->SetBranchAddress("c1", &cc[1]);
  Ap->SetBranchAddress("c2", &cc[2]);
  Ap->SetBranchAddress("c3", &cc[3]);
  Ap->SetBranchAddress("c4", &cc[4]);
  Ap->SetBranchAddress("c5", &cc[5]);
  Am->SetBranchAddress("ha", &ha);
  Am->SetBranchAddress("hb", &hb);
  Am->SetBranchAddress("c0", &cc[0]);
  Am->SetBranchAddress("c1", &cc[1]);
  Am->SetBranchAddress("c2", &cc[2]);
  Am->SetBranchAddress("c3", &cc[3]);
  Am->SetBranchAddress("c4", &cc[4]);
  Am->SetBranchAddress("c5", &cc[5]);
  TFile * fs = new TFile(savefile, "RECREATE");
  TTree * Sp = new TTree("Rp", "Rp");
  Sp->SetDirectory(fs);
  TTree * Sm = new TTree("Rm", "Rm");
  Sm->SetDirectory(fs);
  double Asym[3], ErrRelH[3], ErrRelS[3];
  Sp->Branch("BinNumber", &BinNumber, "BinNumber/D");
  Sp->Branch("ErrRelH1", &ErrRelH[0], "ErrRelH1/D");
  Sp->Branch("ErrRelH2", &ErrRelH[1], "ErrRelH2/D");
  Sp->Branch("ErrRelH3", &ErrRelH[2], "ErrRelH3/D");
  Sp->Branch("ErrRelS1", &ErrRelS[0], "ErrRelS1/D");
  Sp->Branch("ErrRelS2", &ErrRelS[1], "ErrRelS2/D");
  Sp->Branch("ErrRelS3", &ErrRelS[2], "ErrRelS3/D");
  Sm->Branch("BinNumber", &BinNumber, "BinNumber/D");
  Sm->Branch("ErrRelH1", &ErrRelH[0], "ErrRelH1/D");
  Sm->Branch("ErrRelH2", &ErrRelH[1], "ErrRelH2/D");
  Sm->Branch("ErrRelH3", &ErrRelH[2], "ErrRelH3/D");
  Sm->Branch("ErrRelS1", &ErrRelS[0], "ErrRelS1/D");
  Sm->Branch("ErrRelS2", &ErrRelS[1], "ErrRelS2/D");
  Sm->Branch("ErrRelS3", &ErrRelS[2], "ErrRelS3/D");
  TF1 * f1 = new TF1("f1", "[0]*(1.0+[1]*cos(x)+[2]*cos(2.0*x)+[3]*cos(3.0*x)+[4]*cos(4.0*x)+[5]*cos(5.0*x))", -M_PI, M_PI);
  TF2 * f2 = new TF2("f2", "[0]*sin(x-y)+[1]*sin(x+y)+[2]*sin(3.0*x-y)", -M_PI, M_PI, 0.0, M_PI);
  double phih, phiS, sum1, sum2;
  TF1 * rs;
  std::cout << "Resolution uncertianty pi+:" << std::endl;
  //for (int nb = 55; nb < 60; nb++){
  for (int nb = 0; nb < Nbinp; nb++){
    Tp->GetEntry(nb);
    Rp->GetEntry(nb);
    Ap->GetEntry(nb);
    std::cout << "#" << nb << " in " << Nbinp << std::endl;
    f1->SetParameters(1.0, cc[1], cc[2], cc[3], cc[4], cc[5]);
    f2->SetParameters(A1, A2, A3);
    TH2D * ht1 = new TH2D("ht1", "ht1", 720, -M_PI, M_PI, 360, 0.0, M_PI);
    TH2D * ht2 = new TH2D("ht2", "ht2", 720, -M_PI, M_PI, 360, 0.0, M_PI);
    TH2D * hh1 = new TH2D("hh1", "hh1", 720, -M_PI, M_PI, 360, 0.0, M_PI);
    TH2D * hh2 = new TH2D("hh2", "hh2", 720, -M_PI, M_PI, 360, 0.0, M_PI);
    TH2D * hs1 = new TH2D("hs1", "hs1", 720, -M_PI, M_PI, 360, 0.0, M_PI);
    TH2D * hs2 = new TH2D("hs2", "hs2", 720, -M_PI, M_PI, 360, 0.0, M_PI);
    TH2D * hh = new TH2D("hh", "hh", 720, -M_PI, M_PI, 360, 0.0, M_PI);
    TH2D * hs = new TH2D("hs", "hs", 720, -M_PI, M_PI, 360, 0.0, M_PI);
    for (int i = 0; i <= 720; i++){
      for (int j = 1; j <= 360; j++){
	phih = hs->GetXaxis()->GetBinCenter(i);
	if (std::abs(phih) < ha || std::abs(phih) > hb) continue;
	phiS = hs->GetYaxis()->GetBinCenter(j);
	ht1->SetBinContent(i, j, (1.0+f2->Eval(phih, phiS))*f1->Eval(phih));
	ht2->SetBinContent(i, j, (1.0-f2->Eval(phih, phiS))*f1->Eval(phih));	
      }
    }
    for (int i = 0; i <= 720; i++){
      for (int j = 1; j <=360; j++){
	phih = hs->GetXaxis()->GetBinCenter(i);
	if (std::abs(phih) < ha || std::abs(phih) > hb) continue;
	phiS = hs->GetYaxis()->GetBinCenter(j);
	sum1 = 0.0;
	sum2 = 0.0;
	for (int k = 1; k <= 720; k++){
	  sum1 = sum1 + ht1->GetBinContent(k, j) * TMath::Gaus(phih, ht1->GetXaxis()->GetBinCenter(k), dphih);
	  sum2 = sum2 + ht2->GetBinContent(k, j) * TMath::Gaus(phih, ht2->GetXaxis()->GetBinCenter(k), dphih);
	}
	hh1->SetBinContent(i, j, sum1);
	hh2->SetBinContent(i, j, sum2);
      }
    }
    for (int i = 1; i <= 720; i++){
      for (int j = 1; j <=360; j++){
	phih = hs->GetXaxis()->GetBinCenter(i);
	if (std::abs(phih) < ha || std::abs(phih) > hb) continue;
	phiS = hs->GetYaxis()->GetBinCenter(j);
	sum1 = 0.0;
	sum2 = 0.0;
	for (int k = 1; k <= 360; k++){
	  sum1 = sum1 + ht1->GetBinContent(i, k) * TMath::Gaus(phiS, ht1->GetYaxis()->GetBinCenter(k), dphiS);
	  sum2 = sum2 + ht2->GetBinContent(i, k) * TMath::Gaus(phiS, ht2->GetYaxis()->GetBinCenter(k), dphiS);
	}
	hs1->SetBinContent(i, j, sum1);
	hs2->SetBinContent(i, j, sum2);
      }
    }
    for (int i = 1; i <= 720; i++){
      for (int j = 1; j <=360; j++){
	phih = hs->GetXaxis()->GetBinCenter(i);
	if (std::abs(phih) < ha || std::abs(phih) > hb) continue;
	hh->SetBinContent(i, j, (hh1->GetBinContent(i,j)-hh2->GetBinContent(i,j))/(hh1->GetBinContent(i,j)+hh2->GetBinContent(i,j)));
	hs->SetBinContent(i, j, (hs1->GetBinContent(i,j)-hs2->GetBinContent(i,j))/(hs1->GetBinContent(i,j)+hs2->GetBinContent(i,j)));
      }
    }
    hh->Fit("f2","QO", "", ha, hb);
    rs = hh->GetFunction("f2");
    Asym[0] = rs->GetParameter(0);
    Asym[1] = rs->GetParameter(1);
    Asym[2] = rs->GetParameter(2);
    ErrRelH[0] = std::abs((A1-Asym[0]) / A1);
    ErrRelH[1] = std::abs((A2-Asym[1]) / A2);
    ErrRelH[2] = std::abs((A3-Asym[2]) / A3);
    hs->Fit("f2","QO", "", ha, hb);
    rs = hs->GetFunction("f2");
    Asym[0] = rs->GetParameter(0);
    Asym[1] = rs->GetParameter(1);
    Asym[2] = rs->GetParameter(2);
    ErrRelS[0] = std::abs((A1-Asym[0]) / A1);
    ErrRelS[1] = std::abs((A2-Asym[1]) / A2);
    ErrRelS[2] = std::abs((A3-Asym[2]) / A3);
    Sp->Fill();
    ht1->Delete();
    ht2->Delete();
    hh1->Delete();
    hh2->Delete();
    hs1->Delete();
    hs2->Delete();
    hh->Delete();
    hs->Delete();
  }
  std::cout << "Resolution uncertianty pi-:" << std::endl;
  //for (int nb = 55; nb < 60; nb++){
  for (int nb = 0; nb < Nbinm; nb++){
    Tm->GetEntry(nb);
    Rm->GetEntry(nb);
    Am->GetEntry(nb);
    std::cout << "#" << nb << " in " << Nbinm << std::endl;
    f1->SetParameters(1.0, cc[1], cc[2], cc[3], cc[4], cc[5]);
    f2->SetParameters(A1, A2, A3);
    TH2D * ht1 = new TH2D("ht1", "ht1", 720, -M_PI, M_PI, 360, 0.0, M_PI);
    TH2D * ht2 = new TH2D("ht2", "ht2", 720, -M_PI, M_PI, 360, 0.0, M_PI);
    TH2D * hh1 = new TH2D("hh1", "hh1", 720, -M_PI, M_PI, 360, 0.0, M_PI);
    TH2D * hh2 = new TH2D("hh2", "hh2", 720, -M_PI, M_PI, 360, 0.0, M_PI);
    TH2D * hs1 = new TH2D("hs1", "hs1", 720, -M_PI, M_PI, 360, 0.0, M_PI);
    TH2D * hs2 = new TH2D("hs2", "hs2", 720, -M_PI, M_PI, 360, 0.0, M_PI);
    TH2D * hh = new TH2D("hh", "hh", 720, -M_PI, M_PI, 360, 0.0, M_PI);
    TH2D * hs = new TH2D("hs", "hs", 720, -M_PI, M_PI, 360, 0.0, M_PI);
    for (int i = 0; i <= 720; i++){
      for (int j = 1; j <= 360; j++){
	phih = hs->GetXaxis()->GetBinCenter(i);
	if (std::abs(phih) < ha || std::abs(phih) > hb) continue;
	phiS = hs->GetYaxis()->GetBinCenter(j);
	ht1->SetBinContent(i, j, (1.0+f2->Eval(phih, phiS))*f1->Eval(phih));
	ht2->SetBinContent(i, j, (1.0-f2->Eval(phih, phiS))*f1->Eval(phih));	
      }
    }
    for (int i = 0; i <= 720; i++){
      for (int j = 1; j <=360; j++){
	phih = hs->GetXaxis()->GetBinCenter(i);
	if (std::abs(phih) < ha || std::abs(phih) > hb) continue;
	phiS = hs->GetYaxis()->GetBinCenter(j);
	sum1 = 0.0;
	sum2 = 0.0;
	for (int k = 1; k <= 720; k++){
	  sum1 = sum1 + ht1->GetBinContent(k, j) * TMath::Gaus(phih, ht1->GetXaxis()->GetBinCenter(k), dphih);
	  sum2 = sum2 + ht2->GetBinContent(k, j) * TMath::Gaus(phih, ht2->GetXaxis()->GetBinCenter(k), dphih);
	}
	hh1->SetBinContent(i, j, sum1);
	hh2->SetBinContent(i, j, sum2);
      }
    }
    for (int i = 1; i <= 720; i++){
      for (int j = 1; j <=360; j++){
	phih = hs->GetXaxis()->GetBinCenter(i);
	if (std::abs(phih) < ha || std::abs(phih) > hb) continue;
	phiS = hs->GetYaxis()->GetBinCenter(j);
	sum1 = 0.0;
	sum2 = 0.0;
	for (int k = 1; k <= 360; k++){
	  sum1 = sum1 + ht1->GetBinContent(i, k) * TMath::Gaus(phiS, ht1->GetYaxis()->GetBinCenter(k), dphiS);
	  sum2 = sum2 + ht2->GetBinContent(i, k) * TMath::Gaus(phiS, ht2->GetYaxis()->GetBinCenter(k), dphiS);
	}
	hs1->SetBinContent(i, j, sum1);
	hs2->SetBinContent(i, j, sum2);
      }
    }
    for (int i = 1; i <= 720; i++){
      for (int j = 1; j <=360; j++){
	phih = hs->GetXaxis()->GetBinCenter(i);
	if (std::abs(phih) < ha || std::abs(phih) > hb) continue;
	hh->SetBinContent(i, j, (hh1->GetBinContent(i,j)-hh2->GetBinContent(i,j))/(hh1->GetBinContent(i,j)+hh2->GetBinContent(i,j)));
	hs->SetBinContent(i, j, (hs1->GetBinContent(i,j)-hs2->GetBinContent(i,j))/(hs1->GetBinContent(i,j)+hs2->GetBinContent(i,j)));
      }
    }
    hh->Fit("f2","QO", "", ha, hb);
    rs = hh->GetFunction("f2");
    Asym[0] = rs->GetParameter(0);
    Asym[1] = rs->GetParameter(1);
    Asym[2] = rs->GetParameter(2);
    ErrRelH[0] = std::abs((A1-Asym[0]) / A1);
    ErrRelH[1] = std::abs((A2-Asym[1]) / A2);
    ErrRelH[2] = std::abs((A3-Asym[2]) / A3);
    hs->Fit("f2","QO", "", ha, hb);
    rs = hs->GetFunction("f2");
    Asym[0] = rs->GetParameter(0);
    Asym[1] = rs->GetParameter(1);
    Asym[2] = rs->GetParameter(2);
    ErrRelS[0] = std::abs((A1-Asym[0]) / A1);
    ErrRelS[1] = std::abs((A2-Asym[1]) / A2);
    ErrRelS[2] = std::abs((A3-Asym[2]) / A3);
    Sm->Fill();
    ht1->Delete();
    ht2->Delete();
    hh1->Delete();
    hh2->Delete();
    hs1->Delete();
    hs2->Delete();
    hh->Delete();
    hs->Delete();
  }
  fs->Write();
  return 0;
}

int Lanalysis::ENuclearPDF(const char bintree[], const char savefile[]){
  TFile * fb = new TFile(bintree, "r");
  TTree * Tp = (TTree *) fb->Get("binplus");
  double Nbinp = Tp->GetEntries();
  TTree * Tm = (TTree *) fb->Get("binminus");
  double Nbinm = Tm->GetEntries();
  double BinNumber;
  double zm, Q2m, xl, xu, Ptl, Ptu;
  double Nacc;
  Tp->SetBranchAddress("BinNumber", &BinNumber);
  Tp->SetBranchAddress("zm", &zm);
  Tp->SetBranchAddress("Q2m", &Q2m);
  Tp->SetBranchAddress("xl", &xl);
  Tp->SetBranchAddress("xu", &xu);
  Tp->SetBranchAddress("Ptl", &Ptl);
  Tp->SetBranchAddress("Ptu", &Ptu);
  Tp->SetBranchAddress("Nacc", &Nacc);
  Tm->SetBranchAddress("BinNumber", &BinNumber);
  Tm->SetBranchAddress("zm", &zm);
  Tm->SetBranchAddress("Q2m", &Q2m);
  Tm->SetBranchAddress("xl", &xl);
  Tm->SetBranchAddress("xu", &xu);
  Tm->SetBranchAddress("Ptl", &Ptl);
  Tm->SetBranchAddress("Ptu", &Ptu);
  Tm->SetBranchAddress("Nacc", &Nacc);
  TFile * fs = new TFile(savefile, "RECREATE");
  TTree * Np = new TTree("npdfplus", "npdfplus");
  Np->SetDirectory(fs);
  TTree * Nm = new TTree("npdfminus", "npdfminus");
  Nm->SetDirectory(fs);
  double Npfree, Nnfree, NNfree;
  double Npbound, Nnbound, NNbound;
  double ErrRelp, ErrReln, ErrRelN;
  Np->Branch("BinNumber", &BinNumber, "BinNumber/D");
  Np->Branch("Nacc", &Nacc, "Nacc/D");
  Np->Branch("Npfree", &Npfree, "Npfree/D");
  Np->Branch("Npbound", &Npbound, "Npbound/D");
  Np->Branch("Nnfree", &Nnfree, "Nnfree/D");
  Np->Branch("Nnbound", &Nnbound, "Nnbound/D");
  Np->Branch("NNfree", &NNfree, "NNfree/D");
  Np->Branch("NNbound", &NNbound, "NNbound/D");
  Np->Branch("ErrRelp", &ErrRelp, "ErrRelp/D");
  Np->Branch("ErrReln", &ErrReln, "ErrReln/D");
  Np->Branch("ErrRelN", &ErrRelN, "ErrRelN/D");
  Nm->Branch("BinNumber", &BinNumber, "BinNumber/D");
  Nm->Branch("Nacc", &Nacc, "Nacc/D");
  Nm->Branch("Npfree", &Npfree, "Npfree/D");
  Nm->Branch("Npbound", &Npbound, "Npbound/D");
  Nm->Branch("Nnfree", &Nnfree, "Nnfree/D");
  Nm->Branch("Nnbound", &Nnbound, "Nnbound/D");
  Nm->Branch("NNfree", &NNfree, "NNfree/D");
  Nm->Branch("NNbound", &NNbound, "NNbound/D");
  Nm->Branch("ErrRelp", &ErrRelp, "ErrRelp/D");
  Nm->Branch("ErrReln", &ErrReln, "ErrReln/D");
  Nm->Branch("ErrRelN", &ErrRelN, "ErrRelN/D");
  double acc_ele, acc_pion[2];
  double sigma[2];
  double AZ[4] = {2, 1, -0.028, 0.86};//total
  TString nq, nz;
  TString selectedfile;
  Long64_t Nevent;
  double lab[7];
  double x, Pt;
  std::cout << "Nuclear PDF: pi+" << std::endl;
  //for (int nb = 12; nb < 13; nb++){
  for (int nb = 0; nb < Nbinp; nb++){
    Tp->GetEntry(nb);
    std::cout << "#" << nb << " in " << Nbinp << std::endl;
    if (Q2m >= 1.0 && Q2m < 2.0) nq = "1";
    else if (Q2m >= 2.0 && Q2m < 3.0) nq = "2";
    else if (Q2m >= 3.0 && Q2m < 4.0) nq = "3";
    else if (Q2m >= 4.0 && Q2m < 5.0) nq = "4";
    else if (Q2m >= 5.0 && Q2m < 6.0) nq = "5";
    else if (Q2m >= 6.0) nq = "6";
    else continue;
    if (zm >= 0.3 && zm < 0.35) nz = "30";
    else if (zm >= 0.35 && zm < 0.40) nz = "35";
    else if (zm >= 0.40 && zm < 0.45) nz = "40";
    else if (zm >= 0.45 && zm < 0.50) nz = "45";
    else if (zm >= 0.50 && zm < 0.55) nz = "50";
    else if (zm >= 0.55 && zm < 0.60) nz = "55";
    else if (zm >= 0.60 && zm < 0.65) nz = "60";
    else if (zm >= 0.65) nz = "65";
    else continue;
    selectedfile = "sidis_select"+nq+nz+".root";
    TChain * Tdata = new TChain("T", "T");
    for (int it = 0; it < _Ntree; it++){
      Tdata->Add(Form(_datadir+"/out%.2d/Selected/"+selectedfile, it));
    }
    Nevent = Tdata->GetEntries();
    Tdata->SetBranchAddress("x", &x);
    Tdata->SetBranchAddress("Pt", &Pt);
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
    TH1D * hfp = new TH1D("hfp", "hfp", 1, 0.0, 2.0);
    TH1D * hbp = new TH1D("hbp", "hbp", 1, 0.0, 2.0);
    TH1D * hfn = new TH1D("hfn", "hfn", 1, 0.0, 2.0);
    TH1D * hbn = new TH1D("hbn", "hbn", 1, 0.0, 2.0);
    TH1D * hfN = new TH1D("hfN", "hfN", 1, 0.0, 2.0);
    TH1D * hbN = new TH1D("hbN", "hbN", 1, 0.0, 2.0);
    for (int ie = 0; ie < Nevent; ie++){
      Tdata->GetEntry(ie);
      if (Pt < Ptl || Pt > Ptu) continue;
      if (x < xl || x > xu) continue;
      Lstructure::sigmaUUTp_bound(lab, sigma, "hydrogen-1");
      hfp->Fill(1.0, sigma[0]*acc_ele*acc_pion[0]);
      Lstructure::sigmaUUTp_bound(lab, sigma, "helium-3");
      hbp->Fill(1.0, sigma[0]*acc_ele*acc_pion[0]);
      Lstructure::sigmaUUTn_bound(lab, sigma, "hydrogen-1");
      hfn->Fill(1.0, sigma[0]*acc_ele*acc_pion[0]);
      Lstructure::sigmaUUTn_bound(lab, sigma, "helium-3");
      hbn->Fill(1.0, sigma[0]*acc_ele*acc_pion[0]);
      Lstructure::sigmaUUT_bound(AZ, lab, sigma, "hydrogen-1");
      hfN->Fill(1.0, sigma[0]*acc_ele*acc_pion[0]);
      Lstructure::sigmaUUT_bound(AZ, lab, sigma, "helium-3");
      hbN->Fill(1.0, sigma[0]*acc_ele*acc_pion[0]);
    }
    Npfree = 2.0 * (hfp->Integral(1, -1)) * _eff * _lumi * _days * 24.0 * 3600.0 / _simdensity;
    Npbound = 2.0 * (hbp->Integral(1, -1)) * _eff * _lumi * _days * 24.0 * 3600.0 / _simdensity;
    Nnfree = (hfn->Integral(1, -1)) * _eff * _lumi * _days * 24.0 * 3600.0 / _simdensity;
    Nnbound = (hbn->Integral(1, -1)) * _eff * _lumi * _days * 24.0 * 3600.0 / _simdensity;
    NNfree = (hfN->Integral(1, -1)) * _eff * _lumi * _days * 24.0 * 3600.0 / _simdensity;
    NNbound = (hbN->Integral(1, -1)) * _eff * _lumi * _days * 24.0 * 3600.0 / _simdensity;
    ErrRelp = std::abs(Npbound - Npfree) / Npbound;
    ErrReln = std::abs(Nnbound - Nnfree) / Nnbound;
    ErrRelN = std::abs(NNbound - NNfree) / NNbound;
    Np->Fill();
    hfp->Delete();
    hbp->Delete();
    hfn->Delete();
    hbn->Delete();
    hfN->Delete();
    hbN->Delete();
    Tdata->Delete();
  }
  std::cout << "Nuclear PDF: pi-" << std::endl;
  //for (int nb = 12; nb < 13; nb++){
  for (int nb = 0; nb < Nbinm; nb++){
    Tm->GetEntry(nb);
    std::cout << "#" << nb << " in " << Nbinm << std::endl;
    if (Q2m >= 1.0 && Q2m < 2.0) nq = "1";
    else if (Q2m >= 2.0 && Q2m < 3.0) nq = "2";
    else if (Q2m >= 3.0 && Q2m < 4.0) nq = "3";
    else if (Q2m >= 4.0 && Q2m < 5.0) nq = "4";
    else if (Q2m >= 5.0 && Q2m < 6.0) nq = "5";
    else if (Q2m >= 6.0) nq = "6";
    else continue;
    if (zm >= 0.3 && zm < 0.35) nz = "30";
    else if (zm >= 0.35 && zm < 0.40) nz = "35";
    else if (zm >= 0.40 && zm < 0.45) nz = "40";
    else if (zm >= 0.45 && zm < 0.50) nz = "45";
    else if (zm >= 0.50 && zm < 0.55) nz = "50";
    else if (zm >= 0.55 && zm < 0.60) nz = "55";
    else if (zm >= 0.60 && zm < 0.65) nz = "60";
    else if (zm >= 0.65) nz = "65";
    else continue;
    selectedfile = "sidis_select"+nq+nz+".root";
    TChain * Tdata = new TChain("T", "T");
    for (int it = 0; it < _Ntree; it++){
      Tdata->Add(Form(_datadir+"/out%.2d/Selected/"+selectedfile, it));
    }
    Nevent = Tdata->GetEntries();
    Tdata->SetBranchAddress("x", &x);
    Tdata->SetBranchAddress("Pt", &Pt);
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
    TH1D * hfp = new TH1D("hfp", "hfp", 1, 0.0, 2.0);
    TH1D * hbp = new TH1D("hbp", "hbp", 1, 0.0, 2.0);
    TH1D * hfn = new TH1D("hfn", "hfn", 1, 0.0, 2.0);
    TH1D * hbn = new TH1D("hbn", "hbn", 1, 0.0, 2.0);
    TH1D * hfN = new TH1D("hfN", "hfN", 1, 0.0, 2.0);
    TH1D * hbN = new TH1D("hbN", "hbN", 1, 0.0, 2.0);
    for (int ie = 0; ie < Nevent; ie++){
      Tdata->GetEntry(ie);
      if (Pt < Ptl || Pt > Ptu) continue;
      if (x < xl || x > xu) continue;
      Lstructure::sigmaUUTp_bound(lab, sigma, "hydrogen-1");
      hfp->Fill(1.0, sigma[1]*acc_ele*acc_pion[1]);
      Lstructure::sigmaUUTp_bound(lab, sigma, "helium-3");
      hbp->Fill(1.0, sigma[1]*acc_ele*acc_pion[1]);
      Lstructure::sigmaUUTn_bound(lab, sigma, "hydrogen-1");
      hfn->Fill(1.0, sigma[1]*acc_ele*acc_pion[1]);
      Lstructure::sigmaUUTn_bound(lab, sigma, "helium-3");
      hbn->Fill(1.0, sigma[1]*acc_ele*acc_pion[1]);
      Lstructure::sigmaUUT_bound(AZ, lab, sigma, "hydrogen-1");
      hfN->Fill(1.0, sigma[1]*acc_ele*acc_pion[1]);
      Lstructure::sigmaUUT_bound(AZ, lab, sigma, "helium-3");
      hbN->Fill(1.0, sigma[1]*acc_ele*acc_pion[1]);
    }
    Npfree = 2.0 * (hfp->Integral(1, -1)) * _eff * _lumi * _days * 24.0 * 3600.0 / _simdensity;
    Npbound = 2.0 * (hbp->Integral(1, -1)) * _eff * _lumi * _days * 24.0 * 3600.0 / _simdensity;
    Nnfree = (hfn->Integral(1, -1)) * _eff * _lumi * _days * 24.0 * 3600.0 / _simdensity;
    Nnbound = (hbn->Integral(1, -1)) * _eff * _lumi * _days * 24.0 * 3600.0 / _simdensity;
    NNfree = (hfN->Integral(1, -1)) * _eff * _lumi * _days * 24.0 * 3600.0 / _simdensity;
    NNbound = (hbN->Integral(1, -1)) * _eff * _lumi * _days * 24.0 * 3600.0 / _simdensity;
    ErrRelp = std::abs(Npbound - Npfree) / Npbound;
    ErrReln = std::abs(Nnbound - Nnfree) / Nnbound;
    ErrRelN = std::abs(NNbound - NNfree) / NNbound;
    Nm->Fill();
    hfp->Delete();
    hbp->Delete();
    hfn->Delete();
    hbn->Delete();
    hfN->Delete();
    hbN->Delete();
    Tdata->Delete();
  }
  fs->Write();
  return 0;
}

int Lanalysis::ENuclearNeutron(const char bintree[], const char savefile[]){
  double BinNumber;
  double At[3], Ap[3], An[3];
  double fn, Nacc;
  double ftime = _days;
  TFile * fb = new TFile(bintree, "r");
  TTree * Tp = (TTree *) fb->Get("binplus");
  double Nbinp = Tp->GetEntries();
  TTree * Tm = (TTree *) fb->Get("binminus");
  double Nbinm = Tm->GetEntries();
  Tp->SetBranchAddress("BinNumber", &BinNumber);
  Tp->SetBranchAddress("A1", &At[0]);
  Tp->SetBranchAddress("A2", &At[1]);
  Tp->SetBranchAddress("A3", &At[2]);
  Tp->SetBranchAddress("A1p", &Ap[0]);
  Tp->SetBranchAddress("A2p", &Ap[1]);
  Tp->SetBranchAddress("A3p", &Ap[2]);
  Tp->SetBranchAddress("A1n", &An[0]);
  Tp->SetBranchAddress("A2n", &An[1]);
  Tp->SetBranchAddress("A3n", &An[2]);
  Tp->SetBranchAddress("fn", &fn);
  Tp->SetBranchAddress("Nacc", &Nacc);
  Tm->SetBranchAddress("BinNumber", &BinNumber);
  Tm->SetBranchAddress("A1", &At[0]);
  Tm->SetBranchAddress("A2", &At[1]);
  Tm->SetBranchAddress("A3", &At[2]);
  Tm->SetBranchAddress("A1p", &Ap[0]);
  Tm->SetBranchAddress("A2p", &Ap[1]);
  Tm->SetBranchAddress("A3p", &Ap[2]);
  Tm->SetBranchAddress("A1n", &An[0]);
  Tm->SetBranchAddress("A2n", &An[1]);
  Tm->SetBranchAddress("A3n", &An[2]);
  Tm->SetBranchAddress("fn", &fn);
  Tm->SetBranchAddress("Nacc", &Nacc);
  TFile * fs = new TFile(savefile, "RECREATE");
  TTree * Np = new TTree("nuclplus", "nuclplus");
  Np->SetDirectory(fs);
  TTree * Nm = new TTree("nuclminus", "nuclminus");
  Nm->SetDirectory(fs);
  double ErrRel[3];
  double dfn, dpp, dpn;
  Np->Branch("BinNumber", &BinNumber, "BinNumber/D");
  Np->Branch("ErrRel1", &ErrRel[0], "ErrRel1/D");
  Np->Branch("ErrRel2", &ErrRel[1], "ErrRel2/D");
  Np->Branch("ErrRel3", &ErrRel[2], "ErrRel3/D");
  Nm->Branch("BinNumber", &BinNumber, "BinNumber/D");
  Nm->Branch("ErrRel1", &ErrRel[0], "ErrRel1/D");
  Nm->Branch("ErrRel2", &ErrRel[1], "ErrRel2/D");
  Nm->Branch("ErrRel3", &ErrRel[2], "ErrRel3/D");
  std::cout << "Nuclear effect: pi+" << std::endl;
  for (int i = 0; i < Nbinp; i++){
    Tp->GetEntry(i);
    std::cout << "#" << i << " in " << Nbinp << std::endl;
    dfn = sqrt((1.0 - fn) * Nacc * ftime) / Nacc / fn;
    dpn = 0.036 / 0.86;
    for (int j = 0; j < 3; j++){
      dpp = sqrt(1.0/((1.0 - fn) * Nacc * ftime) + pow(0.0065/0.028,2) + 0.01) * ((1.0 - fn) * 0.028 * Ap[j]) / (An[j]);
      ErrRel[j] = sqrt(dfn*dfn + dpn*dpn + dpp*dpp);
    }
    Np->Fill();
  }
  std::cout << "Nuclear effect: pi-" << std::endl;
  for (int i = 0; i < Nbinm; i++){
    Tm->GetEntry(i);
    std::cout << "#" << i << " in " << Nbinm << std::endl;
    dfn = sqrt((1.0 - fn) * Nacc * ftime) / Nacc / fn;
    dpn = 0.036 / 0.86;
    for (int j = 0; j < 3; j++){
      dpp = sqrt(1.0/((1.0 - fn) * Nacc * ftime) + pow(0.0065/0.028, 2) + 0.01) * ((1.0 - fn) * 0.028 * Ap[j]) / (An[j]);
      ErrRel[j] = sqrt(dfn*dfn + dpn*dpn + dpp*dpp);
    }
    Nm->Fill();
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
  ThreetermMatrix(hr, M3inv);
  int Nb = h0->GetNbinsX();
  double average = h0->Integral(1, Nb) / Nb;
  double Nacc = average * Nb;
  double phih;
  double Es[3] = {0.0, 0.0, 0.0};
  double gh;
  for (int i = 1; i <= Nb; i++){
    phih = h0->GetBinCenter(i);
    gh = h0->GetBinContent(i) / average;
    Es[0] = Es[0] + Threeterm1(&M3inv[0], &phih) * gh * 2.0 * M_PI * M_PI / Nacc;
    Es[1] = Es[1] + Threeterm1(&M3inv[3], &phih) * gh * 2.0 * M_PI * M_PI / Nacc;
    Es[2] = Es[2] + Threeterm1(&M3inv[6], &phih) * gh * 2.0 * M_PI * M_PI / Nacc;
  }
  double width = h0->GetBinWidth(1);
  Estat[0] = sqrt(2.0 * Es[0] * width) / _ST;
  Estat[1] = sqrt(2.0 * Es[1] * width) / _ST;
  Estat[2] = sqrt(2.0 * Es[2] * width) / _ST;
  return 0;
}

int Lanalysis::EStatisticsUT3(TH1D * h0, double * Estat){
  int Nb = h0->GetNbinsX();
  double average = h0->Integral(1, Nb) / Nb;
  double Nacc = h0->Integral(1, Nb) *_eff * _lumi * _days * 24.0 * 3600.0 / _simdensity;
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
  Estat[0] = sqrt(2.0*M_PI*M_PI/Nacc * (pow(M3(0,0),2) + pow(M3(0,1),2) + pow(M3(0,2),2)) * M_PI * M_PI);
  Estat[1] = sqrt(2.0*M_PI*M_PI/Nacc * (pow(M3(1,0),2) + pow(M3(1,1),2) + pow(M3(1,2),2)) * M_PI * M_PI);
  Estat[2] = sqrt(2.0*M_PI*M_PI/Nacc * (pow(M3(2,0),2) + pow(M3(2,1),2) + pow(M3(2,2),2)) * M_PI * M_PI);
  return 0;
};

int Lanalysis::EStatisticsUT3(TH2D * hs0, double * Estat){
  int NX = hs0->GetNbinsX();
  int NY = hs0->GetNbinsY();
  double average = hs0->Integral(1, NX, 1, NY) / NX / NY;
  double Nacc = hs0->Integral(1, NX, 1, NY) * _eff * _lumi * _days *24.0 * 3600.0 / _simdensity;
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
  Estat[0] = sqrt(2.0*M_PI*M_PI/Nacc * (pow(M3(0,0),2) + pow(M3(0,1),2) + pow(M3(0,2),2)) * M_PI * M_PI);
  Estat[1] = sqrt(2.0*M_PI*M_PI/Nacc * (pow(M3(1,0),2) + pow(M3(1,1),2) + pow(M3(1,2),2)) * M_PI * M_PI);
  Estat[2] = sqrt(2.0*M_PI*M_PI/Nacc * (pow(M3(2,0),2) + pow(M3(2,1),2) + pow(M3(2,2),2)) * M_PI * M_PI);
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

int Lanalysis::GetResolutionFiles(){
  if (_res_ctrl){
    return 1;
  }
  _file_e = new TFile("/var/phy/project/mepg/tl190/SoLID-SIDIS/SIDIS_electron_resolution_2d.root","r");
  _theta_e = (TH2D *) _file_e->GetObjectChecked("theta_resolution","TH2D");
  _phi_e = (TH2D *) _file_e->GetObjectChecked("phi_resolution","TH2D");
  _p_e = (TH2D *) _file_e->GetObjectChecked("p_resolution","TH2D");
  _z_e = (TH2D *) _file_e->GetObjectChecked("vertexz_resolution","TH2D");
  _file_pi = new TFile("/var/phy/project/mepg/tl190/SoLID-SIDIS/SIDIS_pim_resolution_2d.root","r");
  _theta_pi = (TH2D *) _file_e->GetObjectChecked("theta_resolution","TH2D");
  _phi_pi = (TH2D *) _file_e->GetObjectChecked("phi_resolution","TH2D");
  _p_pi = (TH2D *) _file_e->GetObjectChecked("p_resolution","TH2D");
  _z_pi = (TH2D *) _file_e->GetObjectChecked("vertexz_resolution","TH2D");
  _res_ctrl = 1;
  return 0;
}

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
  Lstructure::CalcVariables(lab, phys);
  double newlab[7], newphys[9], d2[7];
  double accum[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (int n = 0; n < calls; n++){
    for (int i = 0; i < 7; i++){
      newlab[i] = lab[i] + gRandom->Gaus(0.0, dw[i]);
    }
    Lstructure::CalcVariables(newlab, newphys);
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
  ThreetermMatrix(hr, M3inv);
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
      dAUT2 = ThreetermDAUT2H(Asym, phi);
      separate2[0] = Threeterm2(&M3inv[0], phi);
      separate2[1] = Threeterm2(&M3inv[3], phi);
      separate2[2] = Threeterm2(&M3inv[6], phi);
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
      dAUT2 = ThreetermDAUT2H(Asym, phi);
      separate2[0] = Threeterm2(&M3inv[0], phi);
      separate2[1] = Threeterm2(&M3inv[3], phi);
      separate2[2] = Threeterm2(&M3inv[6], phi);
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
  ThreetermMatrix(hr, M3inv);
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
      dAUT2 = ThreetermDAUT2S(Asym, phi);
      separate2[0] = Threeterm2(&M3inv[0], phi);
      separate2[1] = Threeterm2(&M3inv[3], phi);
      separate2[2] = Threeterm2(&M3inv[6], phi);
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
      dAUT2 = ThreetermDAUT2S(Asym, phi);
      separate2[0] = Threeterm2(&M3inv[0], phi);
      separate2[1] = Threeterm2(&M3inv[3], phi);
      separate2[2] = Threeterm2(&M3inv[6], phi);
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

int Lanalysis::RandomCoincidenceSigmaN(const double * AZ, const double * lab, double * xs){
  const double degtorad = M_PI / 180.0;
  double radlen;
  if (lab[0] == 11.0) radlen = 0.04896;
  if (lab[0] == 8.8) radlen = 0.04793;
  double sigma_pip = wiser_sigma(lab[0], lab[4], lab[5] * degtorad, radlen, 0) / 3.89379e+5;
  double sigma_pim = wiser_sigma(lab[0], lab[4], lab[5] * degtorad, radlen, 1) / 3.89379e+5;
  double sigma_e;
  Lstructure::sigmaDISN(AZ, lab, &sigma_e);
  xs[0] = sigma_e * sigma_pip * _lumi;
  xs[1] = sigma_e * sigma_pim * _lumi;
  return 0;
}


  


#endif
