#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>

#include "TROOT.h"
#include "TApplication.h"
#include "Rtypes.h"
#include "TMath.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TEventList.h"

#include "Lanalysis.h"

using namespace std;

int main(int argc, char* argv[]){
  TString sidisbinfile = "sidisbin_8.root";
  TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0229E8Proton";
  TString asymfilename = "sidisasym_8.root";
  int Ntree = 5;
  TFile * fbin = new TFile(sidisbinfile, "r");
  TTree * Fp = (TTree *) fbin->GetObjectChecked("binplus", "TTree");
  TTree * Fm = (TTree *) fbin->GetObjectChecked("binminus", "TTree");
  double Ebeam;
  double xl, xu, xm, ym, Q2l, Q2u, Q2m, zl, zu, zm, Ptl, Ptu, Ptm, hr[2];
  double Asym[3], Asymp[3], AsymNH3[3];
  double fp, fNH3, Nacc;
  Fp->SetBranchAddress("Ebeam", &Ebeam);
  Fp->SetBranchAddress("xl", &xl);
  Fp->SetBranchAddress("xu", &xu);
  Fp->SetBranchAddress("xm", &xm);
  Fp->SetBranchAddress("ym", &ym);
  Fp->SetBranchAddress("Q2l", &Q2l);
  Fp->SetBranchAddress("Q2u", &Q2u);
  Fp->SetBranchAddress("Q2m", &Q2m);
  Fp->SetBranchAddress("zl", &zl);
  Fp->SetBranchAddress("zu", &zu);
  Fp->SetBranchAddress("zm", &zm);
  Fp->SetBranchAddress("Ptl", &Ptl);
  Fp->SetBranchAddress("Ptu", &Ptu);
  Fp->SetBranchAddress("Ptm", &Ptm);
  Fp->SetBranchAddress("ha", &hr[0]);
  Fp->SetBranchAddress("hb", &hr[1]);
  Fp->SetBranchAddress("A1", &Asym[0]);
  Fp->SetBranchAddress("A2", &Asym[1]);
  Fp->SetBranchAddress("A3", &Asym[2]);
  Fp->SetBranchAddress("A1p", &Asymp[0]);
  Fp->SetBranchAddress("A2p", &Asymp[1]);
  Fp->SetBranchAddress("A3p", &Asymp[2]);
  Fp->SetBranchAddress("A1NH3", &AsymNH3[0]);
  Fp->SetBranchAddress("A2NH3", &AsymNH3[1]);
  Fp->SetBranchAddress("A3NH3", &AsymNH3[2]);
  Fp->SetBranchAddress("fp", &fp);
  Fp->SetBranchAddress("fNH3", &fNH3);
  Fp->SetBranchAddress("Nacc", &Nacc);
  Fm->SetBranchAddress("Ebeam", &Ebeam);
  Fm->SetBranchAddress("xl", &xl);
  Fm->SetBranchAddress("xu", &xu);
  Fm->SetBranchAddress("xm", &xm);
  Fm->SetBranchAddress("ym", &ym);
  Fm->SetBranchAddress("Q2l", &Q2l);
  Fm->SetBranchAddress("Q2u", &Q2u);
  Fm->SetBranchAddress("Q2m", &Q2m);
  Fm->SetBranchAddress("zl", &zl);
  Fm->SetBranchAddress("zu", &zu);
  Fm->SetBranchAddress("zm", &zm);
  Fm->SetBranchAddress("Ptl", &Ptl);
  Fm->SetBranchAddress("Ptu", &Ptu);
  Fm->SetBranchAddress("Ptm", &Ptm);
  Fm->SetBranchAddress("ha", &hr[0]);
  Fm->SetBranchAddress("hb", &hr[1]);
  Fm->SetBranchAddress("A1", &Asym[0]);
  Fm->SetBranchAddress("A2", &Asym[1]);
  Fm->SetBranchAddress("A3", &Asym[2]);
  Fm->SetBranchAddress("A1p", &Asymp[0]);
  Fm->SetBranchAddress("A2p", &Asymp[1]);
  Fm->SetBranchAddress("A3p", &Asymp[2]);
  Fm->SetBranchAddress("A1NH3", &AsymNH3[0]);
  Fm->SetBranchAddress("A2NH3", &AsymNH3[1]);
  Fm->SetBranchAddress("A3NH3", &AsymNH3[2]);
  Fm->SetBranchAddress("fp", &fp);
  Fm->SetBranchAddress("fNH3", &fNH3);
  Fm->SetBranchAddress("Nacc", &Nacc);

  TFile * fs = new TFile(asymfilename,"RECREATE");
  TTree * Tp = new TTree("asymplus", "asymplus");
  Tp->SetDirectory(fs);
  TTree * Tm = new TTree("asymminus", "asymminus");
  Tm->SetDirectory(fs);
  double Estat[3], Etarpol[3], EresH[3], EresS[3], Ecoin[3], Edilu[3], Erho[3], Erad[3];
  double Eraw, tmp0;
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
  Tp->Branch("ha", &hr[0], "ha/D");
  Tp->Branch("hb", &hr[1], "hb/D");
  Tp->Branch("A1", &Asym[0], "A1/D");
  Tp->Branch("A2", &Asym[1], "A2/D");
  Tp->Branch("A3", &Asym[2], "A3/D");
  Tp->Branch("A1p", &Asymp[0], "A1p/D");
  Tp->Branch("A2p", &Asymp[1], "A2p/D");
  Tp->Branch("A3p", &Asymp[2], "A3p/D");
  Tp->Branch("A1NH3", &AsymNH3[0], "A1NH3/D");
  Tp->Branch("A2NH3", &AsymNH3[1], "A2NH3/D");
  Tp->Branch("A3NH3", &AsymNH3[2], "A3NH3/D");
  Tp->Branch("fp", &fp, "fp/D");
  Tp->Branch("fNH3", &fNH3, "fNH3/D");
  Tp->Branch("Nacc", &Nacc, "Nacc/D");
  Tp->Branch("E1stat", &Estat[0], "E1stat/D");
  Tp->Branch("E2stat", &Estat[1], "E2stat/D");
  Tp->Branch("E3stat", &Estat[2], "E3stat/D");
  Tp->Branch("Eraw", &Eraw, "Eraw/D");
  Tp->Branch("E1tarpol", &Etarpol[0], "E1tarpol/D");
  Tp->Branch("E2tarpol", &Etarpol[1], "E2tarpol/D");
  Tp->Branch("E3tarpol", &Etarpol[2], "E3tarpol/D");
  Tp->Branch("E1resH", &EresH[0], "E1resH/D");
  Tp->Branch("E2resH", &EresH[1], "E2resH/D");
  Tp->Branch("E3resH", &EresH[2], "E3resH/D");
  Tp->Branch("E1resS", &EresS[0], "E1resS/D");
  Tp->Branch("E2resS", &EresS[1], "E2resS/D");
  Tp->Branch("E3resS", &EresS[2], "E3resS/D");
  Tp->Branch("E1coin", &Ecoin[0], "E1coin/D");
  Tp->Branch("E2coin", &Ecoin[1], "E2coin/D");
  Tp->Branch("E3coin", &Ecoin[2], "E3coin/D");
  Tp->Branch("E1dilu", &Edilu[0], "E1dilu/D");
  Tp->Branch("E2dilu", &Edilu[1], "E2dilu/D");
  Tp->Branch("E3dilu", &Edilu[2], "E3dilu/D");
  Tp->Branch("E1rho", &Erho[0], "E1rho/D");
  Tp->Branch("E2rho", &Erho[1], "E2rho/D");
  Tp->Branch("E3rho", &Erho[2], "E3rho/D");
  Tp->Branch("E1rad", &Erad[0], "E1rad/D");
  Tp->Branch("E2rad", &Erad[1], "E2rad/D");
  Tp->Branch("E3rad", &Erad[2], "E3rad/D");

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
  Tm->Branch("ha", &hr[0], "ha/D");
  Tm->Branch("hb", &hr[1], "hb/D");
  Tm->Branch("A1", &Asym[0], "A1/D");
  Tm->Branch("A2", &Asym[1], "A2/D");
  Tm->Branch("A3", &Asym[2], "A3/D");
  Tm->Branch("A1p", &Asymp[0], "A1p/D");
  Tm->Branch("A2p", &Asymp[1], "A2p/D");
  Tm->Branch("A3p", &Asymp[2], "A3p/D");
  Tm->Branch("A1NH3", &AsymNH3[0], "A1NH3/D");
  Tm->Branch("A2NH3", &AsymNH3[1], "A2NH3/D");
  Tm->Branch("A3NH3", &AsymNH3[2], "A3NH3/D");
  Tm->Branch("fp", &fp, "fp/D");
  Tm->Branch("fNH3", &fNH3, "fNH3/D");
  Tm->Branch("Nacc", &Nacc, "Nacc/D");
  Tm->Branch("E1stat", &Estat[0], "E1stat/D");
  Tm->Branch("E2stat", &Estat[1], "E2stat/D");
  Tm->Branch("E3stat", &Estat[2], "E3stat/D");
  Tm->Branch("Eraw", &Eraw, "Eraw/D");
  Tm->Branch("E1tarpol", &Etarpol[0], "E1tarpol/D");
  Tm->Branch("E2tarpol", &Etarpol[1], "E2tarpol/D");
  Tm->Branch("E3tarpol", &Etarpol[2], "E3tarpol/D");
  Tm->Branch("E1resH", &EresH[0], "E1resH/D");
  Tm->Branch("E2resH", &EresH[1], "E2resH/D");
  Tm->Branch("E3resH", &EresH[2], "E3resH/D");
  Tm->Branch("E1resS", &EresS[0], "E1resS/D");
  Tm->Branch("E2resS", &EresS[1], "E2resS/D");
  Tm->Branch("E3resS", &EresS[2], "E3resS/D");
  Tm->Branch("E1coin", &Ecoin[0], "E1coin/D");
  Tm->Branch("E2coin", &Ecoin[1], "E2coin/D");
  Tm->Branch("E3coin", &Ecoin[2], "E3coin/D");
  Tm->Branch("E1dilu", &Edilu[0], "E1dilu/D");
  Tm->Branch("E2dilu", &Edilu[1], "E2dilu/D");
  Tm->Branch("E3dilu", &Edilu[2], "E3dilu/D");
  Tm->Branch("E1rho", &Erho[0], "E1rho/D");
  Tm->Branch("E2rho", &Erho[1], "E2rho/D");
  Tm->Branch("E3rho", &Erho[2], "E3rho/D");
  Tm->Branch("E1rad", &Erad[0], "E1rad/D");
  Tm->Branch("E2rad", &Erad[1], "E2rad/D");
  Tm->Branch("E3rad", &Erad[2], "E3rad/D");
 
  double acc_ele, acc_pion[2];
  double sigma[2];
  TString nq, nz, npt;
  TString selectedfile;
  Long64_t Nevent;
  double AZ[4] = {3.705, 2.703, 1/3.705, 0.0};
  double kin[7], lab[7], rms[7], sigmacoin;
  const double lumi = 0.1e+10 * pow(0.197327,2);//Polarized Proton
  double Nsim, genvol, days;
  lan->getsimulationinfo(datadir, &Nsim, &genvol, &days);

  int Nbin = Fp->GetEntries();
  //for (int nb = 24; nb < 25; nb++){
  for (int nb = 0; nb < Nbin; nb++){
    Fp->GetEntry(nb);
    std::cout << "pi+ #" << nb << " in " << Nbin << std::endl;
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
      Tdata->Add(Form(datadir+"/out%.2d/Selected/"+selectedfile, it));
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

    TH1D * h0 = new TH1D("h0", "h0", 180, -M_PI, M_PI);
    TH2D * hs0 = new TH2D("hs0", "hs0", 120, -M_PI, M_PI, 60, 0.0, M_PI);
    TH2D * hsH = new TH2D("hsH", "hsH", 120, -M_PI, M_PI, 60, 0.0, M_PI);
    TH2D * hsS = new TH2D("hsS", "hsS", 120, -M_PI, M_PI, 60, 0.0, M_PI);\
    TH1D * hsb = new TH1D("hsb", "hsb", 1000, 0.0, 200.0);
    for (int ie = 0; ie < Nevent; ie++){
      Tdata->GetEntry(ie);
      if (kin[4] < Ptl || kin[4] > Ptu) continue;
      if (kin[0] < xl || kin[0] > xu) continue;
      lst->sigmaUUT(AZ, lab, sigma);
      sigmacoin = lan->RandomCoincidentSigma(lab, 0);
      //lan->CalRMS(lab, rms);
      h0->Fill(abs(kin[5]), sigma[0]*acc_ele*acc_pion[0]);    
      hs0->Fill(kin[5], abs(kin[6]), sigma[0]*acc_ele*acc_pion[0]);
      //hsH->Fill(kin[5], abs(kin[6]), rms[5]*sigma[0]*acc_ele*acc_pion[0]);
      //hsS->Fill(kin[5], abs(kin[6]), rms[6]*sigma[0]*acc_ele*acc_pion[0]);
      hsb->Fill(sigma[0]/sigmacoin, acc_ele*acc_pion[0]);
    }
    h0->Scale(lumi*days*24.0*3600.0*genvol/Nsim/Ntree);
    hs0->Scale(lumi*days*24.0*3600.0*genvol/Nsim/Ntree);
    //Statistical, Target Polarization, Resolution H&S, Coincidence
    //lan->ThreetermStatistics(hr, h0, Estat);
    lan->EStatisticsUT3(hs0, Estat);
    lan->ETargetPolarization(Asym, Etarpol);
    // hsH->Divide(hs0);
    // lan->EResolutionH(hr, Asym, hsH, EresH);
    // hsS->Divide(hs0);
    // lan->EResolutionS(hr, Asym, hsS, EresS);
    lan->ECoincidence(Asym, hsb, Ecoin);
    for (int i = 0; i < 3; i++){
      Estat[i] = Estat[i] / fp;
      Etarpol[i] = Etarpol[i] / fp;
      EresH[i] = EresH[i] / fp;
      EresS[i] = EresS[i] / fp;
      Ecoin[i] = Ecoin[i] / fp;
    }
    //Dilution factor assume 1% uncertainty from C to N
    tmp0 = sqrt(pow(0.01 * (fNH3 - fp), 2) + (1.0 - fNH3) * 0.45 * 55.0 / Nacc) / fp;
    Edilu[0] = abs(tmp0 * Asymp[0]);
    Edilu[1] = abs(tmp0 * Asymp[1]);
    Edilu[2] = abs(tmp0 * Asymp[2]);
    //Diffractive rho production 3% for pi+, 2% for pi-
    Erho[0] = abs(0.03 * Asymp[0]);
    Erho[1] = abs(0.03 * Asymp[1]);
    Erho[2] = abs(0.03 * Asymp[2]);
    //Radiative correction 2%
    Erad[0] = abs(0.02 * Asymp[0]);
    Erad[1] = abs(0.02 * Asymp[1]);
    Erad[2] = abs(0.02 * Asymp[2]);
    //Raw asymmetry
    if (Ebeam == 11.0) Eraw = 7.8E-4 / (0.7 * fp);
    if (Ebeam == 8.8)  Eraw = 1.1E-3 / (0.7 * fp);

    Tp->Fill();

    h0->Delete();
    hs0->Delete();
    hsH->Delete();
    hsS->Delete();
    hsb->Delete();
    Tdata->Delete();
  }

  Nbin = Fm->GetEntries();
  //for (int nb = 3; nb < 3; nb++){
  for (int nb = 0; nb < Nbin; nb++){
    Fm->GetEntry(nb);
    std::cout << "pi- #" << nb << " in " << Nbin << std::endl;
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
      Tdata->Add(Form(datadir+"/out%.2d/Selected/"+selectedfile, it));
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

    TH1D * h0 = new TH1D("h0", "h0", 180, -M_PI, M_PI);
    TH2D * hs0 = new TH2D("hs0", "hs0", 120, -M_PI, M_PI, 60, 0.0, M_PI);
    TH2D * hsH = new TH2D("hsH", "hsH", 120, -M_PI, M_PI, 60, 0.0, M_PI);
    TH2D * hsS = new TH2D("hsS", "hsS", 120, -M_PI, M_PI, 60, 0.0, M_PI);\
    TH1D * hsb = new TH1D("hsb", "hsb", 1000, 0.0, 200.0);
    for (int ie = 0; ie < Nevent; ie++){
      Tdata->GetEntry(ie);
      if (kin[4] < Ptl || kin[4] > Ptu) continue;
      if (kin[0] < xl || kin[0] > xu) continue;
      lst->sigmaUUT(AZ, lab, sigma);
      sigmacoin = lan->RandomCoincidentSigma(lab, 1);
      //lan->CalRMS(lab, rms);
      h0->Fill(abs(kin[5]), sigma[1]*acc_ele*acc_pion[1]);
      hs0->Fill(kin[5], abs(kin[6]), sigma[1]*acc_ele*acc_pion[1]);
      //hsH->Fill(kin[5], abs(kin[6]), rms[5]*sigma[1]*acc_ele*acc_pion[1]);
      //hsS->Fill(kin[5], abs(kin[6]), rms[6]*sigma[1]*acc_ele*acc_pion[1]);
      hsb->Fill(sigma[1]/sigmacoin, acc_ele*acc_pion[1]);
    }
    h0->Scale(lumi*days*24.0*3600.0*genvol/Nsim/Ntree);
    hs0->Scale(lumi*days*24.0*3600.0*genvol/Nsim/Ntree);
    //Statistical, Target Polarization, Resolution H&S, Coincidence
    //lan->ThreetermStatistics(hr, h0, Estat);
    lan->EStatisticsUT3(hs0, Estat);
    lan->ETargetPolarization(Asym, Etarpol);
    // hsH->Divide(hs0);
    // lan->EResolutionH(hr, Asym, hsH, EresH);
    // hsS->Divide(hs0);
    // lan->EResolutionS(hr, Asym, hsS, EresS);
    lan->ECoincidence(Asym, hsb, Ecoin);
    for (int i = 0; i < 3; i++){
      Estat[i] = Estat[i] / fp;
      Etarpol[i] = Etarpol[i] / fp;
      EresH[i] = EresH[i] / fp;
      EresS[i] = EresS[i] / fp;
      Ecoin[i] = Ecoin[i] / fp;
    }
    //Dilution factor assume 1% uncertainty from C to N
    tmp0 = sqrt(pow(0.01 * (fNH3 - fp), 2) + (1.0 - fNH3) * 0.45 * 55.0 / Nacc) / fp;
    Edilu[0] = abs(tmp0 * Asymp[0]);
    Edilu[1] = abs(tmp0 * Asymp[1]);
    Edilu[2] = abs(tmp0 * Asymp[2]);
    //Diffractive rho production 3% for pi+, 2% for pi-
    Erho[0] = abs(0.02 * Asymp[0]);
    Erho[1] = abs(0.02 * Asymp[1]);
    Erho[2] = abs(0.02 * Asymp[2]);
    //Radiative correction 2%
    Erad[0] = abs(0.02 * Asymp[0]);
    Erad[1] = abs(0.02 * Asymp[1]);
    Erad[2] = abs(0.02 * Asymp[2]);
    //Raw asymmetry
    if (Ebeam == 11.0) Eraw = 7.8E-4 / (0.7 * fp);
    if (Ebeam == 8.8)  Eraw = 1.1E-3 / (0.7 * fp);

    Tm->Fill();

    h0->Delete();
    hs0->Delete();
    hsH->Delete();
    hsS->Delete();
    hsb->Delete();
    Tdata->Delete();
  }

  fs->Write();

  return 0;
}
