#ifndef _LFIT_H_
#define _LFIT_H_

#include <iostream>
#include <fstream>
#include <cmath>
#include "TRandom.h"
#include "TRandom3.h"
#include "TEventList.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"
#include "Lstructure.h"

#define _NDATA_ 4000

int _Atype;
int _Ndata;
double _Nucleon[_NDATA_];
double _hadron[_NDATA_];
double _kin[_NDATA_][5];
double _A0[_NDATA_];
double _A[_NDATA_];
double _E[_NDATA_];
double _para[12];

class Lfit{
 protected:
  /* static int _Atype; */
  /* static int _Ndata; */
  /* static double _Nucleon[_NDATA_]; */
  /* static double _hadron[_NDATA_]; */
  /* static double _kin[_NDATA_][5]; */
  /* static double _A0[_NDATA_]; */
  /* static double _A[_NDATA_]; */
  /* static double _E[_NDATA_]; */
  /* static double _para[12]; */
 public:
  Lfit(double a);
  int InitialA1para();
  int InitialA2para();
  int InitialA3para();
  int GetPara(double * para);
  int GetModelPara(double * para);
  int ReadDataIntoTree0();
  int ResetData();
  int PrintDataStatus();
  int GetNeutronPlusData(const int type);
  int GetNeutronMinusData(const int type);
  int GetNeutronData(const int type);
  int GetProtonPlusData(const int type);
  int GetProtonMinusData(const int type);
  int GetProtonData(const int type);
  int GetData(const int type);
  int PrintData(const int n);
  int SmearData();
  static double ChiSquare0(const double * para);
  static double ChiSquare(const double * para);
  static double A1FitFunction(const double * fitpara);
  static double A2FitFunction(const double * fitpara);
  static double A3FitFunction(const double * fitpara);
  int A1Hessian(const double * ds, const double * bestfit);
  int A2Hessian(const double * ds, const double * bestfit);
  int A3Hessian(const double * ds, const double * bestfit);
  int A1Minimizer(const char * minName, const char * algoName);
  int A2Minimizer(const char * minName, const char * algoName);
  int A3Minimizer(const char * minName, const char * algoName);
  int A1MinimizerFree(const char * minName, const char * algoName);
  int A2MinimizerFree(const char * minName, const char * algoName);
  int A3MinimizerFree(const char * minName, const char * algoName);
  int FitSivers(const char * filename, const int loops);
  int FitCollins(const char * filename, const int loops);
  int FitPretzelosity(const char * filename, const int loops);
  int FitSiversFree(const char * filename, const int loops);
  int FitCollinsFree(const char * filename, const int loops);
  int FitPretzelosityFree(const char * filename, const int loops);
  int GetFitf1tM(const double * var, double f1tM[3][6], const double dchi2, const char * file);
  int GetFitf1tkt(const double step, double ktm[3][6], const double dchi2, const char * file);
  int GetFittmdf1t(const double * var, double tmdf1t[3][6], const double dchi2, const char * file);
  int GetFitf1tMcharge(const double Q2, double tc[3][8], const double dchi2, const char * file);
  int GetFith1(const double * var, double h1[3][6], const double dchi2, const char * file);
  int GetFith1charge(const double Q2, double tc[3][6], const double dchi2, const char * file);
  int GetFith1tM(const double * var, double h1t[3][6], const double dchi2, const char * file);
  int GetFith1tMcharge(const double Q2, double tc[3][6], const double dchi2, const char * file);
};

Lfit::Lfit(double a){
  _Atype = 0;
  _Ndata = 0;
}

int Lfit::InitialA1para(){
  //Sivers
  _para[0] = 0.35;//Nu
  _para[1] = -0.90;//Nd
  _para[2] = -0.24;//Ns
  _para[3] = 0.04;//Nubar
  _para[4] = -0.40;//Ndbar
  _para[5] = 1.0;//Nsbar
  _para[6] = 0.73;//au
  _para[7] = 1.08;//ad
  _para[8] = 0.79;//asea
  _para[9] = 3.46;//b
  _para[10] = 0.34;//MS2
  return 0;
}

int Lfit::InitialA2para(){
  //transversity
  _para[0] = 0.46;//Nu
  _para[1] = -1.0;//Nd
  _para[2] = 1.11;//a
  _para[3] = 3.64;//b
  _para[4] = 0.25;//kt2
  //Collins
  _para[5] = 0.49;//Nfav
  _para[6] = -1.0;//Ndis
  _para[7] = 1.06;//c
  _para[8] = 0.07;//d
  _para[9] = 1.5;//MC2
  return 0;
}

int Lfit::InitialA3para(){
  //pretzelosity
  _para[0] = 1.0;//Nu
  _para[1] = -1.0;//Nd
  _para[2] = 2.5;//a
  _para[3] = 2.0;//b
  _para[4] = 0.18;//Mh2
  //Collins
  _para[5] = 0.49;//Nfav
  _para[6] = -1.0;//Ndis
  _para[7] = 1.06;//c
  _para[8] = 0.07;//d
  _para[9] = 1.5;//MC2
  return 0;
}

int Lfit::GetPara(double * para){
  for (int i = 0; i < 12; i++){
    para[i] = _para[i];
  }
  return 0;
}

int Lfit::GetModelPara(double * para){
  if (_Atype == 1){
    for (int i = 0; i < 11; i++){
      para[i] = _para[i];
    }
  }
  else if (_Atype == 2 || _Atype == 3){
    para[0] = _para[0];
    para[1] = _para[1];
    para[2] = _para[2];
    para[3] = _para[2];
    para[4] = _para[3];
    para[5] = _para[3];
    para[6] = _para[4];
    para[7] = _para[5];
    para[8] = _para[6];
    para[9] = _para[7];
    para[10] = _para[7];
    para[11] = _para[8];
    para[12] = _para[8];
    para[13] = _para[9];
  }
  else {
    std::cout << "Lfit::GetModelPara: Unknown asymmetry type!" << std::endl;
    return 1;
  }
  return 0;
}

int Lfit::ReadDataIntoTree0(){
  TString neutronfile11 = "/var/phy/project/mepg/tl190/newSIDIS/sidisasym_11.root";
  TString neutronfile8 = "/var/phy/project/mepg/tl190/newSIDIS/sidisasym_8.root";
  TString protonfile11 = "/var/phy/project/mepg/tl190/protonSIDIS/sidisasym_11.root";
  TString protonfile8 = "/var/phy/project/mepg/tl190/protonSIDIS/sidisasym_8.root";
  double Nucleon, hadron;
  double kin[5];//x, y, z, Q2, Pt
  double A0[3], E0[3];
  double Estat[3], Eraw, Etarpol[3], EresH[3], EresS[3], Ecoin[3], Edilu[3], Erho[3], Erad[3], EPn[3], EPp[3];
  double Nentry;
  TTree * Tr;
  //neutron data
  TFile * fn = new TFile("asymdata_neutron.root","RECREATE");
  TTree * Tn = new TTree("asymdata", "asymdata");
  Tn->SetDirectory(fn);
  Tn->Branch("Nucleon", &Nucleon, "Nucleon/D");
  Tn->Branch("hadron", &hadron, "hadron/D");
  Tn->Branch("x", &kin[0], "x/D");
  Tn->Branch("y", &kin[1], "y/D");
  Tn->Branch("z", &kin[2], "z/D");
  Tn->Branch("Q2", &kin[3], "Q2/D");
  Tn->Branch("Pt", &kin[4], "Pt/D");
  Tn->Branch("A1", &A0[0], "A1/D");
  Tn->Branch("A2", &A0[1], "A2/D");
  Tn->Branch("A3", &A0[2], "A3/D");
  Tn->Branch("E1stat", &Estat[0], "E1stat/D");
  Tn->Branch("E2stat", &Estat[1], "E2stat/D");
  Tn->Branch("E3stat", &Estat[2], "E3stat/D");
  Tn->Branch("E1", &E0[0], "E1/D");
  Tn->Branch("E2", &E0[1], "E2/D");
  Tn->Branch("E3", &E0[2], "E3/D");
  TFile * fn11 = new TFile(neutronfile11,"r");
  TFile * fn8 = new TFile(neutronfile8,"r");
  Tr = (TTree *) fn11->GetObjectChecked("asymplus", "TTree");
  Tr->SetBranchAddress("xm", &kin[0]);
  Tr->SetBranchAddress("ym", &kin[1]);
  Tr->SetBranchAddress("zm", &kin[2]);
  Tr->SetBranchAddress("Q2m", &kin[3]);
  Tr->SetBranchAddress("Ptm", &kin[4]);
  Tr->SetBranchAddress("A1n", &A0[0]);
  Tr->SetBranchAddress("A2n", &A0[1]);
  Tr->SetBranchAddress("A3n", &A0[2]);
  Tr->SetBranchAddress("E1stat", &Estat[0]);
  Tr->SetBranchAddress("E2stat", &Estat[1]);
  Tr->SetBranchAddress("E3stat", &Estat[2]);
  Tr->SetBranchAddress("Eraw", &Eraw);
  Tr->SetBranchAddress("E1tarpol", &Etarpol[0]);
  Tr->SetBranchAddress("E2tarpol", &Etarpol[1]);
  Tr->SetBranchAddress("E3tarpol", &Etarpol[2]);
  Tr->SetBranchAddress("E1resH", &EresH[0]);
  Tr->SetBranchAddress("E2resH", &EresH[1]);
  Tr->SetBranchAddress("E3resH", &EresH[2]);
  Tr->SetBranchAddress("E1resS", &EresS[0]);
  Tr->SetBranchAddress("E2resS", &EresS[1]);
  Tr->SetBranchAddress("E3resS", &EresS[2]);
  Tr->SetBranchAddress("E1coin", &Ecoin[0]);
  Tr->SetBranchAddress("E2coin", &Ecoin[1]);
  Tr->SetBranchAddress("E3coin", &Ecoin[2]);
  Tr->SetBranchAddress("E1dilu", &Edilu[0]);
  Tr->SetBranchAddress("E2dilu", &Edilu[1]);
  Tr->SetBranchAddress("E3dilu", &Edilu[2]);
  Tr->SetBranchAddress("E1rho", &Erho[0]);
  Tr->SetBranchAddress("E2rho", &Erho[1]);
  Tr->SetBranchAddress("E3rho", &Erho[2]);
  Tr->SetBranchAddress("E1rad", &Erad[0]);
  Tr->SetBranchAddress("E2rad", &Erad[1]);
  Tr->SetBranchAddress("E3rad", &Erad[2]);
  Tr->SetBranchAddress("E1Pn", &EPn[0]);
  Tr->SetBranchAddress("E2Pn", &EPn[1]);
  Tr->SetBranchAddress("E3Pn", &EPn[2]);
  Tr->SetBranchAddress("E1Pp", &EPp[0]);
  Tr->SetBranchAddress("E2Pp", &EPp[1]);
  Tr->SetBranchAddress("E3Pp", &EPp[2]);
  Nentry = Tr->GetEntries();
  Nucleon = 0;
  hadron = 0;
  for (int i = 0; i < Nentry; i++){
    Tr->GetEntry(i);
    for (int j = 0; j < 3; j++){
      E0[j] = sqrt(Estat[j]*Estat[j] + Eraw*Eraw + Etarpol[j]*Etarpol[j] + EresH[j]*EresH[j] + EresS[j]*EresS[j] + Ecoin[j]*Ecoin[j] + Edilu[j]*Edilu[j] + Erho[j]*Erho[j] + Erad[j]*Erad[j]+ EPn[j]*EPn[j] + EPp[j]*EPp[j]);
    }
    Tn->Fill();
  }
  Tr = (TTree *) fn8->GetObjectChecked("asymplus", "TTree");
  Tr->SetBranchAddress("xm", &kin[0]);
  Tr->SetBranchAddress("ym", &kin[1]);
  Tr->SetBranchAddress("zm", &kin[2]);
  Tr->SetBranchAddress("Q2m", &kin[3]);
  Tr->SetBranchAddress("Ptm", &kin[4]);
  Tr->SetBranchAddress("A1n", &A0[0]);
  Tr->SetBranchAddress("A2n", &A0[1]);
  Tr->SetBranchAddress("A3n", &A0[2]);
  Tr->SetBranchAddress("E1stat", &Estat[0]);
  Tr->SetBranchAddress("E2stat", &Estat[1]);
  Tr->SetBranchAddress("E3stat", &Estat[2]);
  Tr->SetBranchAddress("Eraw", &Eraw);
  Tr->SetBranchAddress("E1tarpol", &Etarpol[0]);
  Tr->SetBranchAddress("E2tarpol", &Etarpol[1]);
  Tr->SetBranchAddress("E3tarpol", &Etarpol[2]);
  Tr->SetBranchAddress("E1resH", &EresH[0]);
  Tr->SetBranchAddress("E2resH", &EresH[1]);
  Tr->SetBranchAddress("E3resH", &EresH[2]);
  Tr->SetBranchAddress("E1resS", &EresS[0]);
  Tr->SetBranchAddress("E2resS", &EresS[1]);
  Tr->SetBranchAddress("E3resS", &EresS[2]);
  Tr->SetBranchAddress("E1coin", &Ecoin[0]);
  Tr->SetBranchAddress("E2coin", &Ecoin[1]);
  Tr->SetBranchAddress("E3coin", &Ecoin[2]);
  Tr->SetBranchAddress("E1dilu", &Edilu[0]);
  Tr->SetBranchAddress("E2dilu", &Edilu[1]);
  Tr->SetBranchAddress("E3dilu", &Edilu[2]);
  Tr->SetBranchAddress("E1rho", &Erho[0]);
  Tr->SetBranchAddress("E2rho", &Erho[1]);
  Tr->SetBranchAddress("E3rho", &Erho[2]);
  Tr->SetBranchAddress("E1rad", &Erad[0]);
  Tr->SetBranchAddress("E2rad", &Erad[1]);
  Tr->SetBranchAddress("E3rad", &Erad[2]);
  Tr->SetBranchAddress("E1Pn", &EPn[0]);
  Tr->SetBranchAddress("E2Pn", &EPn[1]);
  Tr->SetBranchAddress("E3Pn", &EPn[2]);
  Tr->SetBranchAddress("E1Pp", &EPp[0]);
  Tr->SetBranchAddress("E2Pp", &EPp[1]);
  Tr->SetBranchAddress("E3Pp", &EPp[2]);
  Nentry = Tr->GetEntries();
  Nucleon = 0;
  hadron = 0;
  for (int i = 0; i < Nentry; i++){
    Tr->GetEntry(i);
    for (int j = 0; j < 3; j++){
      E0[j] = sqrt(Estat[j]*Estat[j] + Eraw*Eraw + Etarpol[j]*Etarpol[j] + EresH[j]*EresH[j] + EresS[j]*EresS[j] + Ecoin[j]*Ecoin[j] + Edilu[j]*Edilu[j] + Erho[j]*Erho[j] + Erad[j]*Erad[j]+ EPn[j]*EPn[j] + EPp[j]*EPp[j]);
    }
    Tn->Fill();
  }
  Tr = (TTree *) fn11->GetObjectChecked("asymminus", "TTree");
  Tr->SetBranchAddress("xm", &kin[0]);
  Tr->SetBranchAddress("ym", &kin[1]);
  Tr->SetBranchAddress("zm", &kin[2]);
  Tr->SetBranchAddress("Q2m", &kin[3]);
  Tr->SetBranchAddress("Ptm", &kin[4]);
  Tr->SetBranchAddress("A1n", &A0[0]);
  Tr->SetBranchAddress("A2n", &A0[1]);
  Tr->SetBranchAddress("A3n", &A0[2]);
  Tr->SetBranchAddress("E1stat", &Estat[0]);
  Tr->SetBranchAddress("E2stat", &Estat[1]);
  Tr->SetBranchAddress("E3stat", &Estat[2]);
  Tr->SetBranchAddress("Eraw", &Eraw);
  Tr->SetBranchAddress("E1tarpol", &Etarpol[0]);
  Tr->SetBranchAddress("E2tarpol", &Etarpol[1]);
  Tr->SetBranchAddress("E3tarpol", &Etarpol[2]);
  Tr->SetBranchAddress("E1resH", &EresH[0]);
  Tr->SetBranchAddress("E2resH", &EresH[1]);
  Tr->SetBranchAddress("E3resH", &EresH[2]);
  Tr->SetBranchAddress("E1resS", &EresS[0]);
  Tr->SetBranchAddress("E2resS", &EresS[1]);
  Tr->SetBranchAddress("E3resS", &EresS[2]);
  Tr->SetBranchAddress("E1coin", &Ecoin[0]);
  Tr->SetBranchAddress("E2coin", &Ecoin[1]);
  Tr->SetBranchAddress("E3coin", &Ecoin[2]);
  Tr->SetBranchAddress("E1dilu", &Edilu[0]);
  Tr->SetBranchAddress("E2dilu", &Edilu[1]);
  Tr->SetBranchAddress("E3dilu", &Edilu[2]);
  Tr->SetBranchAddress("E1rho", &Erho[0]);
  Tr->SetBranchAddress("E2rho", &Erho[1]);
  Tr->SetBranchAddress("E3rho", &Erho[2]);
  Tr->SetBranchAddress("E1rad", &Erad[0]);
  Tr->SetBranchAddress("E2rad", &Erad[1]);
  Tr->SetBranchAddress("E3rad", &Erad[2]);
  Tr->SetBranchAddress("E1Pn", &EPn[0]);
  Tr->SetBranchAddress("E2Pn", &EPn[1]);
  Tr->SetBranchAddress("E3Pn", &EPn[2]);
  Tr->SetBranchAddress("E1Pp", &EPp[0]);
  Tr->SetBranchAddress("E2Pp", &EPp[1]);
  Tr->SetBranchAddress("E3Pp", &EPp[2]);
  Nentry = Tr->GetEntries();
  Nucleon = 0;
  hadron = 1;
  for (int i = 0; i < Nentry; i++){
    Tr->GetEntry(i);
    for (int j = 0; j < 3; j++){
      E0[j] = sqrt(Estat[j]*Estat[j] + Eraw*Eraw + Etarpol[j]*Etarpol[j] + EresH[j]*EresH[j] + EresS[j]*EresS[j] + Ecoin[j]*Ecoin[j] + Edilu[j]*Edilu[j] + Erho[j]*Erho[j] + Erad[j]*Erad[j]+ EPn[j]*EPn[j] + EPp[j]*EPp[j]);
    }
    Tn->Fill();
  }
  Tr = (TTree *) fn8->GetObjectChecked("asymminus", "TTree");
  Tr->SetBranchAddress("xm", &kin[0]);
  Tr->SetBranchAddress("ym", &kin[1]);
  Tr->SetBranchAddress("zm", &kin[2]);
  Tr->SetBranchAddress("Q2m", &kin[3]);
  Tr->SetBranchAddress("Ptm", &kin[4]);
  Tr->SetBranchAddress("A1n", &A0[0]);
  Tr->SetBranchAddress("A2n", &A0[1]);
  Tr->SetBranchAddress("A3n", &A0[2]);
  Tr->SetBranchAddress("E1stat", &Estat[0]);
  Tr->SetBranchAddress("E2stat", &Estat[1]);
  Tr->SetBranchAddress("E3stat", &Estat[2]);
  Tr->SetBranchAddress("Eraw", &Eraw);
  Tr->SetBranchAddress("E1tarpol", &Etarpol[0]);
  Tr->SetBranchAddress("E2tarpol", &Etarpol[1]);
  Tr->SetBranchAddress("E3tarpol", &Etarpol[2]);
  Tr->SetBranchAddress("E1resH", &EresH[0]);
  Tr->SetBranchAddress("E2resH", &EresH[1]);
  Tr->SetBranchAddress("E3resH", &EresH[2]);
  Tr->SetBranchAddress("E1resS", &EresS[0]);
  Tr->SetBranchAddress("E2resS", &EresS[1]);
  Tr->SetBranchAddress("E3resS", &EresS[2]);
  Tr->SetBranchAddress("E1coin", &Ecoin[0]);
  Tr->SetBranchAddress("E2coin", &Ecoin[1]);
  Tr->SetBranchAddress("E3coin", &Ecoin[2]);
  Tr->SetBranchAddress("E1dilu", &Edilu[0]);
  Tr->SetBranchAddress("E2dilu", &Edilu[1]);
  Tr->SetBranchAddress("E3dilu", &Edilu[2]);
  Tr->SetBranchAddress("E1rho", &Erho[0]);
  Tr->SetBranchAddress("E2rho", &Erho[1]);
  Tr->SetBranchAddress("E3rho", &Erho[2]);
  Tr->SetBranchAddress("E1rad", &Erad[0]);
  Tr->SetBranchAddress("E2rad", &Erad[1]);
  Tr->SetBranchAddress("E3rad", &Erad[2]);
  Tr->SetBranchAddress("E1Pn", &EPn[0]);
  Tr->SetBranchAddress("E2Pn", &EPn[1]);
  Tr->SetBranchAddress("E3Pn", &EPn[2]);
  Tr->SetBranchAddress("E1Pp", &EPp[0]);
  Tr->SetBranchAddress("E2Pp", &EPp[1]);
  Tr->SetBranchAddress("E3Pp", &EPp[2]);
  Nentry = Tr->GetEntries();
  Nucleon = 0;
  hadron = 1;
  for (int i = 0; i < Nentry; i++){
    Tr->GetEntry(i);
    for (int j = 0; j < 3; j++){
      E0[j] = sqrt(Estat[j]*Estat[j] + Eraw*Eraw + Etarpol[j]*Etarpol[j] + EresH[j]*EresH[j] + EresS[j]*EresS[j] + Ecoin[j]*Ecoin[j] + Edilu[j]*Edilu[j] + Erho[j]*Erho[j] + Erad[j]*Erad[j]+ EPn[j]*EPn[j] + EPp[j]*EPp[j]);
    }
    Tn->Fill();
  }
  fn->Write();
  //proton data
  TFile * fp = new TFile("asymdata_proton.root","RECREATE");
  TTree * Tp = new TTree("asymdata", "asymdata");
  Tp->SetDirectory(fp);
  Tp->Branch("Nucleon", &Nucleon, "Nucleon/D");
  Tp->Branch("hadron", &hadron, "hadron/D");
  Tp->Branch("x", &kin[0], "x/D");
  Tp->Branch("y", &kin[1], "y/D");
  Tp->Branch("z", &kin[2], "z/D");
  Tp->Branch("Q2", &kin[3], "Q2/D");
  Tp->Branch("Pt", &kin[4], "Pt/D");
  Tp->Branch("A1", &A0[0], "A1/D");
  Tp->Branch("A2", &A0[1], "A2/D");
  Tp->Branch("A3", &A0[2], "A3/D");
  Tp->Branch("E1", &E0[0], "E1/D");
  Tp->Branch("E2", &E0[1], "E2/D");
  Tp->Branch("E3", &E0[2], "E3/D");
  TFile * fp11 = new TFile(protonfile11,"r");
  TFile * fp8 = new TFile(protonfile8,"r");
  Tr = (TTree *) fp11->GetObjectChecked("asymplus", "TTree");
  Tr->SetBranchAddress("xm", &kin[0]);
  Tr->SetBranchAddress("ym", &kin[1]);
  Tr->SetBranchAddress("zm", &kin[2]);
  Tr->SetBranchAddress("Q2m", &kin[3]);
  Tr->SetBranchAddress("Ptm", &kin[4]);
  Tr->SetBranchAddress("A1p", &A0[0]);
  Tr->SetBranchAddress("A2p", &A0[1]);
  Tr->SetBranchAddress("A3p", &A0[2]);
  Tr->SetBranchAddress("E1stat", &Estat[0]);
  Tr->SetBranchAddress("E2stat", &Estat[1]);
  Tr->SetBranchAddress("E3stat", &Estat[2]);
  Tr->SetBranchAddress("Eraw", &Eraw);
  Tr->SetBranchAddress("E1tarpol", &Etarpol[0]);
  Tr->SetBranchAddress("E2tarpol", &Etarpol[1]);
  Tr->SetBranchAddress("E3tarpol", &Etarpol[2]);
  Tr->SetBranchAddress("E1resH", &EresH[0]);
  Tr->SetBranchAddress("E2resH", &EresH[1]);
  Tr->SetBranchAddress("E3resH", &EresH[2]);
  Tr->SetBranchAddress("E1resS", &EresS[0]);
  Tr->SetBranchAddress("E2resS", &EresS[1]);
  Tr->SetBranchAddress("E3resS", &EresS[2]);
  Tr->SetBranchAddress("E1coin", &Ecoin[0]);
  Tr->SetBranchAddress("E2coin", &Ecoin[1]);
  Tr->SetBranchAddress("E3coin", &Ecoin[2]);
  Tr->SetBranchAddress("E1dilu", &Edilu[0]);
  Tr->SetBranchAddress("E2dilu", &Edilu[1]);
  Tr->SetBranchAddress("E3dilu", &Edilu[2]);
  Tr->SetBranchAddress("E1rho", &Erho[0]);
  Tr->SetBranchAddress("E2rho", &Erho[1]);
  Tr->SetBranchAddress("E3rho", &Erho[2]);
  Tr->SetBranchAddress("E1rad", &Erad[0]);
  Tr->SetBranchAddress("E2rad", &Erad[1]);
  Tr->SetBranchAddress("E3rad", &Erad[2]);
  Nentry = Tr->GetEntries();
  Nucleon = 1;
  hadron = 0;
  for (int i = 0; i < Nentry; i++){
    Tr->GetEntry(i);
    for (int j = 0; j < 3; j++){
      E0[j] = sqrt(Estat[j]*Estat[j] + Eraw*Eraw + Etarpol[j]*Etarpol[j] + Ecoin[j]*Ecoin[j] + Edilu[j]*Edilu[j] + Erho[j]*Erho[j] + Erad[j]*Erad[j]);
    }
    Tp->Fill();
  }
  Tr = (TTree *) fp8->GetObjectChecked("asymplus", "TTree");
  Tr->SetBranchAddress("xm", &kin[0]);
  Tr->SetBranchAddress("ym", &kin[1]);
  Tr->SetBranchAddress("zm", &kin[2]);
  Tr->SetBranchAddress("Q2m", &kin[3]);
  Tr->SetBranchAddress("Ptm", &kin[4]);
  Tr->SetBranchAddress("A1p", &A0[0]);
  Tr->SetBranchAddress("A2p", &A0[1]);
  Tr->SetBranchAddress("A3p", &A0[2]);
  Tr->SetBranchAddress("E1stat", &Estat[0]);
  Tr->SetBranchAddress("E2stat", &Estat[1]);
  Tr->SetBranchAddress("E3stat", &Estat[2]);
  Tr->SetBranchAddress("Eraw", &Eraw);
  Tr->SetBranchAddress("E1tarpol", &Etarpol[0]);
  Tr->SetBranchAddress("E2tarpol", &Etarpol[1]);
  Tr->SetBranchAddress("E3tarpol", &Etarpol[2]);
  Tr->SetBranchAddress("E1resH", &EresH[0]);
  Tr->SetBranchAddress("E2resH", &EresH[1]);
  Tr->SetBranchAddress("E3resH", &EresH[2]);
  Tr->SetBranchAddress("E1resS", &EresS[0]);
  Tr->SetBranchAddress("E2resS", &EresS[1]);
  Tr->SetBranchAddress("E3resS", &EresS[2]);
  Tr->SetBranchAddress("E1coin", &Ecoin[0]);
  Tr->SetBranchAddress("E2coin", &Ecoin[1]);
  Tr->SetBranchAddress("E3coin", &Ecoin[2]);
  Tr->SetBranchAddress("E1dilu", &Edilu[0]);
  Tr->SetBranchAddress("E2dilu", &Edilu[1]);
  Tr->SetBranchAddress("E3dilu", &Edilu[2]);
  Tr->SetBranchAddress("E1rho", &Erho[0]);
  Tr->SetBranchAddress("E2rho", &Erho[1]);
  Tr->SetBranchAddress("E3rho", &Erho[2]);
  Tr->SetBranchAddress("E1rad", &Erad[0]);
  Tr->SetBranchAddress("E2rad", &Erad[1]);
  Tr->SetBranchAddress("E3rad", &Erad[2]);
  Nentry = Tr->GetEntries();
  Nucleon = 1;
  hadron = 0;
  for (int i = 0; i < Nentry; i++){
    Tr->GetEntry(i);
    for (int j = 0; j < 3; j++){
      E0[j] = sqrt(Estat[j]*Estat[j] + Eraw*Eraw + Etarpol[j]*Etarpol[j] + Ecoin[j]*Ecoin[j] + Edilu[j]*Edilu[j] + Erho[j]*Erho[j] + Erad[j]*Erad[j]);
    }
    Tp->Fill();
  }
  Tr = (TTree *) fp11->GetObjectChecked("asymminus", "TTree");
  Tr->SetBranchAddress("xm", &kin[0]);
  Tr->SetBranchAddress("ym", &kin[1]);
  Tr->SetBranchAddress("zm", &kin[2]);
  Tr->SetBranchAddress("Q2m", &kin[3]);
  Tr->SetBranchAddress("Ptm", &kin[4]);
  Tr->SetBranchAddress("A1p", &A0[0]);
  Tr->SetBranchAddress("A2p", &A0[1]);
  Tr->SetBranchAddress("A3p", &A0[2]);
  Tr->SetBranchAddress("E1stat", &Estat[0]);
  Tr->SetBranchAddress("E2stat", &Estat[1]);
  Tr->SetBranchAddress("E3stat", &Estat[2]);
  Tr->SetBranchAddress("Eraw", &Eraw);
  Tr->SetBranchAddress("E1tarpol", &Etarpol[0]);
  Tr->SetBranchAddress("E2tarpol", &Etarpol[1]);
  Tr->SetBranchAddress("E3tarpol", &Etarpol[2]);
  Tr->SetBranchAddress("E1resH", &EresH[0]);
  Tr->SetBranchAddress("E2resH", &EresH[1]);
  Tr->SetBranchAddress("E3resH", &EresH[2]);
  Tr->SetBranchAddress("E1resS", &EresS[0]);
  Tr->SetBranchAddress("E2resS", &EresS[1]);
  Tr->SetBranchAddress("E3resS", &EresS[2]);
  Tr->SetBranchAddress("E1coin", &Ecoin[0]);
  Tr->SetBranchAddress("E2coin", &Ecoin[1]);
  Tr->SetBranchAddress("E3coin", &Ecoin[2]);
  Tr->SetBranchAddress("E1dilu", &Edilu[0]);
  Tr->SetBranchAddress("E2dilu", &Edilu[1]);
  Tr->SetBranchAddress("E3dilu", &Edilu[2]);
  Tr->SetBranchAddress("E1rho", &Erho[0]);
  Tr->SetBranchAddress("E2rho", &Erho[1]);
  Tr->SetBranchAddress("E3rho", &Erho[2]);
  Tr->SetBranchAddress("E1rad", &Erad[0]);
  Tr->SetBranchAddress("E2rad", &Erad[1]);
  Tr->SetBranchAddress("E3rad", &Erad[2]);
  Nentry = Tr->GetEntries();
  Nucleon = 1;
  hadron = 1;
  for (int i = 0; i < Nentry; i++){
    Tr->GetEntry(i);
    for (int j = 0; j < 3; j++){
      E0[j] = sqrt(Estat[j]*Estat[j] + Eraw*Eraw + Etarpol[j]*Etarpol[j] + Ecoin[j]*Ecoin[j] + Edilu[j]*Edilu[j] + Erho[j]*Erho[j] + Erad[j]*Erad[j]);
    }
    Tp->Fill();
  }
  Tr = (TTree *) fp8->GetObjectChecked("asymminus", "TTree");
  Tr->SetBranchAddress("xm", &kin[0]);
  Tr->SetBranchAddress("ym", &kin[1]);
  Tr->SetBranchAddress("zm", &kin[2]);
  Tr->SetBranchAddress("Q2m", &kin[3]);
  Tr->SetBranchAddress("Ptm", &kin[4]);
  Tr->SetBranchAddress("A1p", &A0[0]);
  Tr->SetBranchAddress("A2p", &A0[1]);
  Tr->SetBranchAddress("A3p", &A0[2]);
  Tr->SetBranchAddress("E1stat", &Estat[0]);
  Tr->SetBranchAddress("E2stat", &Estat[1]);
  Tr->SetBranchAddress("E3stat", &Estat[2]);
  Tr->SetBranchAddress("Eraw", &Eraw);
  Tr->SetBranchAddress("E1tarpol", &Etarpol[0]);
  Tr->SetBranchAddress("E2tarpol", &Etarpol[1]);
  Tr->SetBranchAddress("E3tarpol", &Etarpol[2]);
  Tr->SetBranchAddress("E1resH", &EresH[0]);
  Tr->SetBranchAddress("E2resH", &EresH[1]);
  Tr->SetBranchAddress("E3resH", &EresH[2]);
  Tr->SetBranchAddress("E1resS", &EresS[0]);
  Tr->SetBranchAddress("E2resS", &EresS[1]);
  Tr->SetBranchAddress("E3resS", &EresS[2]);
  Tr->SetBranchAddress("E1coin", &Ecoin[0]);
  Tr->SetBranchAddress("E2coin", &Ecoin[1]);
  Tr->SetBranchAddress("E3coin", &Ecoin[2]);
  Tr->SetBranchAddress("E1dilu", &Edilu[0]);
  Tr->SetBranchAddress("E2dilu", &Edilu[1]);
  Tr->SetBranchAddress("E3dilu", &Edilu[2]);
  Tr->SetBranchAddress("E1rho", &Erho[0]);
  Tr->SetBranchAddress("E2rho", &Erho[1]);
  Tr->SetBranchAddress("E3rho", &Erho[2]);
  Tr->SetBranchAddress("E1rad", &Erad[0]);
  Tr->SetBranchAddress("E2rad", &Erad[1]);
  Tr->SetBranchAddress("E3rad", &Erad[2]);
  Nentry = Tr->GetEntries();
  Nucleon = 1;
  hadron = 1;
  for (int i = 0; i < Nentry; i++){
    Tr->GetEntry(i);
    for (int j = 0; j < 3; j++){
      E0[j] = sqrt(Estat[j]*Estat[j] + Eraw*Eraw + Etarpol[j]*Etarpol[j] + Ecoin[j]*Ecoin[j] + Edilu[j]*Edilu[j] + Erho[j]*Erho[j] + Erad[j]*Erad[j]);
    }
    Tp->Fill();
  }
  fp->Write();
  return 0;
}

int Lfit::ResetData(){
  _Atype = 0;
  _Ndata = 0;
  return 0;
}

int Lfit::PrintDataStatus(){
  std::cout << "Data status: " << std::endl;
  std::cout << "Number of Data: " << _Ndata << "    Asymmetry type: " << _Atype << std::endl;
  return 0;
}

int Lfit::GetNeutronPlusData(const int type){
  if (type != 1 && type != 2 && type != 3){
    std::cout << "Lfit::GetNeutronData: No such asymmetry type!" << std::endl;
    return 1;
  }
  _Atype = type;
  TFile * fa = new TFile("asymdata_neutron.root", "r");
  TTree * Ta = (TTree *) fa->GetObjectChecked("asymdata", "TTree");
  double Nentry = Ta->GetEntries();
  double Nucleon, hadron;
  double kin[5];
  double A0, E0;
  Ta->SetBranchAddress("Nucleon", &Nucleon);
  Ta->SetBranchAddress("hadron", &hadron);
  Ta->SetBranchAddress("x", &kin[0]);
  Ta->SetBranchAddress("y", &kin[1]);
  Ta->SetBranchAddress("z", &kin[2]);
  Ta->SetBranchAddress("Q2", &kin[3]);
  Ta->SetBranchAddress("Pt", &kin[4]);
  if (type == 1){
    Ta->SetBranchAddress("A1", &A0);
    Ta->SetBranchAddress("E1", &E0);
  }
  else if (type == 2){
    Ta->SetBranchAddress("A2", &A0);
    Ta->SetBranchAddress("E2", &E0);
  }
  else if (type == 3){
    Ta->SetBranchAddress("A3", &A0);
    Ta->SetBranchAddress("E3", &E0);
  }
  for (int i = 0; i < Nentry; i++){
    Ta->GetEntry(i);
    if (E0 < 10 && hadron == 0) {
      _Nucleon[_Ndata] = Nucleon;
      _hadron[_Ndata] = hadron;
      for (int j = 0; j < 5; j++){
	_kin[_Ndata][j] = kin[j];
      }
      _A0[_Ndata] = A0;
      _A[_Ndata] = A0;
      _E[_Ndata] = E0;
      _Ndata = _Ndata + 1;
    }
  }
  return 0;
}

int Lfit::GetNeutronMinusData(const int type){
  if (type != 1 && type != 2 && type != 3){
    std::cout << "Lfit::GetNeutronData: No such asymmetry type!" << std::endl;
    return 1;
  }
  _Atype = type;
  TFile * fa = new TFile("asymdata_neutron.root", "r");
  TTree * Ta = (TTree *) fa->GetObjectChecked("asymdata", "TTree");
  double Nentry = Ta->GetEntries();
  double Nucleon, hadron;
  double kin[5];
  double A0, E0;
  Ta->SetBranchAddress("Nucleon", &Nucleon);
  Ta->SetBranchAddress("hadron", &hadron);
  Ta->SetBranchAddress("x", &kin[0]);
  Ta->SetBranchAddress("y", &kin[1]);
  Ta->SetBranchAddress("z", &kin[2]);
  Ta->SetBranchAddress("Q2", &kin[3]);
  Ta->SetBranchAddress("Pt", &kin[4]);
  if (type == 1){
    Ta->SetBranchAddress("A1", &A0);
    Ta->SetBranchAddress("E1", &E0);
  }
  else if (type == 2){
    Ta->SetBranchAddress("A2", &A0);
    Ta->SetBranchAddress("E2", &E0);
  }
  else if (type == 3){
    Ta->SetBranchAddress("A3", &A0);
    Ta->SetBranchAddress("E3", &E0);
  }
  for (int i = 0; i < Nentry; i++){
    Ta->GetEntry(i);
    if (E0 < 10 && hadron == 1) {
      _Nucleon[_Ndata] = Nucleon;
      _hadron[_Ndata] = hadron;
      for (int j = 0; j < 5; j++){
	_kin[_Ndata][j] = kin[j];
      }
      _A0[_Ndata] = A0;
      _A[_Ndata] = A0;
      _E[_Ndata] = E0;
      _Ndata = _Ndata + 1;
    }
  }
  return 0;
}

int Lfit::GetNeutronData(const int type){
  if (type != 1 && type != 2 && type != 3){
    std::cout << "Lfit::GetNeutronData: No such asymmetry type!" << std::endl;
    return 1;
  }
  _Atype = type;
  TFile * fa = new TFile("asymdata_neutron.root", "r");
  TTree * Ta = (TTree *) fa->GetObjectChecked("asymdata", "TTree");
  double Nentry = Ta->GetEntries();
  double Nucleon, hadron;
  double kin[5];
  double A0, E0;
  Ta->SetBranchAddress("Nucleon", &Nucleon);
  Ta->SetBranchAddress("hadron", &hadron);
  Ta->SetBranchAddress("x", &kin[0]);
  Ta->SetBranchAddress("y", &kin[1]);
  Ta->SetBranchAddress("z", &kin[2]);
  Ta->SetBranchAddress("Q2", &kin[3]);
  Ta->SetBranchAddress("Pt", &kin[4]);
  if (type == 1){
    Ta->SetBranchAddress("A1", &A0);
    Ta->SetBranchAddress("E1", &E0);
  }
  else if (type == 2){
    Ta->SetBranchAddress("A2", &A0);
    Ta->SetBranchAddress("E2", &E0);
  }
  else if (type == 3){
    Ta->SetBranchAddress("A3", &A0);
    Ta->SetBranchAddress("E3", &E0);
  }
  for (int i = 0; i < Nentry; i++){
    Ta->GetEntry(i);
    if (E0 < 10) {
      _Nucleon[_Ndata] = Nucleon;
      _hadron[_Ndata] = hadron;
      for (int j = 0; j < 5; j++){
	_kin[_Ndata][j] = kin[j];
      }
      _A0[_Ndata] = A0;
      _A[_Ndata] = A0;
      _E[_Ndata] = E0;
      _Ndata = _Ndata + 1;
    }
  }
  return 0;
}

int Lfit::GetProtonPlusData(const int type){
  if (type != 1 && type != 2 && type != 3){
    std::cout << "Lfit::GetNeutronData: No such asymmetry type!" << std::endl;
    return 1;
  }
  _Atype = type;
  TFile * fa = new TFile("asymdata_proton.root", "r");
  TTree * Ta = (TTree *) fa->GetObjectChecked("asymdata", "TTree");
  double Nentry = Ta->GetEntries();
  double Nucleon, hadron;
  double kin[5];
  double A0, E0;
  Ta->SetBranchAddress("Nucleon", &Nucleon);
  Ta->SetBranchAddress("hadron", &hadron);
  Ta->SetBranchAddress("x", &kin[0]);
  Ta->SetBranchAddress("y", &kin[1]);
  Ta->SetBranchAddress("z", &kin[2]);
  Ta->SetBranchAddress("Q2", &kin[3]);
  Ta->SetBranchAddress("Pt", &kin[4]);
  if (type == 1){
    Ta->SetBranchAddress("A1", &A0);
    Ta->SetBranchAddress("E1", &E0);
  }
  else if (type == 2){
    Ta->SetBranchAddress("A2", &A0);
    Ta->SetBranchAddress("E2", &E0);
  }
  else if (type == 3){
    Ta->SetBranchAddress("A3", &A0);
    Ta->SetBranchAddress("E3", &E0);
  }
  for (int i = 0; i < Nentry; i++){
    Ta->GetEntry(i);
    if (E0 < 10 && hadron == 0){
      _Nucleon[_Ndata] = Nucleon;
      _hadron[_Ndata] = hadron;
      for (int j = 0; j < 5; j++){
	_kin[_Ndata][j] = kin[j];
      }
      _A0[_Ndata] = A0;
      _A[_Ndata] = A0;
      _E[_Ndata] = E0;
      _Ndata = _Ndata + 1;
    }
  }
  return 0;
}

int Lfit::GetProtonMinusData(const int type){
  if (type != 1 && type != 2 && type != 3){
    std::cout << "Lfit::GetNeutronData: No such asymmetry type!" << std::endl;
    return 1;
  }
  _Atype = type;
  TFile * fa = new TFile("asymdata_proton.root", "r");
  TTree * Ta = (TTree *) fa->GetObjectChecked("asymdata", "TTree");
  double Nentry = Ta->GetEntries();
  double Nucleon, hadron;
  double kin[5];
  double A0, E0;
  Ta->SetBranchAddress("Nucleon", &Nucleon);
  Ta->SetBranchAddress("hadron", &hadron);
  Ta->SetBranchAddress("x", &kin[0]);
  Ta->SetBranchAddress("y", &kin[1]);
  Ta->SetBranchAddress("z", &kin[2]);
  Ta->SetBranchAddress("Q2", &kin[3]);
  Ta->SetBranchAddress("Pt", &kin[4]);
  if (type == 1){
    Ta->SetBranchAddress("A1", &A0);
    Ta->SetBranchAddress("E1", &E0);
  }
  else if (type == 2){
    Ta->SetBranchAddress("A2", &A0);
    Ta->SetBranchAddress("E2", &E0);
  }
  else if (type == 3){
    Ta->SetBranchAddress("A3", &A0);
    Ta->SetBranchAddress("E3", &E0);
  }
  for (int i = 0; i < Nentry; i++){
    Ta->GetEntry(i);
    if (E0 < 10 && hadron == 1){
      _Nucleon[_Ndata] = Nucleon;
      _hadron[_Ndata] = hadron;
      for (int j = 0; j < 5; j++){
	_kin[_Ndata][j] = kin[j];
      }
      _A0[_Ndata] = A0;
      _A[_Ndata] = A0;
      _E[_Ndata] = E0;
      _Ndata = _Ndata + 1;
    }
  }
  return 0;
}

int Lfit::GetProtonData(const int type){
  if (type != 1 && type != 2 && type != 3){
    std::cout << "Lfit::GetNeutronData: No such asymmetry type!" << std::endl;
    return 1;
  }
  _Atype = type;
  TFile * fa = new TFile("asymdata_proton.root", "r");
  TTree * Ta = (TTree *) fa->GetObjectChecked("asymdata", "TTree");
  double Nentry = Ta->GetEntries();
  double Nucleon, hadron;
  double kin[5];
  double A0, E0;
  Ta->SetBranchAddress("Nucleon", &Nucleon);
  Ta->SetBranchAddress("hadron", &hadron);
  Ta->SetBranchAddress("x", &kin[0]);
  Ta->SetBranchAddress("y", &kin[1]);
  Ta->SetBranchAddress("z", &kin[2]);
  Ta->SetBranchAddress("Q2", &kin[3]);
  Ta->SetBranchAddress("Pt", &kin[4]);
  if (type == 1){
    Ta->SetBranchAddress("A1", &A0);
    Ta->SetBranchAddress("E1", &E0);
  }
  else if (type == 2){
    Ta->SetBranchAddress("A2", &A0);
    Ta->SetBranchAddress("E2", &E0);
  }
  else if (type == 3){
    Ta->SetBranchAddress("A3", &A0);
    Ta->SetBranchAddress("E3", &E0);
  }
  for (int i = 0; i < Nentry; i++){
    Ta->GetEntry(i);
    if (E0 < 10){
      _Nucleon[_Ndata] = Nucleon;
      _hadron[_Ndata] = hadron;
      for (int j = 0; j < 5; j++){
	_kin[_Ndata][j] = kin[j];
      }
      _A0[_Ndata] = A0;
      _A[_Ndata] = A0;
      _E[_Ndata] = E0;
      _Ndata = _Ndata + 1;
    }
  }
  return 0;
}

int Lfit::GetData(const int type){
  GetNeutronData(type);
  GetProtonData(type);
  return 0;
}

int Lfit::PrintData(const int n){
  if (n >= _Ndata){
    std::cout << "Lfit::PrintData: Overflow data number!" << std::endl;
    return 1;
  }
  std::cout << "Asymmetry type: " << _Atype << std::endl;
  std::cout << "Entry:   " << n << std::endl;
  std::cout << "Nucleon: " << _Nucleon[n] << std::endl;
  std::cout << "hadron:  " << _hadron[n] << std::endl;
  std::cout << "x:   " << _kin[n][0] << std::endl;
  std::cout << "y:   " << _kin[n][1] << std::endl;
  std::cout << "z:   " << _kin[n][2] << std::endl;
  std::cout << "Q2:  " << _kin[n][3] << std::endl;
  std::cout << "Pt:  " << _kin[n][4] << std::endl;
  std::cout << "A0:  " << _A0[n] << std::endl;
  std::cout << "A:   " << _A[n] << std::endl;
  std::cout << "E:   " << _E[n] << std::endl;
  return 0;
}

int Lfit::SmearData(){
  TRandom3 ran;
  ran.SetSeed(0);
  for (int i = 0; i < _Ndata; i++){
    _A[i] = _A0[i] + ran.Gaus(0.0, _E[i]);
  }
  return 0;
}

double Lfit::ChiSquare0(const double * para){
  double Asym[2] = {0, 0};
  double sum = 0;
  for (int i = 0; i < _Ndata; i++){
    if (_Atype == 1){
      if (_Nucleon[i] == 0){
	lst->AsinHmSn(_kin[i], Asym, para);
      }
      else if (_Nucleon[i] == 1){
	lst->AsinHmSp(_kin[i], Asym, para);
      }
      sum = sum + pow((Asym[(int)_hadron[i]] - _A0[i]) / _E[i], 2);
    }
    else if (_Atype == 2){
      if (_Nucleon[i] == 0){
	lst->AsinHpSn(_kin[i], Asym, &para[0], &para[7]);
      }
      else if (_Nucleon[i] == 1){
	lst->AsinHpSp(_kin[i], Asym, &para[0], &para[7]);
      }
      sum = sum + pow((Asym[(int)_hadron[i]] - _A0[i]) / _E[i], 2);
    }
    else if (_Atype == 3){
      if (_Nucleon[i] == 0){
	lst->Asin3HmSn(_kin[i], Asym, &para[0], &para[7]);
      }
      else if (_Nucleon[i] == 1){
	lst->Asin3HmSp(_kin[i], Asym, &para[0], &para[7]);
      }
      sum = sum + pow((Asym[(int)_hadron[i]] - _A0[i]) / _E[i], 2);
    }
    else {
      sum = -1;
    }
  }
  return sum;
}

double Lfit::ChiSquare(const double * para){
  double Asym[2] = {0, 0};
  double sum = 0;
  for (int i = 0; i < _Ndata; i++){
    if (_Atype == 1){
      if (_Nucleon[i] == 0){
	lst->AsinHmSn(_kin[i], Asym, para);
      }
      else if (_Nucleon[i] == 1){
	lst->AsinHmSp(_kin[i], Asym, para);
      }
      sum = sum + pow((Asym[(int)_hadron[i]] - _A[i]) / _E[i], 2);
    }
    else if (_Atype == 2){
      if (_Nucleon[i] == 0){
	lst->AsinHpSn(_kin[i], Asym, &para[0], &para[7]);
      }
      else if (_Nucleon[i] == 1){
	lst->AsinHpSp(_kin[i], Asym, &para[0], &para[7]);
      }
      sum = sum + pow((Asym[(int)_hadron[i]] - _A[i]) / _E[i], 2);
    }
    else if (_Atype == 3){
      if (_Nucleon[i] == 0){
	lst->Asin3HmSn(_kin[i], Asym, &para[0], &para[7]);
      }
      else if (_Nucleon[i] == 1){
	lst->Asin3HmSp(_kin[i], Asym, &para[0], &para[7]);
      }
      sum = sum + pow((Asym[(int)_hadron[i]] - _A[i]) / _E[i], 2);
    }
    else {
      sum = -1;
    }
  }
  return sum;
}

double Lfit::A1FitFunction(const double * fitpara){
  if (_Atype != 1){
    std::cout << "Lfit::A1FitFunction: Wrong asymmetry type!" << std::endl;
    return -1;
  }
  double para[11];
  for (int i = 0; i < 11; i++){
    para[i] = fitpara[i];
  }
  return ChiSquare(para);
}

double Lfit::A2FitFunction(const double * fitpara){
  if (_Atype != 2){
    std::cout << "Lfit::A2FitFunction: Wrong asymmetry type!" << std::endl;
    return -1;
  }
  double para[14];
  para[0] = fitpara[0];
  para[1] = fitpara[1];
  para[2] = fitpara[2];
  para[3] = fitpara[2];
  para[4] = fitpara[3];
  para[5] = fitpara[3];
  para[6] = fitpara[4];
  para[7] = fitpara[5];
  para[8] = fitpara[6];
  para[9] = fitpara[7];
  para[10] = fitpara[7];
  para[11] = fitpara[8];
  para[12] = fitpara[8];
  para[13] = fitpara[9];
  return ChiSquare(para);
}
  
double Lfit::A3FitFunction(const double * fitpara){
  if (_Atype != 3){
    std::cout << "Lfit::A3FitFunction: Wrong asymmetry type!" << std::endl;
    return -1;
  }
  double para[14];
  para[0] = fitpara[0];
  para[1] = fitpara[1];
  para[2] = fitpara[2];
  para[3] = fitpara[2];
  para[4] = fitpara[3];
  para[5] = fitpara[3];
  para[6] = fitpara[4];
  para[7] = fitpara[5];
  para[8] = fitpara[6];
  para[9] = fitpara[7];
  para[10] = fitpara[7];
  para[11] = fitpara[8];
  para[12] = fitpara[8];
  para[13] = fitpara[9];
  return ChiSquare(para);
}

int Lfit::A2Hessian(const double * ds, const double * bestfit = _para){
  if (_Atype != 2){
    std::cout << "Lfit::A2Hessian: Wrong asymmetry type!" << std::endl;
    return 1;
  }
  double da = ds[0];
  double db = ds[1];
  double chi2min = A2FitFunction(bestfit);
  double chi2[4];
  std::cout << "Best fit Chi2: " << chi2min << std::endl;
  int nlist[9] = {0, 1, 2, 3, 5, 6, 7, 8, 9};
  double para[10];
  TMatrixD hessian(9,9);
  //first derivatives check
  /*std::cout << "First derivatives check: " << std::endl;
  for (int i = 0; i < 9; i++){
    for (int s = 0; s < 10; s++) para[s] = bestfit[s];
    para[nlist[i]] = para[nlist[i]] + da;
    chi2[0] = A2FitFunction(para);//a+da
    para[nlist[i]] = para[nlist[i]] - 2.0 * da;
    chi2[1] = A2FitFunction(para);//a-da
    std::cout << nlist[i] << ": " << (chi2[1] - chi2[0]) / (2.0 * da) << std::endl;
    }*/
  for (int i = 0; i < 9; i++){
    for (int j = 0; j < 9; j++){
      for (int s = 0; s < 10; s++) para[s] = bestfit[s];
      para[nlist[i]] = para[nlist[i]] + da;
      para[nlist[j]] = para[nlist[j]] + db;
      chi2[0] = A2FitFunction(para);//a+da, b+db
      para[nlist[j]] = para[nlist[j]] - 2.0 * db;
      chi2[1] = A2FitFunction(para);//a+da, b-db
      para[nlist[i]] = para[nlist[i]] - 2.0 * da;
      chi2[3] = A2FitFunction(para);//a-da, b-db
      para[nlist[j]] = para[nlist[j]] + 2.0 * db;
      chi2[2] = A2FitFunction(para);//a-da, b+db
      if (i == j) hessian(i, i) = (chi2[0] + chi2[3] - 2.0 * chi2min) / (2.0 * pow(da+db, 2));
      else hessian(i, j) = (chi2[0] + chi2[3] - chi2[1] - chi2[2]) / (8.0 * da * db);
    }
  }
  hessian.Print();
  return 0;
}


int Lfit::A1Minimizer(const char * minName = "Minuit", const char * algoName = ""){
  if (_Atype != 1){
    std::cout << "Lfit::A1Minimizer: Wrong asymmetry type!" << std::endl;
    return 1;
  }
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(100000); // for Minuit
  min->SetMaxIterations(10000);  // for GSL 
  min->SetTolerance(0.1);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&Lfit::A1FitFunction, 11);//Fitting function
  double step[11], variable[11];
  InitialA1para();
  for (int i = 0; i < 11; i++){
    step[i] = 0.001;
    variable[i] = _para[i];
  }
  min->SetFunction(f);
  //Set variables to be minimized
  min->SetLimitedVariable(0, "Nu", variable[0], step[0], -1, 1);
  min->SetLimitedVariable(1,"Nd",variable[1],step[1], -1, 1);
  min->SetFixedVariable(2,"Ns",variable[2]);
  min->SetLimitedVariable(3,"Nubar",variable[3],step[3], -1, 1);
  //min->SetFixedVariable(3, "Nubar", variable[3]);
  min->SetLimitedVariable(4,"Ndbar",variable[4],step[4], -1, 1);
  //min->SetFixedVariable(4, "Ndbar", variable[4]);
  min->SetFixedVariable(5,"Nsbar",variable[5]);
  min->SetVariable(6,"au",variable[6],step[6]);
  min->SetVariable(7,"ad",variable[7],step[7]);
  min->SetVariable(8,"asea",variable[8],step[8]);
  min->SetVariable(9,"b",variable[9],step[9]);
  min->SetVariable(10,"Ms2",variable[10],step[10]);
  // do the minimization
  min->Minimize();
  // get fitting result
  const double * xs = min->X();
  for (int i = 0; i < 11; i++){
    _para[i] = xs[i];
  }
  double para[11];
  GetModelPara(para);
  _para[11] = ChiSquare0(para);
  return 0;
}

int Lfit::A2Minimizer(const char * minName = "Minuit", const char * algoName = ""){
  if (_Atype != 2){
    std::cout << "Lfit::A2Minimizer: Wrong asymmetry type!" << std::endl;
    return 1;
  }
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(100000); // for Minuit
  min->SetMaxIterations(10000);  // for GSL 
  min->SetTolerance(0.1);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&Lfit::A2FitFunction, 10);//Fitting function
  double step[10], variable[10];
  InitialA2para();
  for (int i = 0; i < 10; i++){
    step[i] = 0.001;
    variable[i] = _para[i];
  }
  min->SetFunction(f);
  //Set variables to be minimized
  min->SetLimitedVariable(0, "Nu", variable[0], step[0], -1, 1);
  min->SetLimitedVariable(1,"Nd",variable[1],step[1], -1, 1);
  min->SetVariable(2,"a",variable[2],step[2]);
  min->SetVariable(3,"b",variable[3],step[3]);
  min->SetFixedVariable(4,"kt2",variable[4]);
  min->SetLimitedVariable(5,"Nfav",variable[5], step[5], -1, 1);
  min->SetLimitedVariable(6,"Ndis",variable[6],step[6], -1, 1);
  min->SetVariable(7,"c",variable[7],step[7]);
  min->SetVariable(8,"d",variable[8],step[8]);
  min->SetVariable(9,"Mc2",variable[9],step[9]);
  // do the minimization
  min->Minimize();
  // get fitting result
  const double * xs = min->X();
  for (int i = 0; i < 10; i++){
    _para[i] = xs[i];
  }
  double para[11];
  GetModelPara(para);
  _para[11] = ChiSquare0(para);
  return 0;
}

int Lfit::A3Minimizer(const char * minName = "Minuit", const char * algoName = ""){
  if (_Atype != 3){
    std::cout << "Lfit::A3Minimizer: Wrong asymmetry type!" << std::endl;
    return 1;
  }
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(100000); // for Minuit
  min->SetMaxIterations(10000);  // for GSL 
  min->SetTolerance(0.1);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&Lfit::A3FitFunction, 10);//Fitting function
  double step[10], variable[10];
  InitialA3para();
  for (int i = 0; i < 10; i++){
    step[i] = 0.001;
    variable[i] = _para[i];
  }
  min->SetFunction(f);
  //Set variables to be minimized
  min->SetLimitedVariable(0, "Nu", variable[0], step[0], -1, 1);
  min->SetLimitedVariable(1,"Nd",variable[1],step[1], -1, 1);
  min->SetVariable(2,"a",variable[2],step[2]);
  min->SetFixedVariable(3,"b",variable[3]);
  min->SetVariable(4,"Mt2",variable[4], step[4]);
  min->SetFixedVariable(5,"Nfav",variable[5]);
  min->SetFixedVariable(6,"Ndis",variable[6]);
  min->SetFixedVariable(7,"c",variable[7]);
  min->SetFixedVariable(8,"d",variable[8]);
  min->SetFixedVariable(9,"Mc2",variable[9]);
  // do the minimization
  min->Minimize();
  // get fitting result
  const double * xs = min->X();
  for (int i = 0; i < 10; i++){
    _para[i] = xs[i];
  }
  double para[11];
  GetModelPara(para);
  _para[11] = ChiSquare0(para);
  return 0;
}

int Lfit::A1MinimizerFree(const char * minName = "Minuit", const char * algoName = ""){
  if (_Atype != 1){
    std::cout << "Lfit::A1Minimizer: Wrong asymmetry type!" << std::endl;
    return 1;
  }
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(100000); // for Minuit
  min->SetMaxIterations(10000);  // for GSL 
  min->SetTolerance(0.1);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&Lfit::A1FitFunction, 11);//Fitting function
  double step[11], variable[11];
  InitialA1para();
  for (int i = 0; i < 11; i++){
    step[i] = 0.001;
    variable[i] = _para[i];
  }
  min->SetFunction(f);
  //Set variables to be minimized
  min->SetVariable(0, "Nu", variable[0], step[0]);
  min->SetVariable(1,"Nd",variable[1],step[1]);
  min->SetFixedVariable(2,"Ns",variable[2]);
  min->SetVariable(3,"Nubar",variable[3],step[3]);
  min->SetVariable(4,"Ndbar",variable[4],step[4]);
  min->SetFixedVariable(5,"Nsbar",variable[5]);
  min->SetVariable(6,"au",variable[6],step[6]);
  min->SetVariable(7,"ad",variable[7],step[7]);
  min->SetVariable(8,"asea",variable[8],step[8]);
  min->SetVariable(9,"b",variable[9],step[9]);
  min->SetVariable(10,"Ms2",variable[10],step[10]);
  // do the minimization
  min->Minimize();
  // get fitting result
  const double * xs = min->X();
  for (int i = 0; i < 11; i++){
    _para[i] = xs[i];
  }
  double para[11];
  GetModelPara(para);
  _para[11] = ChiSquare0(para);
  return 0;
}

int Lfit::A2MinimizerFree(const char * minName = "Minuit", const char * algoName = ""){
  if (_Atype != 2){
    std::cout << "Lfit::A2Minimizer: Wrong asymmetry type!" << std::endl;
    return 1;
  }
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(100000); // for Minuit
  min->SetMaxIterations(10000);  // for GSL 
  min->SetTolerance(0.1);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&Lfit::A2FitFunction, 10);//Fitting function
  double step[10], variable[10];
  InitialA2para();
  for (int i = 0; i < 10; i++){
    step[i] = 0.001;
    variable[i] = _para[i];
  }
  min->SetFunction(f);
  //Set variables to be minimized
  min->SetVariable(0, "Nu", variable[0], step[0]);
  min->SetVariable(1,"Nd",variable[1],step[1]);
  min->SetVariable(2,"a",variable[2],step[2]);
  min->SetVariable(3,"b",variable[3],step[3]);
  min->SetFixedVariable(4,"kt2",variable[4]);
  min->SetVariable(5,"Nfav",variable[5], step[5]);
  min->SetVariable(6,"Ndis",variable[6],step[6]);
  min->SetVariable(7,"c",variable[7],step[7]);
  min->SetVariable(8,"d",variable[8],step[8]);
  min->SetVariable(9,"Mc2",variable[9],step[9]);
  // do the minimization
  min->Minimize();
  // get fitting result
  const double * xs = min->X();
  for (int i = 0; i < 10; i++){
    _para[i] = xs[i];
  }
  double para[11];
  GetModelPara(para);
  _para[11] = ChiSquare0(para);
  return 0;
}

int Lfit::A3MinimizerFree(const char * minName = "Minuit", const char * algoName = ""){
  if (_Atype != 3){
    std::cout << "Lfit::A3Minimizer: Wrong asymmetry type!" << std::endl;
    return 1;
  }
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(100000); // for Minuit
  min->SetMaxIterations(10000);  // for GSL 
  min->SetTolerance(0.1);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&Lfit::A3FitFunction, 10);//Fitting function
  double step[10], variable[10];
  InitialA3para();
  for (int i = 0; i < 10; i++){
    step[i] = 0.001;
    variable[i] = _para[i];
  }
  min->SetFunction(f);
  //Set variables to be minimized
  min->SetVariable(0, "Nu", variable[0], step[0]);
  min->SetVariable(1,"Nd",variable[1],step[1]);
  min->SetVariable(2,"a",variable[2],step[2]);
  min->SetFixedVariable(3,"b",variable[3]);
  min->SetVariable(4,"Mt2",variable[4], step[4]);
  min->SetFixedVariable(5,"Nfav",variable[5]);
  min->SetFixedVariable(6,"Ndis",variable[6]);
  min->SetFixedVariable(7,"c",variable[7]);
  min->SetFixedVariable(8,"d",variable[8]);
  min->SetFixedVariable(9,"Mc2",variable[9]);
  // do the minimization
  min->Minimize();
  // get fitting result
  const double * xs = min->X();
  for (int i = 0; i < 10; i++){
    _para[i] = xs[i];
  }
  double para[11];
  GetModelPara(para);
  _para[11] = ChiSquare0(para);
  return 0;
}

int Lfit::FitSivers(const char * filename = "fitsivers.root", const int loops = 1000){
  if (_Atype != 1){
    std::cout << "Lfit::FitSivers: Wrong asymmetry type!" << std::endl;
    return 1;
  }
  TFile * fs = new TFile(filename, "RECREATE");
  TTree * Ts = new TTree("fitpara", "fitpara");
  Ts->SetDirectory(fs);
  double fitpara[11];
  double chi2;
  Ts->Branch("Nu", &fitpara[0], "Nu/D");
  Ts->Branch("Nd", &fitpara[1], "Nd/D");
  Ts->Branch("Ns", &fitpara[2], "Ns/D");
  Ts->Branch("Nubar", &fitpara[3], "Nubar/D");
  Ts->Branch("Ndbar", &fitpara[4], "Ndbar/D");
  Ts->Branch("Nsbar", &fitpara[5], "Nsbar/D");
  Ts->Branch("au", &fitpara[6], "au/D");
  Ts->Branch("ad", &fitpara[7], "ad/D");
  Ts->Branch("asea", &fitpara[8], "asea/D");
  Ts->Branch("b", &fitpara[9], "b/D");
  Ts->Branch("Ms2", &fitpara[10], "Ms2/D");
  Ts->Branch("Chi2", &chi2, "Chi2/D");
  for (int i = 0; i < loops; i++){
    std::cout << "Loop: " << i + 1 << std::endl;
    SmearData();
    A1Minimizer();
    GetModelPara(fitpara);
    chi2 = _para[11];
    Ts->Fill();
  }
  fs->Write();
  return 0;
}

int Lfit::FitCollins(const char * filename = "fitcollins.root", const int loops = 1000){
  if (_Atype != 2){
    std::cout << "Lfit::FitCollins: Wrong asymmetry type!" << std::endl;
    return 1;
  }
  TFile * fs = new TFile(filename, "RECREATE");
  TTree * Ts = new TTree("fitpara", "fitpara");
  Ts->SetDirectory(fs);
  double fitpara[14];
  double chi2;
  Ts->Branch("Nu", &fitpara[0], "Nu/D");
  Ts->Branch("Nd", &fitpara[1], "Nd/D");
  Ts->Branch("au", &fitpara[2], "au/D");
  Ts->Branch("ad", &fitpara[3], "ad/D");
  Ts->Branch("bu", &fitpara[4], "bu/D");
  Ts->Branch("bd", &fitpara[5], "bd/D");
  Ts->Branch("kt2", &fitpara[6], "kt2/D");
  Ts->Branch("Nfav", &fitpara[7], "Nfav/D");
  Ts->Branch("Ndis", &fitpara[8], "Ndis/D");
  Ts->Branch("cu", &fitpara[9], "cu/D");
  Ts->Branch("cd", &fitpara[10], "cd/D");
  Ts->Branch("du", &fitpara[11], "du/D");
  Ts->Branch("dd", &fitpara[12], "dd/D");
  Ts->Branch("Mc2", &fitpara[13], "Mc2/D");
  Ts->Branch("Chi2", &chi2, "Chi2/D"); 
  for (int i = 0; i < loops; i++){
    std::cout << "Loop: " << i + 1 << std::endl;
    SmearData();
    A2Minimizer();
    GetModelPara(fitpara);
    chi2 = _para[11];
    Ts->Fill();
  }
  fs->Write();
  return 0;
}

int Lfit::FitPretzelosity(const char * filename = "fitpretzelosity.root", const int loops = 1000){
  if (_Atype != 3){
    std::cout << "Lfit::FitPretzelosity: Wrong asymmetry type!" << std::endl;
    return 1;
  }
  TFile * fs = new TFile(filename, "RECREATE");
  TTree * Ts = new TTree("fitpara", "fitpara");
  Ts->SetDirectory(fs);
  double fitpara[14];
  double chi2;
  Ts->Branch("Nu", &fitpara[0], "Nu/D");
  Ts->Branch("Nd", &fitpara[1], "Nd/D");
  Ts->Branch("au", &fitpara[2], "au/D");
  Ts->Branch("ad", &fitpara[3], "ad/D");
  Ts->Branch("bu", &fitpara[4], "bu/D");
  Ts->Branch("bd", &fitpara[5], "bd/D");
  Ts->Branch("Mt2", &fitpara[6], "Mt2/D");
  Ts->Branch("Nfav", &fitpara[7], "Nfav/D");
  Ts->Branch("Ndis", &fitpara[8], "Ndis/D");
  Ts->Branch("cu", &fitpara[9], "cu/D");
  Ts->Branch("cd", &fitpara[10], "cd/D");
  Ts->Branch("du", &fitpara[11], "du/D");
  Ts->Branch("dd", &fitpara[12], "dd/D");
  Ts->Branch("Mc2", &fitpara[13], "Mc2/D");
  Ts->Branch("Chi2", &chi2, "Chi2/D"); 
  for (int i = 0; i < loops; i++){
    std::cout << "Loop: " << i + 1 << std::endl;
    SmearData();
    A3Minimizer();
    GetModelPara(fitpara);
    chi2 = _para[11];
    Ts->Fill();
  }
  fs->Write();
  return 0;
}

int Lfit::FitSiversFree(const char * filename = "fitsiversfree.root", const int loops = 1000){
  if (_Atype != 1){
    std::cout << "Lfit::FitSivers: Wrong asymmetry type!" << std::endl;
    return 1;
  }
  TFile * fs = new TFile(filename, "RECREATE");
  TTree * Ts = new TTree("fitpara", "fitpara");
  Ts->SetDirectory(fs);
  double fitpara[11];
  double chi2;
  Ts->Branch("Nu", &fitpara[0], "Nu/D");
  Ts->Branch("Nd", &fitpara[1], "Nd/D");
  Ts->Branch("Ns", &fitpara[2], "Ns/D");
  Ts->Branch("Nubar", &fitpara[3], "Nubar/D");
  Ts->Branch("Ndbar", &fitpara[4], "Ndbar/D");
  Ts->Branch("Nsbar", &fitpara[5], "Nsbar/D");
  Ts->Branch("au", &fitpara[6], "au/D");
  Ts->Branch("ad", &fitpara[7], "ad/D");
  Ts->Branch("asea", &fitpara[8], "asea/D");
  Ts->Branch("b", &fitpara[9], "b/D");
  Ts->Branch("Ms2", &fitpara[10], "Ms2/D");
  Ts->Branch("Chi2", &chi2, "Chi2/D");
  for (int i = 0; i < loops; i++){
    std::cout << "Loop: " << i + 1 << std::endl;
    SmearData();
    A1MinimizerFree();
    GetModelPara(fitpara);
    //if (fitpara[0]*fitpara[0] > 1.001 || fitpara[1]*fitpara[1] > 1.001) continue;
    //if (fitpara[3]*fitpara[3] > 1.001 || fitpara[4]*fitpara[4] > 1.001) continue;
    chi2 = _para[11];
    Ts->Fill();
  }
  fs->Write();
  return 0;
}

int Lfit::FitCollinsFree(const char * filename = "fitcollinsfree.root", const int loops = 1000){
  if (_Atype != 2){
    std::cout << "Lfit::FitCollins: Wrong asymmetry type!" << std::endl;
    return 1;
  }
  TFile * fs = new TFile(filename, "RECREATE");
  TTree * Ts = new TTree("fitpara", "fitpara");
  Ts->SetDirectory(fs);
  double fitpara[14];
  double chi2;
  Ts->Branch("Nu", &fitpara[0], "Nu/D");
  Ts->Branch("Nd", &fitpara[1], "Nd/D");
  Ts->Branch("au", &fitpara[2], "au/D");
  Ts->Branch("ad", &fitpara[3], "ad/D");
  Ts->Branch("bu", &fitpara[4], "bu/D");
  Ts->Branch("bd", &fitpara[5], "bd/D");
  Ts->Branch("kt2", &fitpara[6], "kt2/D");
  Ts->Branch("Nfav", &fitpara[7], "Nfav/D");
  Ts->Branch("Ndis", &fitpara[8], "Ndis/D");
  Ts->Branch("cu", &fitpara[9], "cu/D");
  Ts->Branch("cd", &fitpara[10], "cd/D");
  Ts->Branch("du", &fitpara[11], "du/D");
  Ts->Branch("dd", &fitpara[12], "dd/D");
  Ts->Branch("Mc2", &fitpara[13], "Mc2/D");
  Ts->Branch("Chi2", &chi2, "Chi2/D"); 
  for (int i = 0; i < loops; i++){
    std::cout << "Loop: " << i + 1 << std::endl;
    SmearData();
    A2MinimizerFree();
    GetModelPara(fitpara);
    //if (fitpara[0]*fitpara[0] > 1.001 || fitpara[1]*fitpara[1] > 1.001) continue;
    //if (fitpara[7]*fitpara[7] > 1.001 || fitpara[8]*fitpara[8] > 1.001) continue;
    chi2 = _para[11];
    Ts->Fill();
  }
  fs->Write();
  return 0;
}

int Lfit::FitPretzelosityFree(const char * filename = "fitpretzelosityfree.root", const int loops = 1000){
  if (_Atype != 3){
    std::cout << "Lfit::FitPretzelosity: Wrong asymmetry type!" << std::endl;
    return 1;
  }
  TFile * fs = new TFile(filename, "RECREATE");
  TTree * Ts = new TTree("fitpara", "fitpara");
  Ts->SetDirectory(fs);
  double fitpara[14];
  double chi2;
  Ts->Branch("Nu", &fitpara[0], "Nu/D");
  Ts->Branch("Nd", &fitpara[1], "Nd/D");
  Ts->Branch("au", &fitpara[2], "au/D");
  Ts->Branch("ad", &fitpara[3], "ad/D");
  Ts->Branch("bu", &fitpara[4], "bu/D");
  Ts->Branch("bd", &fitpara[5], "bd/D");
  Ts->Branch("Mt2", &fitpara[6], "Mt2/D");
  Ts->Branch("Nfav", &fitpara[7], "Nfav/D");
  Ts->Branch("Ndis", &fitpara[8], "Ndis/D");
  Ts->Branch("cu", &fitpara[9], "cu/D");
  Ts->Branch("cd", &fitpara[10], "cd/D");
  Ts->Branch("du", &fitpara[11], "du/D");
  Ts->Branch("dd", &fitpara[12], "dd/D");
  Ts->Branch("Mc2", &fitpara[13], "Mc2/D");
  Ts->Branch("Chi2", &chi2, "Chi2/D"); 
  for (int i = 0; i < loops; i++){
    std::cout << "Loop: " << i + 1 << std::endl;
    SmearData();
    A3MinimizerFree();
    GetModelPara(fitpara);
    //if (fitpara[0]*fitpara[0] > 1.001 || fitpara[1]*fitpara[1] > 1.001) continue;
    chi2 = _para[11];
    Ts->Fill();
  }
  fs->Write();
  return 0;
}

/* dchi2: 95.45%
   1:   4.0000
   2:   6.1801
   3:   8.0249
   4:   9.7156
   5:  11.3139
   6:  12.8489
   7:  14.3371
   8:  15.7891
   9:  17.2118
   10: 18.6104
   11: 19.9884
 */


int Lfit::GetFitf1tM(const double * var, double f1tM[3][6], const double dchi2 = 17.2118, const char * file = "fitsivers.root"){
  TFile * ff = new TFile(file,"r");
  TTree * Tf = (TTree *) ff->GetObjectChecked("fitpara","TTree");
  int Nentry = Tf->GetEntries();
  double fpara[11];
  double chi2;
  Tf->SetBranchAddress("Nu", &fpara[0]);
  Tf->SetBranchAddress("Nd", &fpara[1]);
  Tf->SetBranchAddress("Ns", &fpara[2]);
  Tf->SetBranchAddress("Nubar", &fpara[3]);
  Tf->SetBranchAddress("Ndbar", &fpara[4]);
  Tf->SetBranchAddress("Nsbar", &fpara[5]);
  Tf->SetBranchAddress("au", &fpara[6]);
  Tf->SetBranchAddress("ad", &fpara[7]);
  Tf->SetBranchAddress("asea", &fpara[8]);
  Tf->SetBranchAddress("b", &fpara[9]);
  Tf->SetBranchAddress("Ms2", &fpara[10]);
  Tf->SetBranchAddress("Chi2", &chi2);
  lst->f1tM_std(var, f1tM[0]);//center values
  for(int j = 0; j < 6; j++){
    f1tM[1][j] = f1tM[0][j];//upper values
    f1tM[2][j] = f1tM[0][j];//lower values
  }
  double temp[6];
  for (int i = 0; i < Nentry; i++){
    Tf->GetEntry(i);
    if (chi2 > dchi2) continue;
    lst->f1tM_std(var, temp, fpara);
    for (int j = 0; j < 6; j++){
      if (temp[j] > f1tM[1][j]) f1tM[1][j] = temp[j];
      if (temp[j] < f1tM[2][j]) f1tM[2][j] = temp[j];
    }
  }
  return 0;
}

int Lfit::GetFittmdf1t(const double * var, double tmdf1t[3][6], const double dchi2 = 17.2118, const char * file = "fitsivers.root"){
  TFile * ff = new TFile(file,"r");
  TTree * Tf = (TTree *) ff->GetObjectChecked("fitpara","TTree");
  int Nentry = Tf->GetEntries();
  double fpara[11];
  double chi2;
  Tf->SetBranchAddress("Nu", &fpara[0]);
  Tf->SetBranchAddress("Nd", &fpara[1]);
  Tf->SetBranchAddress("Ns", &fpara[2]);
  Tf->SetBranchAddress("Nubar", &fpara[3]);
  Tf->SetBranchAddress("Ndbar", &fpara[4]);
  Tf->SetBranchAddress("Nsbar", &fpara[5]);
  Tf->SetBranchAddress("au", &fpara[6]);
  Tf->SetBranchAddress("ad", &fpara[7]);
  Tf->SetBranchAddress("asea", &fpara[8]);
  Tf->SetBranchAddress("b", &fpara[9]);
  Tf->SetBranchAddress("Ms2", &fpara[10]);
  Tf->SetBranchAddress("Chi2", &chi2);
  lst->tmd_f1t_std(var, tmdf1t[0]);//center values
  for(int j = 0; j < 6; j++){
    tmdf1t[1][j] = tmdf1t[0][j];//upper values
    tmdf1t[2][j] = tmdf1t[0][j];//lower values
  }
  double temp[6];
  for (int i = 0; i < Nentry; i++){
    Tf->GetEntry(i);
    if (chi2 > dchi2) continue;
    lst->tmd_f1t_std(var, temp, fpara);
    for (int j = 0; j < 6; j++){
      if (temp[j] > tmdf1t[1][j]) tmdf1t[1][j] = temp[j];
      if (temp[j] < tmdf1t[2][j]) tmdf1t[2][j] = temp[j];
    }
  }
  return 0;
}

int Lfit::GetFitf1tMcharge(const double Q2, double tc[3][8], const double dchi2 = 17.2118, const char * file = "fitsivers.root"){
  TFile * ff = new TFile(file,"r");
  TTree * Tf = (TTree *) ff->GetObjectChecked("fitpara","TTree");
  int Nentry = Tf->GetEntries();
  double fpara[11];
  double chi2;
  Tf->SetBranchAddress("Nu", &fpara[0]);
  Tf->SetBranchAddress("Nd", &fpara[1]);
  Tf->SetBranchAddress("Ns", &fpara[2]);
  Tf->SetBranchAddress("Nubar", &fpara[3]);
  Tf->SetBranchAddress("Ndbar", &fpara[4]);
  Tf->SetBranchAddress("Nsbar", &fpara[5]);
  Tf->SetBranchAddress("au", &fpara[6]);
  Tf->SetBranchAddress("ad", &fpara[7]);
  Tf->SetBranchAddress("asea", &fpara[8]);
  Tf->SetBranchAddress("b", &fpara[9]);
  Tf->SetBranchAddress("Ms2", &fpara[10]);
  Tf->SetBranchAddress("Chi2", &chi2);
  lst->f1tMcharge_std(Q2, tc[0]);//center values
  for(int j = 0; j < 8; j++){
    tc[1][j] = tc[0][j];//upper values
    tc[2][j] = tc[0][j];//lower values
  }
  double temp[8];
  for (int i = 0; i < Nentry; i++){
    Tf->GetEntry(i);
    if (chi2 > dchi2) continue;
    lst->f1tMcharge_std(Q2, temp, fpara);
    for (int j = 0; j < 8; j++){
      if (temp[j] > tc[1][j]) tc[1][j] = temp[j];
      if (temp[j] < tc[2][j]) tc[2][j] = temp[j];
    }
  }
  return 0;
}
  
int Lfit::GetFith1(const double * var, double h1[3][6], const double dchi2 = 17.2118, const char * file = "fitcollins.root"){
  TFile * ff = new TFile(file,"r");
  TTree * Tf = (TTree *) ff->GetObjectChecked("fitpara","TTree");
  int Nentry = Tf->GetEntries();
  double hpara[7];
  double chi2;
  Tf->SetBranchAddress("Nu", &hpara[0]);
  Tf->SetBranchAddress("Nd", &hpara[1]);
  Tf->SetBranchAddress("au", &hpara[2]);
  Tf->SetBranchAddress("ad", &hpara[3]);
  Tf->SetBranchAddress("bu", &hpara[4]);
  Tf->SetBranchAddress("bd", &hpara[5]);
  Tf->SetBranchAddress("kt2", &hpara[6]);
  Tf->SetBranchAddress("Chi2", &chi2);
  lst->h1_std(var, h1[0]);//center values
  for(int j = 0; j < 6; j++){
    h1[1][j] = h1[0][j];//upper values
    h1[2][j] = h1[0][j];//lower values
  }
  double temp[6];
  for (int i = 0; i < Nentry; i++){
    Tf->GetEntry(i);
    if (chi2 > dchi2) continue;
    lst->h1_std(var, temp, hpara);
    if (hpara[0]*hpara[0] < 1.001){
      if (temp[0] > h1[1][0]) h1[1][0] = temp[0];
      if (temp[0] < h1[2][0]) h1[2][0] = temp[0];
    }
    if (hpara[1]*hpara[1] < 1.001){
      if (temp[1] > h1[1][1]) h1[1][1] = temp[1];
      if (temp[1] < h1[2][1]) h1[2][1] = temp[1];
    }
  }
  return 0;
}  

int Lfit::GetFith1charge(const double Q2, double tc[3][6], const double dchi2 = 17.2118, const char * file = "fitcollins.root"){
  TFile * ff = new TFile(file,"r");
  TTree * Tf = (TTree *) ff->GetObjectChecked("fitpara","TTree");
  int Nentry = Tf->GetEntries();
  double hpara[7];
  double chi2;
  Tf->SetBranchAddress("Nu", &hpara[0]);
  Tf->SetBranchAddress("Nd", &hpara[1]);
  Tf->SetBranchAddress("au", &hpara[2]);
  Tf->SetBranchAddress("ad", &hpara[3]);
  Tf->SetBranchAddress("bu", &hpara[4]);
  Tf->SetBranchAddress("bd", &hpara[5]);
  Tf->SetBranchAddress("kt2", &hpara[6]);
  Tf->SetBranchAddress("Chi2", &chi2);
  lst->h1charge_std(Q2, tc[0]);
  for(int j = 0; j < 6; j++){
    tc[1][j] = tc[0][j];//upper values
    tc[2][j] = tc[0][j];//lower values
  }
  double temp[6];
  for (int i = 0; i < Nentry; i++){
    Tf->GetEntry(i);
    if (chi2 > dchi2) continue;
    lst->h1charge_std(Q2, temp, hpara);
    if (hpara[0]*hpara[0] < 1.001){
      if (temp[0] > tc[1][0]) tc[1][0] = temp[0];
      if (temp[0] < tc[2][0]) tc[2][0] = temp[0];
    }
    if (hpara[1]*hpara[1] < 1.001){
      if (temp[1] > tc[1][1]) tc[1][1] = temp[1];
      if (temp[1] < tc[2][1]) tc[2][1] = temp[1];
    }
  }
  return 0;
}
  
  
int Lfit::GetFith1tM(const double * var, double h1tM[3][6], const double dchi2 = 9.7156, const char * file = "fitpretzelosity.root"){
  TFile * ff = new TFile(file,"r");
  TTree * Tf = (TTree *) ff->GetObjectChecked("fitpara","TTree");
  int Nentry = Tf->GetEntries();
  double hpara[7];
  double chi2;
  Tf->SetBranchAddress("Nu", &hpara[0]);
  Tf->SetBranchAddress("Nd", &hpara[1]);
  Tf->SetBranchAddress("au", &hpara[2]);
  Tf->SetBranchAddress("ad", &hpara[3]);
  Tf->SetBranchAddress("bu", &hpara[4]);
  Tf->SetBranchAddress("bd", &hpara[5]);
  Tf->SetBranchAddress("Mt2", &hpara[6]);
  Tf->SetBranchAddress("Chi2", &chi2);
  lst->h1tM_std(var, h1tM[0]);//center values
  for(int j = 0; j < 6; j++){
    h1tM[1][j] = h1tM[0][j];//upper values
    h1tM[2][j] = h1tM[0][j];//lower values
  }
  double temp[6];
  for (int i = 0; i < Nentry; i++){
    Tf->GetEntry(i);
    if (chi2 > dchi2) continue;
    lst->h1tM_std(var, temp, hpara);
    if (hpara[0]*hpara[0] < 1.001){
      if (temp[0] > h1tM[1][0]) h1tM[1][0] = temp[0];
      if (temp[0] < h1tM[2][0]) h1tM[2][0] = temp[0];
    }
    if (hpara[1]*hpara[1] < 1.001){
      if (temp[1] > h1tM[1][1]) h1tM[1][1] = temp[1];
      if (temp[1] < h1tM[2][1]) h1tM[2][1] = temp[1];
    }
  }
  return 0;
}

int Lfit::GetFith1tMcharge(const double Q2, double tc[3][6], const double dchi2 = 9.7156, const char * file = "fitpretzelosity.root"){
  TFile * ff = new TFile(file,"r");
  TTree * Tf = (TTree *) ff->GetObjectChecked("fitpara","TTree");
  int Nentry = Tf->GetEntries();
  double hpara[7];
  double chi2;
  Tf->SetBranchAddress("Nu", &hpara[0]);
  Tf->SetBranchAddress("Nd", &hpara[1]);
  Tf->SetBranchAddress("au", &hpara[2]);
  Tf->SetBranchAddress("ad", &hpara[3]);
  Tf->SetBranchAddress("bu", &hpara[4]);
  Tf->SetBranchAddress("bd", &hpara[5]);
  Tf->SetBranchAddress("Mt2", &hpara[6]);
  Tf->SetBranchAddress("Chi2", &chi2);
  lst->h1tMcharge_std(Q2, tc[0]);//center values
  for(int j = 0; j < 6; j++){
    tc[1][j] = tc[0][j];//upper values
    tc[2][j] = tc[0][j];//lower values
  }
  double temp[6];
  for (int i = 0; i < Nentry; i++){
    Tf->GetEntry(i);
    if (chi2 > dchi2) continue;
    lst->h1tMcharge_std(Q2, temp, hpara);
    if (hpara[0]*hpara[0] < 1.001){
      if (temp[0] > tc[1][0]) tc[1][0] = temp[0];
      if (temp[0] < tc[2][0]) tc[2][0] = temp[0];
    }
    if (hpara[1]*hpara[1] < 1.001){
      if (temp[1] > tc[1][1]) tc[1][1] = temp[1];
      if (temp[1] < tc[2][1]) tc[2][1] = temp[1];
    }
  }
  return 0;
}

#endif
