#ifndef _LSAMPLE_H_
#define _LSAMPLE_H_

#include <iostream>
#include <fstream>
#include <cmath>
#include "TRandom.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TH1D.h"
#include "TString.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TTree.h"
#include "TFile.h"
#include "Lstructure.h"

#define _NMAX_ 4000

class Lsample{
 protected:
  int _Atype;//1: Siver, 2: Collins, 3: pretzelosity
  int _Ndata;
  double _Nucleon[_NMAX_];
  double _Hadron[_NMAX_];
  double _kin[_NMAX_][5];
  double _A0[_NMAX_];
  double _E0[_NMAX_];
  double _ET[_NMAX_];
  double _A[_NMAX_];
  int _PSD;//parameter space dimension
  double _hessian[11][11];
  double _smatrix[11][11];
  double _eigenvalues[11];
  double _para[11];
  double _central[11];
  int _err;
 public:
  Lsample(int Atype);
  int InitialA1para();
  int InitialA2para();
  int InitialA3para();
  int GetData(const char * datafile, const int Atype);
  int PrintData(const int n);
  int PrintStatus();
  double Chi2A1(const int err, const double * fitpara);
  double Chi2A2(const int err, const double * fitpara);
  double Chi2A3(const int err, const double * fitpara);
  double SetHessianA1(const int err);
  double SetHessianA2(const int err);
  double SetHessianA3(const int err);
  int LikeliSampleA1(const int _err, const int calls, const char * savefile);
  int LikeliSampleA2(const int _err, const int calls, const char * savefile);
  int LikeliSampleA3(const int _err, const int calls, const char * savefile);
  int UniformSampleA1(const int _err, const int calls, const char * savefile);
  int UniformSampleA2(const int _err, const int calls, const char * savefile);
  int UniformSampleA3(const int _err, const int calls, const char * savefile);
  double FindChi2Limit(const char * paraset, const double CL, const int cut);
};

Lsample::Lsample(int Atype){
  _Atype = Atype;
  _Ndata = 0;
  if (Atype == 1) _PSD = 9;
  else if (Atype == 2) _PSD = 9;
  else if (Atype == 3) _PSD = 4;
  else {
    _PSD = 0;
    std::cout << "Lsample initialize error: invalid asymmetry type!" << std::endl;
  }
}

int Lsample::InitialA1para(){
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
  //central
  _central[0] = _para[0];//Nu
  _central[1] = _para[1];//Nd
  _central[2] = _para[3];//Nubar
  _central[3] = _para[4];//Ndbar
  _central[4] = _para[6];//au
  _central[5] = _para[7];//ad
  _central[6] = _para[8];//asea
  _central[7] = _para[9];//b
  _central[8] = _para[10];//MS2
  return 0;
}

int Lsample::InitialA2para(){
  //transversity
  _para[0] = 0.46;//Nu
  _para[1] = -1.0;//Nd
  _para[2] = 1.11;//a
  _para[3] = 3.64;//b
  //Collins
  _para[4] = 0.49;//Nfav
  _para[5] = -1.0;//Ndis
  _para[6] = 1.06;//c
  _para[7] = 0.07;//d
  _para[8] = 1.5;//MC2
  //central
  _central[0] = _para[0];//Nu
  _central[1] = _para[1];//Nd
  _central[2] = _para[2];//a
  _central[3] = _para[3];//b
  _central[4] = _para[4];//Nfav
  _central[5] = _para[5];//Ndis
  _central[6] = _para[6];//c
  _central[7] = _para[7];//d
  _central[8] = _para[8];//MC2
  return 0;
}

int Lsample::InitialA3para(){
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
  //central
  _central[0] = _para[0];//Nu
  _central[1] = _para[1];//Nd
  _central[2] = _para[2];//a
  _central[3] = _para[4];//Mh2
  return 0;
}
 
int Lsample::GetData(const char * datafile, const int Atype){
  if (Atype != _Atype && _Ndata != 0){
    std::cout << "Lsample::GetData: dismatch Atype!" << std::endl;
    return 1;
  }
  _Atype = Atype;
  TFile * fd = new TFile(datafile, "r");
  TTree * Td = (TTree *) fd->Get("data");
  double ndata = Td->GetEntries();
  double Nucleon, Hadron, Ebeam;
  double kin[5];
  double Asym[3], Estat[3], Eabs[3], Erel[3];
  Td->SetBranchAddress("Nucleon", &Nucleon);
  Td->SetBranchAddress("Hadron", &Hadron);
  Td->SetBranchAddress("Ebeam", &Ebeam);
  Td->SetBranchAddress("x", &kin[0]);
  Td->SetBranchAddress("y", &kin[1]);
  Td->SetBranchAddress("z", &kin[2]);
  Td->SetBranchAddress("Q2", &kin[3]);
  Td->SetBranchAddress("Pt", &kin[4]);
  Td->SetBranchAddress("A1", &Asym[0]);
  Td->SetBranchAddress("A2", &Asym[1]);
  Td->SetBranchAddress("A3", &Asym[2]);
  Td->SetBranchAddress("E1stat", &Estat[0]);
  Td->SetBranchAddress("E2stat", &Estat[1]);
  Td->SetBranchAddress("E3stat", &Estat[2]);
  Td->SetBranchAddress("E1abs", &Eabs[0]);
  Td->SetBranchAddress("E2abs", &Eabs[1]);
  Td->SetBranchAddress("E3abs", &Eabs[2]);
  Td->SetBranchAddress("E1rel", &Erel[0]);
  Td->SetBranchAddress("E2rel", &Erel[1]);
  Td->SetBranchAddress("E3rel", &Erel[2]);
  int skip = 0;
  for (int i = 0; i < ndata; i++){
    Td->GetEntry(i);
    _Nucleon[_Ndata + i - skip] = Nucleon;
    _Hadron[_Ndata + i - skip] = Hadron;
    for (int j = 0; j < 5; j++){
      _kin[_Ndata + i - skip][j] = kin[j];
    }
    _A0[_Ndata + i - skip] = Asym[Atype - 1];
    _E0[_Ndata + i - skip] = Estat[Atype - 1];
    _ET[_Ndata + i - skip] = sqrt(Estat[Atype - 1] * Estat[Atype - 1] + Eabs[Atype - 1] * Eabs[Atype - 1] + Erel[Atype - 1] * Erel[Atype - 1] * Asym[Atype - 1] * Asym[Atype - 1]);
    if (_ET[_Ndata + i - skip] < 1.0e5) continue;
    else skip++;
  }
  _Ndata = _Ndata + ndata - skip;
  std::cout << "skip: " << skip << std::endl;
  return 0;
}

int Lsample::PrintData(const int n){
  if (n >= _Ndata){
    std::cout << "Lsample::PrintData: overflow!" << std::endl;
    return 1;
  }
  printf("   x       y       z       Q2      Pt        A      Estat     Etotal\n");
  printf("%.4f  %.4f  %.4f  %.4f  %.4f  %.6f  %.6f  %.6f\n", _kin[n][0], _kin[n][1], _kin[n][2], _kin[n][3], _kin[n][4], _A0[n], _E0[n], _ET[n]);
  return 0;
}

int Lsample::PrintStatus(){
  std::cout << "Asymmetry type: " << _Atype << std::endl;
  std::cout << "Data number   : " << _Ndata << std::endl;
  return 0;
}

double Lsample::Chi2A2(const int err = 1, const double * fitpara = 0){
  double Tf[7], TD[7];
  if (fitpara == 0){
    Tf[0] = _central[0];//Nu
    Tf[1] = _central[1];//Nd
    Tf[2] = _central[2];//au
    Tf[3] = _central[2];//ad
    Tf[4] = _central[3];//bu
    Tf[5] = _central[3];//bd
    Tf[6] = 0.25;
    TD[0] = _central[4];//Nfav
    TD[1] = _central[5];//Ndis
    TD[2] = _central[6];//c
    TD[3] = _central[6];//c
    TD[4] = _central[7];//d
    TD[5] = _central[7];//d
    TD[6] = _central[8];//MC2
  }
  else {
    Tf[0] = fitpara[0];//Nu
    Tf[1] = fitpara[1];//Nd
    Tf[2] = fitpara[2];//au
    Tf[3] = fitpara[2];//ad
    Tf[4] = fitpara[3];//bu
    Tf[5] = fitpara[3];//bd
    Tf[6] = 0.25;
    TD[0] = fitpara[4];//Nfav
    TD[1] = fitpara[5];//Ndis
    TD[2] = fitpara[6];//c
    TD[3] = fitpara[6];//c
    TD[4] = fitpara[7];//d
    TD[5] = fitpara[7];//d
    TD[6] = fitpara[8];//MC2
  }
  double Asym[2];
  double sum = 0.0;
  for (int i = 0; i < _Ndata; i++){
    if (_Nucleon[i] == 0) Lstructure::AsinHpSn(_kin[i], Asym, Tf, TD);
    else if (_Nucleon[i] == 1) Lstructure::AsinHpSp(_kin[i], Asym, Tf, TD);
    _A[i] = Asym[(int)_Hadron[i]];
    if (err == 0){
      sum = sum + pow( (_A[i] - _A0[i]) / _E0[i], 2);
    }
    else {
      sum = sum + pow( (_A[i] - _A0[i]) / _ET[i], 2);
    }
  }
  return sum;
}

double Lsample::SetHessianA2(const int err){
  _err = err;
  double da = 1.0e-6;
  double db = 1.0e-6;
  double chi2min = Chi2A2(err);
  std::cout << "chi2min: " << chi2min << std::endl;
  double newpara[9];
  double chi2[4];
  TMatrixD hessian(9,9);
  double scale[9] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  //double scale[9] = {5.0e-3, 1.0e-2, 1.0e-2, 5.0e-2, 5.0e-3, 1.0e-2, 1.0e-2, 2.0e-3, 1.0e-2};
  for (int i = 0; i < 9; i++){
    for (int j = 0; j < 9; j++){
      for (int s = 0; s < 9; s++) newpara[s] = _central[s];
      newpara[i] = newpara[i] + da * scale[i];
      newpara[j] = newpara[j] + db * scale[j];
      chi2[0] = Chi2A2(err, newpara);//a+da, b+db
      newpara[j] = newpara[j] - 2.0 * db * scale[j];
      chi2[1] = Chi2A2(err, newpara);//a+da, b-db
      newpara[i] = newpara[i] - 2.0 * da * scale[i];
      chi2[3] = Chi2A2(err, newpara);//a-da, b-db
      newpara[j] = newpara[j] + 2.0 * db * scale[j];
      chi2[2] = Chi2A2(err, newpara);//a-da, b+db
      if (i == j) hessian(i, i) = (chi2[0] + chi2[3] - 2.0 * chi2min) / (2.0 * pow((da+db), 2));
      else hessian(i, j) = (chi2[0] + chi2[3] - chi2[1] - chi2[2]) / (8.0 * da * db);
      _hessian[i][j] = hessian(i, j);
    }
  }
  TMatrixDEigen eigen(hessian);
  TMatrixD ev = eigen.GetEigenValues();
  TMatrixD smatrix = eigen.GetEigenVectors();
  for (int i = 0; i < 9; i++){
    _eigenvalues[i] = ev(i, i);
    for (int j = 0; j < 9; j++){
      _smatrix[i][j] = smatrix(i, j);
    }
  }
  /* hessian.Print(); */
  /* ev.Print(); */
  /* smatrix.Print(); */
  /* std::cout << hessian.Determinant() << std::endl; */
  return 0;
}

int Lsample::LikeliSampleA2(const int err, const int calls, const char * savefile){
  _err = err;
  TFile * fs = new TFile(savefile, "RECREATE");
  TTree * Ts = new TTree("fitpara", "fitpara");
  Ts->SetDirectory(fs);
  TMatrixD as(9, 1);
  TMatrixD bs(9, 1);
  TMatrixD smatrix(9, 9);
  TMatrixD ev(9, 9);
  double para[9];
  double weight, bound, chi2;
  double kt2 = 0.25;
  Ts->Branch("Nu", &para[0], "Nu/D");
  Ts->Branch("Nd", &para[1], "Nd/D");
  Ts->Branch("au", &para[2], "au/D");
  Ts->Branch("ad", &para[2], "ad/D");
  Ts->Branch("bu", &para[3], "bu/D");
  Ts->Branch("bd", &para[3], "bd/D");
  Ts->Branch("kt2", &kt2, "kt2/D");
  Ts->Branch("Nfav", &para[4], "Nfav/D");
  Ts->Branch("Ndis", &para[5], "Ndis/D");
  Ts->Branch("cu", &para[6], "cu/D");
  Ts->Branch("cd", &para[6], "cd/D");
  Ts->Branch("du", &para[7], "du/D");
  Ts->Branch("dd", &para[7], "dd/D");
  Ts->Branch("MC2", &para[8], "MC2/D");
  Ts->Branch("Chi2", &chi2, "Chi2/D");
  Ts->Branch("weight", &weight, "weight/D");
  Ts->Branch("bound", &bound, "bound/D");
  double cup;
  //gRandom->SetSeed(0);
  for (int i = 0; i < 9; i++){
    ev(i, i) = _eigenvalues[i];
    for (int j = 0; j < 9; j++){
      smatrix(i, j) = _smatrix[i][j];
    }
  }
  for (Long64_t n = 0; n < calls; n++){
    std::cout << "#" << n << std::endl;
    cup = 0;
    for (int j = 0; j < 9; j++){
      bs(j,0) = gRandom->Gaus(0.0, 1.0/sqrt(ev(j,j)) );
      cup = cup + bs(j,0)*bs(j,0)*ev(j,j)/2.0;
    }
    bs.Print();
    as = smatrix * bs;
    as.Print();
    as = smatrix * bs;
    as.Print();
    for (int i = 0; i < 9; i++){
      para[i] = _central[i] + as(i, 0);
    }
    if (std::abs(para[0]) > 1.0 || std::abs(para[1]) > 1.0) bound = 0;
    else if (std::abs(para[4]) > 1.0 || std::abs(para[5]) > 1.0) bound = 0;
    else bound = 1;
    chi2 = Chi2A2(_err, para);
    weight = exp(cup-chi2/2.0);
    Ts->Fill();
  }
  fs->Write();
  return 0;
}
    
int Lsample::UniformSampleA2(const int err, const int calls, const char * savefile){
  _err = err;
  TFile * fs = new TFile(savefile, "RECREATE");
  TTree * Ts = new TTree("fitpara", "fitpara");
  Ts->SetDirectory(fs);
  double para[9];
  double shift[9];
  double width[9] = {0.025, 0.05, 0.025, 0.11, 0.03, 0.048, 0.055, 0.07, 0.06};
  double weight, bound, chi2;
  double kt2 = 0.25;
  Ts->Branch("Nu", &para[0], "Nu/D");
  Ts->Branch("Nd", &para[1], "Nd/D");
  Ts->Branch("au", &para[2], "au/D");
  Ts->Branch("ad", &para[2], "ad/D");
  Ts->Branch("bu", &para[3], "bu/D");
  Ts->Branch("bd", &para[3], "bd/D");
  Ts->Branch("kt2", &kt2, "kt2/D");
  Ts->Branch("Nfav", &para[4], "Nfav/D");
  Ts->Branch("Ndis", &para[5], "Ndis/D");
  Ts->Branch("cu", &para[6], "cu/D");
  Ts->Branch("cd", &para[6], "cd/D");
  Ts->Branch("du", &para[7], "du/D");
  Ts->Branch("dd", &para[7], "dd/D");
  Ts->Branch("MC2", &para[8], "MC2/D");
  Ts->Branch("Chi2", &chi2, "Chi2/D");
  Ts->Branch("weight", &weight, "weight/D");
  Ts->Branch("bound", &bound, "bound/D");
  double cup;
  gRandom->SetSeed(0);  
  for (Long64_t n = 0; n < calls; n++){
    if (n%100 == 0) std::cout << "#" << n << std::endl;
    shift[0] = gRandom->Gaus(0.0, width[0]);
    shift[1] = gRandom->Gaus(0.0, width[1]);
    shift[2] = gRandom->Gaus(0.0, width[2]);
    shift[3] = gRandom->Gaus(0.0, width[3]);
    shift[4] = gRandom->Gaus(0.0, width[4]);
    shift[5] = gRandom->Gaus(0.0, width[5]);
    shift[6] = gRandom->Gaus(0.0, width[6]);
    shift[7] = gRandom->Gaus(0.0, width[7]);
    shift[8] = gRandom->Gaus(0.0, width[8]);
    cup = 
      (shift[0]*shift[0]/width[0]/width[0] 
      + shift[1]*shift[1]/width[1]/width[1]
      + shift[2]*shift[2]/width[2]/width[2]
      + shift[3]*shift[3]/width[3]/width[3]
      + shift[4]*shift[4]/width[4]/width[4]
      + shift[5]*shift[5]/width[5]/width[5]
      + shift[6]*shift[6]/width[6]/width[6]
      + shift[7]*shift[7]/width[7]/width[7]
       + shift[8]*shift[8]/width[8]/width[8]) / 2.0;
    for (int i = 0; i < 9; i++){
      para[i] = _central[i] + shift[i];
    }
    if (para[7] < 0.0) continue;
    if (std::abs(para[0]) > 1.0 || std::abs(para[1]) > 1.0) bound = 0;
    else if (std::abs(para[4]) > 1.0 || std::abs(para[5]) > 1.0) bound = 0;
    else bound = 1;
    chi2 = Chi2A2(_err, para);
    weight = exp(cup - chi2/2.0);
    Ts->Fill();
  }
  fs->Write();
  return 0;
} 
  
double Lsample::FindChi2Limit(const char * paraset, const double CL, const int cut){
  TFile * fr = new TFile(paraset, "r");
  TTree * Tr = (TTree *) fr->Get("fitpara");
  TH1D * h0 = new TH1D("h0", "h0", 1, 0.0, 10000.0);
  TH1D * h1 = new TH1D("h1", "h1", 1000, 0.0, 20.0);
  if (cut == 0){
    Tr->Project("h0", "Chi2", "weight*(1.>0)");
    Tr->Project("h1", "Chi2", "weight*(1.>0)");
  }
  else {
    Tr->Project("h0", "Chi2", "weight*(bound>0.1)");
    Tr->Project("h1", "Chi2", "weight*(bound>0.1)");
  }
  double sum = h0->Integral(1, -1);
  double upper = 0;
  for (int i = 1; i <= 1000; i++){
    if (h1->Integral(1, i) / sum < CL) continue;
    upper = h1->GetBinCenter(i);
    break;
  }
  return upper;
}

#endif
