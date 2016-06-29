#ifndef _LSAMPLE_H_
#define _LSAMPLE_H_

#include <iostream>
#include <fstream>
#include <cmath>
#include "TRandom.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TString.h"
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
  double _eigenvectors[11][11];
  double _eigenvalues[11];
  double _para[11];
  double _central[11];
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
  for (int i = 0; i < ndata; i++){
    Td->GetEntry(i);
    _Nucleon[_Ndata + i] = Nucleon;
    _Hadron[_Ndata + i] = Hadron;
    for (int j = 0; j < 5; j++){
      _kin[_Ndata + i][j] = kin[j];
    }
    _A0[_Ndata + i] = Asym[Atype - 1];
    _E0[_Ndata + i] = Estat[Atype - 1];
    _ET[_Ndata + i] = sqrt(Estat[Atype - 1] * Estat[Atype - 1] + Eabs[Atype - 1] * Eabs[Atype - 1] + Erel[Atype - 1] * Erel[Atype - 1] * Asym[Atype - 1] * Asym[Atype - 1]);
  }
  _Ndata = _Ndata + ndata;
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
  


#endif
