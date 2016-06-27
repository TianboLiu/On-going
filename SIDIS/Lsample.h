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
#include "TTree.h"
#include "TFile.h"
#include "Lstructure.h"

#define _NMAX_ 4000

double _par1[11] = {0.35, -0.90, -0.24, 0.04, -0.40, 1.0, 0.73, 1.08, 0.79, 3.46, 0.34};
double _par2[10] = {0.46, -1.0, 1.11, 3.64, 0.25, 0.49, -1.0, 1.06, 0.07, 1.5};
double _par3[10] = {1.0, -1.0, 2.5, 2.0, 0.18, 0.49, -1.0, 1.06, 0.07, 1.5};

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
  double _para0[11];
  double _central[11];
 public:
  Lsample(int Atype);
  int InitialA1para();
  int InitialA2para();
  int InitialA3para();
  int GetNeutronPlusData(const char * datafile);
  int GetNeutronMinusData(const char * datafile);
  int GetProtonPlusData(const char * datafile);
  int GetProtonMinusData(const char * datafile);
};

Lsample::Lsample(int Atype){
  _Atype = Atype;
  if (Atype == 1) _PSD = 9;
  else if (Atype == 2) _PSD = 9;
  else if (Atype == 3) _PSD = 4;
  else {
    _PSD = 0;
    std::cout << "Lsample initialize error: invalid asymmetry type!" << std::endl;
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
  _para[4] = 0.25;//kt2
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
  _central[3] = _para[3];//b
  _central[4] = _para[5];//Nfav
  _central[5] = _para[6];//Ndis
  _central[6] = _para[7];//c
  _central[7] = _para[8];//d
  _central[8] = _para[9];//MC2
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

