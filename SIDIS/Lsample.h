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
  double _para0[11];
 public:
  Lsample(int Atype);
  int SetParameters(const double * para);
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


