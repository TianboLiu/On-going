#ifndef _DIS_H_
#define _DIS_H_

#include <iostream>
#include <fstream>
#include <cmath>

#include "stdlib.h"
#include "TROOT.h"
#include "TApplication.h"
#include "Rtypes.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TObject.h"
#include "TKey.h"
#include "TTree.h"
#include "TString.h"
#include "TEventList.h"

#include "Lstructure.h"

const LHAPDF::PDF * pdf = _pdfs;

extern "C"{
  void f1f2in09_(double * Z, double * A, double * Q2, double * W2, double * F1, double * F2, double * rc);
  void f1f2qe09_(double * Z, double * A, double * Q2, double * W2, double * F1, double * F2);
}

double calin09(double E, double pf, double thf, int Np = 2, int Nn = 1){
  const double alpha = 1.0 / 137.0;
  const double M = 0.938272;
  const double pi = M_PI;
  double Z = (double) Np;
  double A = ( (double) Np) + ( (double) Nn);
  double Q2 = 4.0 * E * pf * sin(thf/2.0) * sin(thf/2.0);
  double W2 = M*M + 2.0*M*(E-pf) - Q2;
  // cout << "W2: " << W2 << endl;
  double F1, F2, rc;
  f1f2in09_(&Z, &A, &Q2, &W2, &F1, &F2, &rc);
  double y = (E-pf) / E;
  double x = Q2 / (2.0 * M * (E - pf));
  
  double sigma = (4.0*pi*alpha*alpha)/(x*y*Q2) * ( x*y*y * F1 + (1.0 - y - pow(x*y*M,2)/Q2) * F2);

  return sigma  * pf / (2.0*pi*M*(E-pf));
}

double calqe09(double E, double pf, double thf, int Np = 2, int Nn = 1){
  const double alpha = 1.0 / 137.0;
  const double M = 0.938272;
  const double pi = M_PI;
  double Z = (double) Np;
  double A = ( (double) Np) + ( (double) Nn);
  double Q2 = 4.0 * E * pf * sin(thf/2.0) * sin(thf/2.0);
  double W2 = M*M + 2.0*M*(E-pf) - Q2;
  double F1, F2;
  f1f2qe09_(&Z, &A, &Q2, &W2, &F1, &F2);
  double y = (E-pf) / E;
  double x = Q2 / (2.0 * M * (E - pf));
  
  double sigma = (4.0*pi*alpha*alpha)/(x*y*Q2) * ( x*y*y * F1 + (1.0 - y - pow(x*y*M,2)/Q2) * F2);

  return sigma * pf / (2.0*pi*M*(E-pf));
}

double cal09(double E, double pf, double thf, int Np = 2, int Nn = 1){
  return calin09(E, pf, thf, Np, Nn) + calqe09(E, pf, thf, Np, Nn);
}

// DIS
double Alpha_S(double Q2){
  return pdf->alphasQ2(Q2);
}

double xpdf_P(int flavor, double x, double Q2){
  return pdf->xfxQ2(flavor, x, Q2);
}

double xpdf_N(int flavor, double x, double Q2){
  if (abs(flavor) == 1 || abs(flavor) == 2){
    return pdf->xfxQ2(2/flavor, x, Q2);
  }
  else {
    return pdf->xfxQ2(flavor, x, Q2);
  }
}

double dis_F2_P(double x, double Q2){
  const double eu = 2.0/3.0;
  const double ed = -1.0/3.0;
  double F2 = eu*eu*( xpdf_P(2,x,Q2) + xpdf_P(-2,x,Q2) + xpdf_P(4,x,Q2) + xpdf_P(-4,x,Q2))
    + ed*ed*( xpdf_P(1,x,Q2) + xpdf_P(-1,x,Q2) + xpdf_P(3,x,Q2) + xpdf_P(-3,x,Q2) + xpdf_P(5,x,Q2) + xpdf_P(-5,x,Q2));
  return F2;
}

double dis_F2_N(double x, double Q2){
  const double eu = 2.0/3.0;
  const double ed = -1.0/3.0;
  double F2 = eu*eu*( xpdf_N(2,x,Q2) + xpdf_N(-2,x,Q2) + xpdf_N(4,x,Q2) + xpdf_N(-4,x,Q2))
    + ed*ed*( xpdf_N(1,x,Q2) + xpdf_N(-1,x,Q2) + xpdf_N(3,x,Q2) + xpdf_N(-3,x,Q2) + xpdf_N(5,x,Q2) + xpdf_N(-5,x,Q2));
  return F2;
}

double dis_sigma_P(double E, double Q2, double W2){
  //d sigma/ dp_e dOmega_e in unit GeV^-2 / GeV.str
  const double alpha_em = 1.0 / 137.0;
  const double M = 0.938272;
  const double pi = M_PI;
  
  double x = Q2 / (Q2 + W2 - M*M);
  double y = Q2 / (2.0 * x * M * E);
  if (x > 1.0 || y > 1.0) { return 0.0;}
  double F2 = dis_F2_P(x, Q2);
  double F1 = F2 / (2.0*x);

  double sigma = (4.0*pi*alpha_em*alpha_em) / (x*y*Q2) * ((1.0 - y - pow(x*y*M,2)/Q2) * F2 + x*y*y*F1);
  return sigma * (1.0 - y) / (2.0 * pi * M * y);
}

double dis_sigma_N(double E, double Q2, double W2){
  //d sigma/ dp_e dOmega_e in unit GeV^-2 / GeV.str
  const double alpha_em = 1.0 / 137.0;
  const double M = 0.9395654;
  const double pi = M_PI;
  
  double x = Q2 / (Q2 + W2 - M*M);
  double y = Q2 / (2.0 * x * M * E);
  if (x > 1.0 || y > 1.0) { return 0.0;}
  double F2 = dis_F2_N(x, Q2);
  double F1 = F2 / (2.0*x);

  double sigma = (4.0*pi*alpha_em*alpha_em) / (x*y*Q2) * ((1.0 - y - pow(x*y*M,2)/Q2) * F2 + x*y*y*F1);
  return sigma * (1.0 - y) / (2.0 * pi * M * y);
}

double dis_sigma(double E, double Q2, double W2, int Np = 2, int Nn = 1){
  double sigma = Np * dis_sigma_P(E, Q2, W2) + Nn * dis_sigma_N(E, Q2, W2);
  return sigma;
}

double calDIS(double E, double p, double th, int Np = 2, int Nn = 1){
  double Q2 = 2.0 * E * p * (1.0 - cos(th));
  const double Mp = 0.938272;
  const double Mn = 0.9395654;
  double W2p = Mp*Mp + 2.0*Mp*(E - p) - Q2;
  double W2n = Mn*Mn + 2.0*Mn*(E - p) - Q2;
  double sigma = Np * dis_sigma_P(E, Q2, W2p) + Nn * dis_sigma_N(E, Q2, W2n);
  return sigma;
}

//combine
double calcombine(double E, double pf, double thf, int Np = 2, int Nn = 1){
  const double M = 0.938272;
  double Q2 = 4.0 * E * pf * sin(thf/2.0) * sin(thf/2.0);
  double W2 = M*M + 2.0*M*(E-pf) - Q2;
  if (W2 < 9.0){
    return calin09(E, pf, thf, Np, Nn);
  }
  else if (W2 >= 9.0){
    return calDIS(E, pf, thf, Np, Nn);
  }
  else {return 0;}
}

#endif
