#ifndef _SIVERS_H_
#define _SIVERS_H_

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

//collinear PDF (CTEQ6L)
const LHAPDF::PDF * _pdfs = LHAPDF::mkPDF("cteq6l1", 0);
int f1_cteq6l(const double * var, double * f1){
  //var: x, Q2
  double x = var[0];
  double Q2 = var[1];
  //f1: u, d, s, ubar, dbar, sbar
  f1[0] = _pdfs->xfxQ2(2, x, Q2) / x;//u
  f1[1] = _pdfs->xfxQ2(1, x, Q2) / x;//d
  f1[2] = _pdfs->xfxQ2(3, x, Q2) / x;//s
  f1[3] = _pdfs->xfxQ2(-2, x, Q2) / x;//ubar
  f1[4] = _pdfs->xfxQ2(-1, x, Q2) / x;//dbar
  f1[5] = _pdfs->xfxQ2(-3, x, Q2) / x;//sbar
  return 0;
}

int f1p(const double * var, double * f1){
  f1_cteq6l(var, f1);
  return 0;
}

int f1n(const double * var, double * f1){
  double f1p[6];
  f1_cteq6l(var, f1p);
  f1[0] = f1p[1];
  f1[1] = f1p[0];
  f1[2] = f1p[2];
  f1[3] = f1p[4];
  f1[4] = f1p[3];
  f1[5] = f1p[5];
  return 0;
}

int f1col(const double * var, double * f1, const char * tar = "proton"){
  if (strcmp(tar, "proton") == 0)
    f1p(var, f1);
  else if (strcmp(tar, "neutron") == 0)
    f1n(var, f1);
  else if (strcmp(tar, "deuteron") == 0){
    double fp[6], fn[6];
    f1p(var, fp);
    f1n(var, fn);
    for (int i = 0; i < 6; i++)
      f1[i] = fp[i] + fn[i];
  }
  else
    std::cerr << "Fail to match target option!" << std::endl;
  return 0;
}

//collinear fragmentation (DSS)
extern "C"{
  extern struct { int FINI;} fraginid_;
}
extern "C"{
  void fdss_(int * IH, int * IC, int * IO, double * Z, double * Q2, double * zU, double * zUBAR, double * zD, double * zDBAR, double * zS, double * zSBAR, double * zC, double * zB, double * zG);
}
int D1col(const double * var, double * D1, const char * had = "pip"){
  //var: z, Q2
  double z = var[0];
  double Q2 = var[1];
  int io = 0;
  int ih = 1;
  int ic = 1;
  if (strcmp(had, "pip") == 0){
    ih = 1; ic = 1;
  }
  else if (strcmp(had, "pim") == 0){
    ih = 1; ic = -1;
  }
  else if (strcmp(had, "pi0") == 0){
    ih = 1; ic = 0;
  }
  else if (strcmp(had, "Kp") == 0){
    ih = 2; ic = 1;
  }
  else if (strcmp(had, "Km") == 0){
    ih = 2; ic = -1;
  }
  else if (strcmp(had, "K0") == 0){
    ih = 2; ic = 0;
  }
  else if (strcmp(had, "proton") == 0){
    ih = 3; ic = 1;
  }
   else if (strcmp(had, "hp") == 0){
    ih = 4; ic = 1;
  }
  else if (strcmp(had, "hm") == 0){
    ih = 4; ic = -1;
  }
  double u, ubar, d, dbar, s, sbar, c, b, g;
  fdss_ (&ih, &ic, &io, &z, &Q2, &u, &ubar, &d, &dbar, &s, &sbar, &c, &b, &g);
  D1[0] = u / z;//u
  D1[1] = d / z;//d
  D1[2] = s / z;//s
  D1[3] = ubar / z;//ubar
  D1[4] = dbar / z;//dbar
  D1[5] = sbar / z;//sbar
  return 0;
}

//Sivers function
int f1tM0p(const double * var, double * f1t, const double * para = 0){
  //var: x, Q2
  double x = var[0];
  //parameter order below
  double Nuv, Ndv, Nubar, Ndbar, auv, adv, buv, bdv, M2;
  if (para == 0){
    Nuv = 0.180119;
    Ndv = -0.520384;
    Nubar = -0.0105975;
    Ndbar = -0.0587188;
    auv = 1.02360;
    adv = 1.91783;
    buv = 6.60053;
    bdv = 10.0849;
    M2 = 0.800441;
  }
  else {
    Nuv = para[0];
    Ndv = para[1];
    Nubar = para[2];
    Ndbar = para[3];
    auv = para[4];
    adv = para[5];
    buv = para[6];
    bdv = para[7];
    M2 = para[8];
  }
  double f1[6];
  f1_cteq6l(var, f1);
  const double Mp = 0.93827;
  const double kt2 = 0.57;
  double M1 = sqrt(M2);
  f1t[0] = - Nuv * pow(x, auv) * pow(1.0-x, buv) * pow(auv+buv, auv+buv) / pow(auv, auv) / pow(buv, buv) * 0.5 * sqrt(2.0*exp(1)) * M1 * Mp / (M2 + kt2) * (f1[0] - f1[3]);
  f1t[1] = - Ndv * pow(x, adv) * pow(1.0-x, bdv) * pow(adv+bdv, adv+bdv) / pow(adv, adv) / pow(bdv, bdv) * 0.5 * sqrt(2.0*exp(1)) * M1 * Mp / (M2 + kt2) * (f1[1] - f1[4]);
  f1t[2] = 0;
  f1t[3] = - Nubar * 0.5 * sqrt(2.0*exp(1)) * M1 * Mp / (M2 + kt2) * f1[3];
  f1t[4] = - Ndbar * 0.5 * sqrt(2.0*exp(1)) * M1 * Mp / (M2 + kt2) * f1[4];
  f1t[5] = 0;
  f1t[0] = f1t[0] + f1t[3];
  f1t[1] = f1t[1] + f1t[4];
  return 0;
}

int f1tM0n(const double * var, double * f1t, const double * para = 0){
  double f1tp[6];
  f1tM0p(var, f1tp, para);
  f1t[0] = f1tp[1];
  f1t[1] = f1tp[0];
  f1t[2] = f1tp[2];
  f1t[3] = f1tp[4];
  f1t[4] = f1tp[3];
  f1t[5] = f1tp[5];
  return 0;
}

int f1tM0(const double * var, double * f1t, const char * tar = "proton", const double * para = 0){
  if (strcmp(tar, "proton") == 0)
    f1tM0p(var, f1t, para);
  else if (strcmp(tar, "neutron") == 0)
    f1tM0n(var, f1t, para);
  else if (strcmp(tar, "deuteron") == 0){
    double f1tp[6], f1tn[6];
    f1tM0p(var, f1tp, para);
    f1tM0n(var, f1tn, para);
    for (int i = 0; i < 6; i++)
      f1t[i] = f1tp[i] + f1tn[i];
  }
  else
    std::cerr << "Fail to match target option!" << std::endl;
  return 0;
}

int f1siv(const double * var, double * siv, const double * para = 0){
  //var: x, Q2, kt
  double kt = var[2];
  double M2;
  if (para == 0){
    M2 = 0.800441;
  }
  else{
    M2 = para[8];
  }
  double ks2 = 0.57 * M2 / (0.57 + M2);
  double f1t[6];
  f1tM0p(var, f1t, para);
  for (int i = 0; i < 6; i++)
    siv[i] = f1t[i] * exp(-kt*kt/ks2) / (M_PI * ks2);
  return 0;
}

double FUU(const double * kin, const char * tar = "proton", const char * had = "pip"){
  //kin: x, Q2, z, Pt
  double x = kin[0];
  double Q2 = kin[1];
  double z = kin[2];
  double Pt = kin[3];
  double f1[6], D1[6];
  double vf[2] = {x, Q2};
  double vD[2] = {z, Q2};
  f1col(vf, f1, tar);
  D1col(vD, D1, had);
  double kt2 = 0.57;
  double pt2 = 0.12;
  double Ph2 = pt2 + z * z * kt2;
  double result = (4.0/9.0) * (f1[0] * D1[0] + f1[3] * D1[3]) + (1.0/9.0) * (f1[1] * D1[1] + f1[2] * D1[2] + f1[4] * D1[4] + f1[5] * D1[5]);
  return result * x * exp(-Pt*Pt/Ph2) / (M_PI * Ph2);
}

double FUT(const double * kin, const char * tar = "proton", const char * had = "pip", const double * para = 0){
  //kin: x, Q2, z, Pt
  double x = kin[0];
  double Q2 = kin[1];
  double z = kin[2];
  double Pt = kin[3];
  double M2;
  if (para == 0){
    M2 = 0.800441;
  }
  else{
    M2 = para[8];
  }
  double ks2 = 0.57 * M2 / (0.57 + M2);
  double pt2 = 0.12;
  double Ph2 = pt2 + z * z * ks2;
  double f1t[6], D1[6];
  double vf[2] = {x, Q2};
  double vD[2] = {z, Q2};
  f1tM0(vf, f1t, tar, para);
  D1col(vD, D1, had);
  const double Mp = 0.93827;
  double result = (4.0/9.0) * (f1t[0] * D1[0] + f1t[3] * D1[3]) + (1.0/9.0) * (f1t[1] * D1[1] + f1t[2] * D1[2] + f1t[4] * D1[4] + f1t[5] * D1[5]);
  return result * x * exp(-Pt*Pt/Ph2) / (M_PI * Ph2) * (-z * ks2 * Pt / Mp / Ph2);
}

double AUT(const double * kin, const char * tar = "proton", const char * had = "pip", const double * para = 0){
  double unpol = FUU(kin, tar, had);
  double pol = FUT(kin, tar, had, para);
  return pol / unpol;
}



#endif
