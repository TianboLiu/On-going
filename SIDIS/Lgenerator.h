#ifndef _LGENERATOR_H_
#define _LGENERATOR_H_

#include <iostream>
#include <fstream>
#include <cmath>

#include "TROOT.h"
#include "TApplication.h"
#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TH3F.h"
#include "TEventList.h"

class Lgenerator{
 protected:
  int _TARGET;// 0: 3He, 1: NH3
  double run[11];
  TFile * file_electron;
  TFile * file_pionp;
  TFile * file_pionm;
  void * accept_ele_FA;
  void * accept_ele_LA;
  void * accept_pion_plus;
  void * accept_pion_minus;
 public:
  Lgenerator(int tar);
  int Print();
  int SetRange(const char file[]);
  int GetNsim(double * Nsim);
  int GetAccElectron(const double * elab, double * accele);
  int GetAccPion(const double * pilab, double * accpion);
  int GetAcc(const double * lab, double * accele, double * accpion);
  int GenerateEvent(double * lab);
  int CalcVariables(const double * lab, double * phys);
  double Jacobian(const double * lab);
  int PreCut(const double * phys);
};

Lgenerator::Lgenerator(int tar){
  _TARGET = tar;
  if (tar != 0 && tar != 1)
    std::cout << "Error: Lgenerator::Lgenerator: Invalid initialization!" << std::endl;
  if (tar == 0){
    file_electron = new TFile("/var/phy/project/mepg/tl190/SoLIDacceptance/acceptance_solid_CLEO_SIDIS_3he_negative_output.root","r");
    file_pionp = new TFile("/var/phy/project/mepg/tl190/SoLIDacceptance/acceptance_solid_CLEO_SIDIS_3he_positive_output.root","r");
    file_pionm = new TFile("/var/phy/project/mepg/tl190/SoLIDacceptance/acceptance_solid_CLEO_SIDIS_3he_negative_output.root","r");
    accept_ele_FA = (TH2F*) file_electron->Get("acceptance_forwardangle");
    accept_ele_LA = (TH2F*) file_electron->Get("acceptance_largeangle");
    accept_pion_plus = (TH2F*) file_pionp->Get("acceptance_forwardangle");
    accept_pion_minus = (TH2F*) file_pionm->Get("acceptance_forwardangle");
  }
  else if (tar == 1){
    file_electron = new TFile("/var/phy/project/mepg/tl190/SoLIDacceptance/acceptance_solid_SIDIS_NH3_electron_output.root","r");
    file_pionp = new TFile("/var/phy/project/mepg/tl190/SoLIDacceptance/acceptance_solid_SIDIS_NH3_pionp_output.root","r");
    file_pionm = new TFile("/var/phy/project/mepg/tl190/SoLIDacceptance/acceptance_solid_SIDIS_NH3_pionm_output.root","r");
    accept_ele_FA = (TH3F*) file_electron->Get("acceptance_ThetaPhiP_forwardangle");
    accept_ele_LA = (TH3F*) file_electron->Get("acceptance_ThetaPhiP_largeangle");
    accept_pion_plus = (TH3F*) file_pionp->Get("acceptance_ThetaPhiP_forwardangle");
    accept_pion_minus = (TH3F*) file_pionm->Get("acceptance_ThetaPhiP_forwardangle");
  }
}

int Lgenerator::Print(){
  if (_TARGET == 0){
    std::cout << "Configuration: 3He" << std::endl;
  }
  else if (_TARGET == 1){
    std::cout << "Configuration: NH3" << std::endl;
  }
  else {
    std::cout << "Error: Invalid configuration!" << std::endl;
    return 1;
  }
  std::cout << "Nsim: " << run[0] << std::endl;
  std::cout << "p_ele_min: " << run[1] << std::endl;
  std::cout << "p_ele_max: " << run[2] << std::endl;
  std::cout << "theta_ele_min: " << run[3] << std::endl;
  std::cout << "theta_ele_max: " << run[4] << std::endl;
  std::cout << "p_pion_min: " << run[5] << std::endl;
  std::cout << "p_pion_max: " << run[6] << std::endl;
  std::cout << "theta_pion_min: " << run[7] << std::endl;
  std::cout << "theta_pion_max: " << run[8] << std::endl;
  std::cout << "Ebeam: " << run[9] << std::endl;
  std::cout << "genvol: " << run[10] << std::endl;
  return 0;
}

int Lgenerator::SetRange(const char file[] = "runfile.dat"){
  const double degtorad = M_PI / 180.0;
  const double twopi = 2.0 * M_PI;
  TString filename = file;
  ifstream infile(filename);
  if (!infile.is_open()){
    std::cout << "Lgenerator::setrange: File does not exist!" << std::endl;
    return 1;
  }
  std::cout << "Run file: " << filename << std::endl;
  TString temp;
  int i = 0;
  while (infile >> temp >> run[i]){
    i++;
  }
  infile.close();
  double p_ele_min = run[1];
  double p_ele_max = run[2];
  double theta_ele_min = run[3];
  double theta_ele_max = run[4];
  double p_pion_min = run[5];
  double p_pion_max = run[6];
  double theta_pion_min = run[7];
  double theta_pion_max = run[8];
  double genvol = (p_ele_max - p_ele_min) * (cos(theta_ele_min * degtorad) - cos(theta_ele_max * degtorad)) * (p_pion_max - p_pion_min) * (cos(theta_pion_min * degtorad) - cos(theta_pion_max * degtorad)) * twopi * twopi;
  run[10] = genvol;
  Print();
  return 0;
}

int Lgenerator::GetNsim(double * Nsim){
  Nsim[0] = run[0];
  return 0;
}

int Lgenerator::GetAccElectron(const double * elab, double * accele){
  //elab: p_ele, theta_ele, phi_ele
  //accele: total, forward, large
  if (_TARGET == 0){
    accele[1] = ((TH2F *) accept_ele_FA)->GetBinContent(((TH2F *) accept_ele_FA)->GetXaxis()->FindBin(elab[1]), ((TH2F *) accept_ele_FA)->GetYaxis()->FindBin(elab[0]));
    accele[2] = ((TH2F *) accept_ele_LA)->GetBinContent(((TH2F *) accept_ele_LA)->GetXaxis()->FindBin(elab[1]), ((TH2F *) accept_ele_LA)->GetYaxis()->FindBin(elab[0]));
    accele[0] = accele[1] + accele[2];
  }
  else if (_TARGET == 1){
    accele[1] = ((TH3F *) accept_ele_FA)->GetBinContent(((TH3F *) accept_ele_FA)->GetXaxis()->FindBin(elab[1]), ((TH3F *) accept_ele_FA)->GetYaxis()->FindBin(elab[2]), ((TH3F *) accept_ele_FA)->GetZaxis()->FindBin(elab[0]));
    accele[2] = ((TH3F *) accept_ele_LA)->GetBinContent(((TH3F *) accept_ele_LA)->GetXaxis()->FindBin(elab[1]), ((TH3F *) accept_ele_LA)->GetYaxis()->FindBin(elab[2]), ((TH3F *) accept_ele_LA)->GetZaxis()->FindBin(elab[0]));
    accele[0] = accele[1] + accele[2];
  }
  else {
    return 1;
  }
  return 0;
}

int Lgenerator::GetAccPion(const double * pilab, double * accpion){
  //pilab: p_pion, theta_pion, phi_pion
  //accpion: pi+, pi-
  if (_TARGET == 0){
    accpion[0] = ((TH2F *) accept_pion_plus)->GetBinContent(((TH2F *) accept_pion_plus)->GetXaxis()->FindBin(pilab[1]), ((TH2F *) accept_pion_plus)->GetYaxis()->FindBin(pilab[0]));
    accpion[1] = ((TH2F *) accept_pion_minus)->GetBinContent(((TH2F *) accept_pion_minus)->GetXaxis()->FindBin(pilab[1]), ((TH2F *) accept_pion_minus)->GetYaxis()->FindBin(pilab[0]));
  }
  else if (_TARGET == 1){
    accpion[0] = ((TH3F *) accept_pion_plus)->GetBinContent(((TH3F *) accept_pion_plus)->GetXaxis()->FindBin(pilab[1]), ((TH3F *) accept_pion_plus)->GetYaxis()->FindBin(pilab[2]), ((TH3F *) accept_pion_plus)->GetZaxis()->FindBin(pilab[0]));
    accpion[1] = ((TH3F *) accept_pion_minus)->GetBinContent(((TH3F *) accept_pion_minus)->GetXaxis()->FindBin(pilab[1]), ((TH3F *) accept_pion_minus)->GetYaxis()->FindBin(pilab[2]), ((TH3F *) accept_pion_minus)->GetZaxis()->FindBin(pilab[0]));
  }
  else {
    return 1;
  }
  return 0;
}

int Lgenerator::GetAcc(const double * lab, double * accele, double * accpion){
  //lab: Ebeam, p_ele, theta_ele, phi_ele, p_pion, theta_pion, phi_pion
  GetAccElectron(&lab[1], accele);
  GetAccPion(&lab[4], accpion);
  return 0;
}

int Lgenerator::GenerateEvent(double * lab){
  //run: Nsim, e_pl, pu, thl, thu, pi_pl, pu, thl, thu, Ebeam, genvol
  //lab: Ebeam, p_ele, theta_ele, phi_ele, p_pion, theta_pion, phi_pion
  const double degtorad = M_PI / 180.0;
  const double radtodeg = 180.0 / M_PI;
  lab[0] = run[9];
  lab[1] = gRandom->Uniform(run[1], run[2]);
  lab[2] = acos(gRandom->Uniform(cos(run[4]*degtorad), cos(run[3]*degtorad))) * radtodeg;
  lab[3] = gRandom->Uniform(-180.0, 180.0);
  lab[4] = gRandom->Uniform(run[5], run[6]);
  lab[5] = acos(gRandom->Uniform(cos(run[8]*degtorad), cos(run[7]*degtorad))) * radtodeg;
  lab[6] = gRandom->Uniform(-180.0, 180.0);
  return 0;
}

int Lgenerator::CalcVariables(const double * lab, double * phys){
  //lab: Ebeam, p_ele, theta_ele, phi_ele, p_pion, theta_pion, phi_pion
  const double degtorad = M_PI / 180.0;
  double Ebeam = lab[0];
  double p_ele = lab[1];
  double theta_ele = lab[2] * degtorad;
  double phi_ele = lab[3] * degtorad;
  double p_pion = lab[4];
  double theta_pion = lab[5] * degtorad;
  double phi_pion = lab[6] * degtorad;
  const double Mnucleon = 0.939;
  const double Mpion = 0.13957;
  double q0 = Ebeam - p_ele;
  double q1 = - p_ele * sin(theta_ele);
  double q3 = Ebeam - p_ele * cos(theta_ele);
  double qt = sqrt(q1*q1 + q3*q3);
  double Q2 = 2.0 * Ebeam * p_ele * (1.0 - cos(theta_ele));
  double x = Q2 / (2.0 * Mnucleon * q0);
  double y = q0 / Ebeam;
  double Epion = sqrt(Mpion * Mpion + p_pion * p_pion);
  double ppion1 = p_pion * sin(theta_pion) * cos(phi_pion - phi_ele);
  double ppion2 = p_pion * sin(theta_pion) * sin(phi_pion - phi_ele);
  double ppion3 = p_pion * cos(theta_pion);
  double z = Epion / q0;
  double W2 = Mnucleon * Mnucleon + (1.0 - x) / x * Q2;
  double Wp2 = W2 + Mpion * Mpion - 2.0 * ((Mnucleon + q0) * Epion - q1 * ppion1 - q3 * ppion3);
  double ct = q3 / qt;
  double st = q1 / qt;
  //double ppionz = ppion3 * ct + ppion1 * st;
  double ppionx = -ppion3 * st + ppion1 * ct;
  double Pt = sqrt(ppionx * ppionx + ppion2 * ppion2);
  double phi_h;
  if (ppion2 >= 0) phi_h = acos(ppionx / Pt);
  else phi_h = -acos(ppionx / Pt);
  double phi_S = - phi_ele;
  //phys: x, y, z, Q2, Pt, phih, phiS, W, Wp
  phys[0] = x;
  phys[1] = y;
  phys[2] = z;
  phys[3] = Q2;
  phys[4] = Pt;
  phys[5] = phi_h;
  phys[6] = phi_S;
  phys[7] = sqrt(W2);
  phys[8] = sqrt(Wp2);
  return 0;
}

double Lgenerator::Jacobian(const double * lab){
  //lab: Ebeam, p_ele, theta_ele, phi_ele, p_pion, theta_pion, phi_pion
  const double degtorad = M_PI / 180.0;
  const double Mnucleon = 0.939;
  const double Mpion = 0.13957;
  double Ebeam = lab[0];
  double p_ele = lab[1];
  double theta_ele = lab[2] * degtorad;
  double phi_ele = lab[3] * degtorad;
  double p_pion = lab[4];
  double theta_pion = lab[5] * degtorad;
  double phi_pion = lab[6] * degtorad;
  double p3 = Ebeam - p_ele * cos(theta_ele);
  double pT = p_pion * sqrt(1.0 + (-pow(p3 * cos(theta_pion), 2) + p3 * p_ele * sin(theta_ele) * sin(2.0 * theta_pion) * cos(phi_ele - phi_pion) - pow(p_ele * sin(theta_ele) * sin(theta_pion) * cos(phi_ele - phi_pion), 2)) / (pow(p3, 2) + pow(p_ele * sin(theta_ele), 2)) );
  double factor = (p_pion * p_ele * sin(theta_ele)) / (Mnucleon * pow(Ebeam - p_ele, 2) * sqrt(pow(p_pion, 2) + pow(Mpion, 2)));
  double jnum = 4.0 * p_pion * sin(theta_pion) * (p_ele * sin(theta_ele) * sin(theta_pion) * cos(phi_ele - phi_pion) - p3 * cos(theta_pion));
  double jden = sqrt(std::abs( -8.0 * pow(p3, 2) * cos(2.0 * theta_pion) + p_ele * (16.0 * p3 * sin(theta_ele) * sin(2.0 * theta_pion) * cos(phi_ele - phi_pion) - p_ele * (8.0 * pow(sin(theta_ele) * sin(theta_pion), 2) * cos(2.0 * (phi_ele - phi_pion)) + cos(2.0 * (theta_ele - theta_pion)) + cos(2.0 * (theta_ele + theta_pion)) + 6.0 * cos(2.0 * theta_ele) - 2.0 * cos(2.0 * theta_pion))) + 8.0 * pow(p3, 2) + 6.0 * pow(p_ele, 2)));
  return std::abs((2.0 * pT * factor * jnum) / (jden * sin(theta_ele) * sin(theta_pion)));
}

int Lgenerator::PreCut(const double * phys){
  //x
  if (phys[0] < 0.01 || phys[0] > 1.0) return 1;
  //y
  if (phys[1] < 0.0 || phys[1] > 1.0) return 1;
  //z
  if (phys[2] < 0.3 || phys[2] > 0.7) return 1;
  //Q2
  if (phys[3] < 1.0 || phys[3] > 10.0) return 1;
  //Pt
  if (phys[4] < 0.0 || phys[4] > 1.6) return 1;
  //W
  if (phys[7] < 2.3) return 1;
  //Wp
  if (phys[8] < 1.6) return 1;
  //
  return 0;
}


#endif
