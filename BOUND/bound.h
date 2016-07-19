#ifndef _BOUND_H_
#define _BOUND_H_

#include <iostream>
#include <fstream>
#include <cmath>

/*********** Global Variables ***********/
double massN;//Nucleon mass
double massphi;//phi meson mass
double AN;//Nucleon number
/******* End of Global Variables ********/


/************** Functions ***************/

int initialize(){
  //Initialize parameters
  massN = (0.938272046 + 0.939565379) / 2.0;//average p/n mass of PDG value
  massphi = 1.019455;//phi meson mass in GeV
  AN = 12.0;
  return 0;
}

double t0_phiN(double Q, double qc){
  //on-shell amplitude of gamma N -> phi N in c.m. frame in unit GeV^-2
  return 0.005;
}

double t_phiN(double Q, double qc){
  //s-wave transition matrix element of gamma N -> phi N in c.m. frame in unit GeV^-2
  return t0_phiN(Q, qc) / (4.0 * M_PI);
}

double sigmaphiproduction_total(double w){
  //total cross section of gamma N -> phi N in c.m. frame
  //w: total energy
  double qc = (w * w - massN * massN) / (2.0 * w);
  double Q = sqrt(pow((w*w + massphi*massphi - massN*massN)/(2.0 * w), 2) - massphi*massphi);
  double rhoqc = M_PI * qc * qc * sqrt(massN*massN + qc*qc) / (qc + sqrt(massN*massN + qc*qc));
  double rhoQ = M_PI * Q * sqrt(massphi*massphi + Q*Q) * sqrt(massN*massN + Q*Q) / (sqrt(massphi*massphi + Q*Q) + sqrt(massN*massN + Q*Q));
  return 4.0 * M_PI * rhoqc * rhoQ / (qc * qc) * pow(t0_phiN(Q, qc), 2);
}

double Vpotential(double r){
  //van de Waals potential between N and phi
  double a = 1.25;
  double mu = 0.6;//in GeV
  return - a * exp(- mu * r) / r;
}

double Cutoff(double Q, double Lambda){
  //cutoff function on relative momentum of N-phi
  return exp(-Lambda*Lambda*Q*Q);
}







/* References
   [1] H. Gao et al., Phys. Rev. C 63, 022201(R) (2001).
 */

#endif
