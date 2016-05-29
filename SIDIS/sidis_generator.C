#include "Lgenerator.h"

using namespace std;
using namespace TMath;

int main (int argc, char* argv[]){
  double lab[7];
  double phys[9];
  double jac;
  double accele[3];
  double accpion[2];
  double Nsim;

  Lgenerator gene(1);
  gene.SetRange(argv[1]);
  gene.GetNsim(&Nsim);
  
  TFile * fdata = new TFile("sidis_data.root","RECREATE");
  TTree * Tdata = new TTree("T","sidis_data_tree");
  Tdata->SetDirectory(fdata);

  Tdata->Branch("Ebeam", &lab[0], "Ebeam/D");
  Tdata->Branch("p_ele", &lab[1], "p_ele/D");
  Tdata->Branch("theta_ele", &lab[2], "theta_ele/D");
  Tdata->Branch("phi_ele", &lab[3], "phi_ele/D");
  Tdata->Branch("p_pion", &lab[4], "p_pion/D");
  Tdata->Branch("theta_pion", &lab[5], "theta_pion/D");
  Tdata->Branch("phi_pion", &lab[6], "phi_pion/D");
  Tdata->Branch("x", &phys[0], "x/D");
  Tdata->Branch("y", &phys[1], "y/D");
  Tdata->Branch("z", &phys[2], "z/D");
  Tdata->Branch("Q2", &phys[3], "Q2/D");
  Tdata->Branch("Pt", &phys[4], "Pt/D");
  Tdata->Branch("phi_h", &phys[5], "phi_h/D");
  Tdata->Branch("phi_S", &phys[6], "phi_S/D");
  Tdata->Branch("W", &phys[7], "W/D");
  Tdata->Branch("Wp", &phys[8], "Wp/D");
  Tdata->Branch("jac", &jac, "jac/D");
  Tdata->Branch("acc_ele", &accele[0], "acc_ele/D");
  Tdata->Branch("acc_ele_FA", &accele[1], "acc_ele_FA/D");
  Tdata->Branch("acc_ele_LA", &accele[2], "acc_ele_LA/D");
  Tdata->Branch("acc_pion_p", &accpion[0], "acc_pion_p/D");
  Tdata->Branch("acc_pion_m", &accpion[1], "acc_pion_m/D");

  Long64_t nsim = 0;

  gRandom->SetSeed(0);
  for (nsim = 0; nsim < Nsim; nsim = nsim + 1){
    if (nsim%10000000 == 0) cout << nsim << endl;
    gene.GenerateEvent(lab);
    gene.CalcVariables(lab, phys);

    if (gene.PreCut(phys)) continue;

    gene.GetAcc(lab, accele, accpion);
    if (lab[1] < 3.5){
      accele[0] = accele[1];//large angle low momentum cut
    }
    if (accele[0] == 0) continue;//acc_ele
    if (accpion[0] == 0 && accpion[1] == 0) continue;//acc_pion
 
    jac = gene.Jacobian(lab);
   
    Tdata->Fill();
  }
  fdata->Write();

  return 0;
}
