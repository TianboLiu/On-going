#include <iostream>
#include <fstream>
using namespace std;

#include "stdlib.h"
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TMath.h"
using namespace TMath;

int main(){

  TFile * fp = new TFile("sidis_data.root");
  TTree * Tp = (TTree *) fp->GetObjectChecked("T", "TTree");

  Long64_t N_entries = Tp->GetEntries();

  Int_t i = 0;

  Double_t x, z, Q2, Pt;
  Tp->SetBranchAddress("x", &x);
  Tp->SetBranchAddress("Q2", &Q2);
  Tp->SetBranchAddress("z", &z);
  Tp->SetBranchAddress("Pt", &Pt);

  // Files and Trees

  // 1 < Q2 < 2
  TFile * fs130 = new TFile("Selected/sidis_select130.root","RECREATE");
  TTree * Ts130 = Tp->CloneTree(0);
  Ts130->SetDirectory(fs130);
  TFile * fs135 = new TFile("Selected/sidis_select135.root","RECREATE");
  TTree * Ts135 = Tp->CloneTree(0);
  Ts135->SetDirectory(fs135);
  TFile * fs140 = new TFile("Selected/sidis_select140.root","RECREATE");
  TTree * Ts140 = Tp->CloneTree(0);
  Ts140->SetDirectory(fs140); 
  TFile * fs145 = new TFile("Selected/sidis_select145.root","RECREATE");
  TTree * Ts145 = Tp->CloneTree(0);
  Ts145->SetDirectory(fs145); 
  TFile * fs150 = new TFile("Selected/sidis_select150.root","RECREATE");
  TTree * Ts150 = Tp->CloneTree(0);
  Ts150->SetDirectory(fs150);
  TFile * fs155 = new TFile("Selected/sidis_select155.root","RECREATE");
  TTree * Ts155 = Tp->CloneTree(0);
  Ts155->SetDirectory(fs155);
  TFile * fs160 = new TFile("Selected/sidis_select160.root","RECREATE");
  TTree * Ts160 = Tp->CloneTree(0);
  Ts160->SetDirectory(fs160);
  TFile * fs165 = new TFile("Selected/sidis_select165.root","RECREATE");
  TTree * Ts165 = Tp->CloneTree(0);
  Ts165->SetDirectory(fs165);
 
   // 2 < Q2 < 3
  TFile * fs230 = new TFile("Selected/sidis_select230.root","RECREATE");
  TTree * Ts230 = Tp->CloneTree(0);
  Ts230->SetDirectory(fs230);
  TFile * fs235 = new TFile("Selected/sidis_select235.root","RECREATE");
  TTree * Ts235 = Tp->CloneTree(0);
  Ts235->SetDirectory(fs235);
  TFile * fs240 = new TFile("Selected/sidis_select240.root","RECREATE");
  TTree * Ts240 = Tp->CloneTree(0);
  Ts240->SetDirectory(fs240); 
  TFile * fs245 = new TFile("Selected/sidis_select245.root","RECREATE");
  TTree * Ts245 = Tp->CloneTree(0);
  Ts245->SetDirectory(fs245); 
  TFile * fs250 = new TFile("Selected/sidis_select250.root","RECREATE");
  TTree * Ts250 = Tp->CloneTree(0);
  Ts250->SetDirectory(fs250);
  TFile * fs255 = new TFile("Selected/sidis_select255.root","RECREATE");
  TTree * Ts255 = Tp->CloneTree(0);
  Ts255->SetDirectory(fs255);
  TFile * fs260 = new TFile("Selected/sidis_select260.root","RECREATE");
  TTree * Ts260 = Tp->CloneTree(0);
  Ts260->SetDirectory(fs260);
  TFile * fs265 = new TFile("Selected/sidis_select265.root","RECREATE");
  TTree * Ts265 = Tp->CloneTree(0);
  Ts265->SetDirectory(fs265);

  // 3 < Q2 < 4
  TFile * fs330 = new TFile("Selected/sidis_select330.root","RECREATE");
  TTree * Ts330 = Tp->CloneTree(0);
  Ts330->SetDirectory(fs330);
  TFile * fs335 = new TFile("Selected/sidis_select335.root","RECREATE");
  TTree * Ts335 = Tp->CloneTree(0);
  Ts335->SetDirectory(fs335);
  TFile * fs340 = new TFile("Selected/sidis_select340.root","RECREATE");
  TTree * Ts340 = Tp->CloneTree(0);
  Ts340->SetDirectory(fs340); 
  TFile * fs345 = new TFile("Selected/sidis_select345.root","RECREATE");
  TTree * Ts345 = Tp->CloneTree(0);
  Ts345->SetDirectory(fs345); 
  TFile * fs350 = new TFile("Selected/sidis_select350.root","RECREATE");
  TTree * Ts350 = Tp->CloneTree(0);
  Ts350->SetDirectory(fs350);
  TFile * fs355 = new TFile("Selected/sidis_select355.root","RECREATE");
  TTree * Ts355 = Tp->CloneTree(0);
  Ts355->SetDirectory(fs355);
  TFile * fs360 = new TFile("Selected/sidis_select360.root","RECREATE");
  TTree * Ts360 = Tp->CloneTree(0);
  Ts360->SetDirectory(fs360);
  TFile * fs365 = new TFile("Selected/sidis_select365.root","RECREATE");
  TTree * Ts365 = Tp->CloneTree(0);
  Ts365->SetDirectory(fs365);

  // 4 < Q2 < 5
  TFile * fs430 = new TFile("Selected/sidis_select430.root","RECREATE");
  TTree * Ts430 = Tp->CloneTree(0);
  Ts430->SetDirectory(fs430);
  TFile * fs435 = new TFile("Selected/sidis_select435.root","RECREATE");
  TTree * Ts435 = Tp->CloneTree(0);
  Ts435->SetDirectory(fs435);
  TFile * fs440 = new TFile("Selected/sidis_select440.root","RECREATE");
  TTree * Ts440 = Tp->CloneTree(0);
  Ts440->SetDirectory(fs440); 
  TFile * fs445 = new TFile("Selected/sidis_select445.root","RECREATE");
  TTree * Ts445 = Tp->CloneTree(0);
  Ts445->SetDirectory(fs445); 
  TFile * fs450 = new TFile("Selected/sidis_select450.root","RECREATE");
  TTree * Ts450 = Tp->CloneTree(0);
  Ts450->SetDirectory(fs450);
  TFile * fs455 = new TFile("Selected/sidis_select455.root","RECREATE");
  TTree * Ts455 = Tp->CloneTree(0);
  Ts455->SetDirectory(fs455);
  TFile * fs460 = new TFile("Selected/sidis_select460.root","RECREATE");
  TTree * Ts460 = Tp->CloneTree(0);
  Ts460->SetDirectory(fs460);
  TFile * fs465 = new TFile("Selected/sidis_select465.root","RECREATE");
  TTree * Ts465 = Tp->CloneTree(0);
  Ts465->SetDirectory(fs465);

  // 5 < Q2 < 6
  TFile * fs530 = new TFile("Selected/sidis_select530.root","RECREATE");
  TTree * Ts530 = Tp->CloneTree(0);
  Ts530->SetDirectory(fs530);
  TFile * fs535 = new TFile("Selected/sidis_select535.root","RECREATE");
  TTree * Ts535 = Tp->CloneTree(0);
  Ts535->SetDirectory(fs535);
  TFile * fs540 = new TFile("Selected/sidis_select540.root","RECREATE");
  TTree * Ts540 = Tp->CloneTree(0);
  Ts540->SetDirectory(fs540); 
  TFile * fs545 = new TFile("Selected/sidis_select545.root","RECREATE");
  TTree * Ts545 = Tp->CloneTree(0);
  Ts545->SetDirectory(fs545); 
  TFile * fs550 = new TFile("Selected/sidis_select550.root","RECREATE");
  TTree * Ts550 = Tp->CloneTree(0);
  Ts550->SetDirectory(fs550);
  TFile * fs555 = new TFile("Selected/sidis_select555.root","RECREATE");
  TTree * Ts555 = Tp->CloneTree(0);
  Ts555->SetDirectory(fs555);
  TFile * fs560 = new TFile("Selected/sidis_select560.root","RECREATE");
  TTree * Ts560 = Tp->CloneTree(0);
  Ts560->SetDirectory(fs560);
  TFile * fs565 = new TFile("Selected/sidis_select565.root","RECREATE");
  TTree * Ts565 = Tp->CloneTree(0);
  Ts565->SetDirectory(fs565);

  // 6 < Q2 < 8
  TFile * fs630 = new TFile("Selected/sidis_select630.root","RECREATE");
  TTree * Ts630 = Tp->CloneTree(0);
  Ts630->SetDirectory(fs630);
  TFile * fs635 = new TFile("Selected/sidis_select635.root","RECREATE");
  TTree * Ts635 = Tp->CloneTree(0);
  Ts635->SetDirectory(fs635);
  TFile * fs640 = new TFile("Selected/sidis_select640.root","RECREATE");
  TTree * Ts640 = Tp->CloneTree(0);
  Ts640->SetDirectory(fs640); 
  TFile * fs645 = new TFile("Selected/sidis_select645.root","RECREATE");
  TTree * Ts645 = Tp->CloneTree(0);
  Ts645->SetDirectory(fs645); 
  TFile * fs650 = new TFile("Selected/sidis_select650.root","RECREATE");
  TTree * Ts650 = Tp->CloneTree(0);
  Ts650->SetDirectory(fs650);
  TFile * fs655 = new TFile("Selected/sidis_select655.root","RECREATE");
  TTree * Ts655 = Tp->CloneTree(0);
  Ts655->SetDirectory(fs655);
  TFile * fs660 = new TFile("Selected/sidis_select660.root","RECREATE");
  TTree * Ts660 = Tp->CloneTree(0);
  Ts660->SetDirectory(fs660);
  TFile * fs665 = new TFile("Selected/sidis_select665.root","RECREATE");
  TTree * Ts665 = Tp->CloneTree(0);
  Ts665->SetDirectory(fs665);

  // Loop the Big Tree and select to small trees
  for (i = 0; i < N_entries; i++){
    Tp->GetEntry(i);
    // if ((x < 0.06 || x > 0.6)) continue;
    // if ((Q2 < 1.0 || Q2 > 8.0)) continue;
    // if ((z < 0.30 || z > 0.70)) continue;
    // if ((Pt > 1.2)) continue;
    
    if (Q2 > 1.0 && Q2 < 2.0) {
      if ( z > 0.3 && z < 0.35) Ts130->Fill();
      else if ( z > 0.35 && z < 0.4) Ts135->Fill();
      else if ( z > 0.4 && z < 0.45) Ts140->Fill();
      else if ( z > 0.45 && z < 0.5) Ts145->Fill();
      else if ( z > 0.5 && z < 0.55) Ts150->Fill();
      else if ( z > 0.55 && z < 0.6) Ts155->Fill();
      else if ( z > 0.6 && z < 0.65) Ts160->Fill();
      else if ( z > 0.65 && z < 0.7) Ts165->Fill();
    }
    else if (Q2 > 2.0 && Q2 < 3.0) {
      if ( z > 0.3 && z < 0.35) Ts230->Fill();
      else if ( z > 0.35 && z < 0.4) Ts235->Fill();
      else if ( z > 0.4 && z < 0.45) Ts240->Fill();
      else if ( z > 0.45 && z < 0.5) Ts245->Fill();
      else if ( z > 0.5 && z < 0.55) Ts250->Fill();
      else if ( z > 0.55 && z < 0.6) Ts255->Fill();
      else if ( z > 0.6 && z < 0.65) Ts260->Fill();
      else if ( z > 0.65 && z < 0.7) Ts265->Fill();
    }
    else if (Q2 > 3.0 && Q2 < 4.0) {
      if ( z > 0.3 && z < 0.35) Ts330->Fill();
      else if ( z > 0.35 && z < 0.4) Ts335->Fill();
      else if ( z > 0.4 && z < 0.45) Ts340->Fill();
      else if ( z > 0.45 && z < 0.5) Ts345->Fill();
      else if ( z > 0.5 && z < 0.55) Ts350->Fill();
      else if ( z > 0.55 && z < 0.6) Ts355->Fill();
      else if ( z > 0.6 && z < 0.65) Ts360->Fill();
      else if ( z > 0.65 && z < 0.7) Ts365->Fill();
    }
    else if (Q2 > 4.0 && Q2 < 5.0) {
      if ( z > 0.3 && z < 0.35) Ts430->Fill();
      else if ( z > 0.35 && z < 0.4) Ts435->Fill();
      else if ( z > 0.4 && z < 0.45) Ts440->Fill();
      else if ( z > 0.45 && z < 0.5) Ts445->Fill();
      else if ( z > 0.5 && z < 0.55) Ts450->Fill();
      else if ( z > 0.55 && z < 0.6) Ts455->Fill();
      else if ( z > 0.6 && z < 0.65) Ts460->Fill();
      else if ( z > 0.65 && z < 0.7) Ts465->Fill();
    }
    else if (Q2 > 5.0 && Q2 < 6.0) {
      if ( z > 0.3 && z < 0.35) Ts530->Fill();
      else if ( z > 0.35 && z < 0.4) Ts535->Fill();
      else if ( z > 0.4 && z < 0.45) Ts540->Fill();
      else if ( z > 0.45 && z < 0.5) Ts545->Fill();
      else if ( z > 0.5 && z < 0.55) Ts550->Fill();
      else if ( z > 0.55 && z < 0.6) Ts555->Fill();
      else if ( z > 0.6 && z < 0.65) Ts560->Fill();
      else if ( z > 0.65 && z < 0.7) Ts565->Fill();
    }
    else if (Q2 > 6.0 && Q2 < 8.0) {
      if ( z > 0.3 && z < 0.35) Ts530->Fill();
      else if ( z > 0.35 && z < 0.4) Ts635->Fill();
      else if ( z > 0.4 && z < 0.45) Ts640->Fill();
      else if ( z > 0.45 && z < 0.5) Ts645->Fill();
      else if ( z > 0.5 && z < 0.55) Ts650->Fill();
      else if ( z > 0.55 && z < 0.6) Ts655->Fill();
      else if ( z > 0.6 && z < 0.65) Ts660->Fill();
      else if ( z > 0.65 && z < 0.7) Ts665->Fill();
    }
  }
   
  fs130->Write();
  fs135->Write();
  fs140->Write();
  fs145->Write();
  fs150->Write();
  fs155->Write();
  fs160->Write();
  fs165->Write();

  fs230->Write();
  fs235->Write();
  fs240->Write();
  fs245->Write();
  fs250->Write();
  fs255->Write();
  fs260->Write();
  fs265->Write();

  fs330->Write();
  fs335->Write();
  fs340->Write();
  fs345->Write();
  fs350->Write();
  fs355->Write();
  fs360->Write();
  fs365->Write();

  fs430->Write();
  fs435->Write();
  fs440->Write();
  fs445->Write();
  fs450->Write();
  fs455->Write();
  fs460->Write();
  fs465->Write();

  fs530->Write();
  fs535->Write();
  fs540->Write();
  fs545->Write();
  fs550->Write();
  fs555->Write();
  fs560->Write();
  fs565->Write();

  fs630->Write();
  fs635->Write();
  fs640->Write();
  fs645->Write();
  fs650->Write();
  fs655->Write();
  fs660->Write();
  fs665->Write();

  return 0;
}
