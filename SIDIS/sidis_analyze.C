#include "Lanalysis.h"

using namespace std;

int main(int argc, char* argv[]){
  //Configuration: N11, N8, P11, P8
  TString Acontrol = argv[1];

  /*********************************/
  if (false){//make bin info tree}
    TString bininfodir = "/var/phy/project/mepg/tl190/SoLID-SIDIS/";
    Lanalysis::MakeBinInfoTree(bininfodir+"bin_info_3he_11p.dat", "Neutron/bininfo_3He_11p.root", 11.0);
    Lanalysis::MakeBinInfoTree(bininfodir+"bin_info_3he_11m.dat", "Neutron/bininfo_3He_11m.root", 11.0);
    Lanalysis::MakeBinInfoTree(bininfodir+"bin_info_3he_8p.dat", "Neutron/bininfo_3He_8p.root", 8.8);
    Lanalysis::MakeBinInfoTree(bininfodir+"bin_info_3he_8m.dat", "Neutron/bininfo_3He_8m.root", 8.8);
    Lanalysis::MakeBinInfoTree(bininfodir+"bin_info_nh3_11p.dat", "Proton/bininfo_NH3_11p.root", 11.0);
    Lanalysis::MakeBinInfoTree(bininfodir+"bin_info_nh3_11m.dat", "Proton/bininfo_NH3_11m.root", 11.0);
    Lanalysis::MakeBinInfoTree(bininfodir+"bin_info_nh3_8p.dat", "Proton/bininfo_NH3_8p.root", 8.8);
    Lanalysis::MakeBinInfoTree(bininfodir+"bin_info_nh3_8m.dat", "Proton/bininfo_NH3_8m.root", 8.8);
  }

  /*********************************/
  if (false){// bin analysis N11
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526N11";
    TString binfilep = "Neutron/bininfo_3He_11p.root";
    TString binfilem = "Neutron/bininfo_3He_11m.root";
    Lanalysis N11(datadir, binfilep, binfilem);
    double lumi = 1.0e+10 * pow(0.197327, 2);
    double days = 48.0;
    double ST = 0.6;
    N11.SetSimInfo(lumi, days, ST, 10);
    N11.BinAnalysisNeutron("Neutron/sidisbin_11.root");
  }
  if (false){// bin analysis N8
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526N8";
    TString binfilep = "Neutron/bininfo_3He_8p.root";
    TString binfilem = "Neutron/bininfo_3He_8m.root";
    Lanalysis N8(datadir, binfilep, binfilem);
    double lumi = 1.0e+10 * pow(0.197327, 2);
    double days = 21.0;
    double ST = 0.6;
    N8.SetSimInfo(lumi, days, ST, 10);
    N8.BinAnalysisNeutron("Neutron/sidisbin_8.root");
  }
  if (false){// bin analysis P11
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526P11";
    TString binfilep = "Proton/bininfo_NH3_11p.root";
    TString binfilem = "Proton/bininfo_NH3_11m.root";
    Lanalysis P11(datadir, binfilep, binfilem);
    double lumi = 1.0e+9 * pow(0.197327, 2);
    double days = 55.0;
    double ST = 0.7;
    P11.SetSimInfo(lumi, days, ST, 10);
    P11.BinAnalysisProton("Proton/sidisbin_11.root");
  }
  if (false){// bin analysis P8
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526P8";
    TString binfilep = "Proton/bininfo_NH3_8p.root";
    TString binfilem = "Proton/bininfo_NH3_8m.root";
    Lanalysis P8(datadir, binfilep, binfilem);
    double lumi = 1.0e+9 * pow(0.197327, 2);
    double days = 27.5;
    double ST = 0.7;
    P8.SetSimInfo(lumi, days, ST, 10);
    P8.BinAnalysisProton("Proton/sidisbin_8.root");
  }

  /*****************************************/
  if (true){// bin analysis triger prescale N11
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526N11";
    TString binfilep = "Neutron/bininfo_3He_11p.root";
    TString binfilem = "Neutron/bininfo_3He_11m.root";
    Lanalysis N11(datadir, binfilep, binfilem);
    double lumi = 1.0e+10 * pow(0.197327, 2);
    double days = 48.0;
    double ST = 0.6;
    N11.SetSimInfo(lumi, days, ST, 10);
    N11.BinAnalysisNeutronTriger("Neutron/sidisbintriger_11.root");
  }
  if (true){// bin analysis triger prescale N8
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526N8";
    TString binfilep = "Neutron/bininfo_3He_8p.root";
    TString binfilem = "Neutron/bininfo_3He_8m.root";
    Lanalysis N8(datadir, binfilep, binfilem);
    double lumi = 1.0e+10 * pow(0.197327, 2);
    double days = 21.0;
    double ST = 0.6;
    N8.SetSimInfo(lumi, days, ST, 10);
    N8.BinAnalysisNeutronTriger("Neutron/sidisbintriger_8.root");
  }
  

  /*****************************************/
  if (false){// bin resolution N11
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526N11";
    TString bintree = "Neutron/sidisbin_11.root";
    Lanalysis N11(datadir);
    double lumi = 1.0e+10 * pow(0.197327, 2);
    double days = 48.0;
    double ST = 0.6;
    N11.SetSimInfo(lumi, days, ST, 1);
    N11.GetResolutionFiles();
    N11.BinResolutionNeutron(bintree, "Neutron/sidisrms_11.root");
  }
  if (false){//bin resolution N8
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526N8";
    TString bintree = "Neutron/sidisbin_8.root";
    Lanalysis N8(datadir);
    double lumi = 1.0e+10 * pow(0.197327, 2);
    double days = 21.0;
    double ST = 0.6;
    N8.SetSimInfo(lumi, days, ST, 1);
    N8.GetResolutionFiles();
    N8.BinResolutionNeutron(bintree, "Neutron/sidisrms_8.root");
  }

  /*****************************************/
  if (false){// Random coincidence N11
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526N11";
    TString bintree = "Neutron/sidisbin_11.root";
    TString rmstree = "Neutron/sidisrms_11.root";
    Lanalysis N11(datadir);
    double lumi = 1.0e+10 * pow(0.197327, 2);
    double days = 48.0;
    double ST = 0.6;
    N11.SetSimInfo(lumi, days, ST, 1);
    N11.ECoincidenceNeutron(bintree, rmstree, "Neutron/sidiscoin_11.root");
  }
  if (false){// Random coincidence N8
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526N8";
    TString bintree = "Neutron/sidisbin_8.root";
    TString rmstree = "Neutron/sidisrms_8.root";
    Lanalysis N8(datadir);
    double lumi = 1.0e+10 * pow(0.197327, 2);
    double days = 21.0;
    double ST = 0.6;
    N8.SetSimInfo(lumi, days, ST, 1);
    N8.ECoincidenceNeutron(bintree, rmstree, "Neutron/sidiscoin_8.root");
  }
  if (false){// Random coincidence P11
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526P11";
    TString bintree = "Proton/sidisbin_11.root";
    Lanalysis P11(datadir);
    double lumi = 1.0e+9 * pow(0.197327, 2);
    double days = 55.0;
    double ST = 0.7;
    P11.SetSimInfo(lumi, days, ST, 1);
    P11.ECoincidenceProton(bintree, "", "Proton/sidiscoin_11.root");
  }
  if (false){// Random coincidence P8
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526P8";
    TString bintree = "Proton/sidisbin_8.root";
    Lanalysis P8(datadir);
    double lumi = 1.0e+9 * pow(0.197327, 2);
    double days = 27.5;
    double ST = 0.7;
    P8.SetSimInfo(lumi, days, ST, 1);
    P8.ECoincidenceProton(bintree, "", "Proton/sidiscoin_8.root");
  }

  /**********************************************/
  if (false){// Bin Acceptance N11
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526N11";
    TString bintree = "Neutron/sidisbin_11.root";
    Lanalysis N11(datadir);
    double lumi = 1.0e+10 * pow(0.197327, 2);
    double days = 48.0;
    double ST = 0.6;
    N11.SetSimInfo(lumi, days, ST, 5);
    N11.BinAcceptanceNeutron(bintree, "Neutron/sidisacc_11.root");
  }
  if (false){// Bin Acceptance N8
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526N8";
    TString bintree = "Neutron/sidisbin_8.root";
    Lanalysis N8(datadir);
    double lumi = 1.0e+10 * pow(0.197327, 2);
    double days = 21.0;
    double ST = 0.6;
    N8.SetSimInfo(lumi, days, ST, 5);
    N8.BinAcceptanceNeutron(bintree, "Neutron/sidisacc_8.root");
  }
  if (false){// Bin Acceptance P11
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526P11";
    TString bintree = "Proton/sidisbin_11.root";
    Lanalysis P11(datadir);
    double lumi = 1.0e+9 * pow(0.197327, 2);
    double days = 55.0;
    double ST = 0.7;
    P11.SetSimInfo(lumi, days, ST, 5);
    P11.BinAcceptanceProton(bintree, "Proton/sidisacc_11.root");
  }
  if (false){// Bin Acceptance P8
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526P8";
    TString bintree = "Proton/sidisbin_8.root";
    Lanalysis P8(datadir);
    double lumi = 1.0e+9 * pow(0.197327, 2);
    double days = 27.5;
    double ST = 0.7;
    P8.SetSimInfo(lumi, days, ST, 5);
    P8.BinAcceptanceNeutron(bintree, "Proton/sidisacc_8.root");
  }

  /************************************************/
  if (false){// E resolution N11
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526N11";
    TString bintree = "Neutron/sidisbin_11.root";
    TString rmstree = "Neutron/sidisrms_11.root";
    TString acctree = "Neutron/sidisacc_11.root";
    Lanalysis N11(datadir);
    double lumi = 1.0e+10 * pow(0.197327, 2);
    double days = 48.0;
    double ST = 0.6;
    N11.SetSimInfo(lumi, days, ST, 5);
    N11.EResolutionNeutron(bintree, rmstree, acctree, "Neutron/sidisres_11.root");
  }
  if (false){// E resolution N8
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526N8";
    TString bintree = "Neutron/sidisbin_8.root";
    TString rmstree = "Neutron/sidisrms_8.root";
    TString acctree = "Neutron/sidisacc_8.root";
    Lanalysis N8(datadir);
    double lumi = 1.0e+10 * pow(0.197327, 2);
    double days = 21.0;
    double ST = 0.6;
    N8.SetSimInfo(lumi, days, ST, 5);
    N8.EResolutionNeutron(bintree, rmstree, acctree, "Neutron/sidisres_8.root");
  }

  /************************************************/
  if (false){//Nuclear PDF N11
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526N11";
    TString bintree = "Neutron/sidisbin_11.root";
    Lanalysis N11(datadir);
    double lumi = 1.0e+10 * pow(0.197327, 2);
    double days = 48.0;
    double ST = 0.6;
    N11.SetSimInfo(lumi, days, ST, 5);
    N11.ENuclearPDF(bintree, "Neutron/sidisnpdf_11.root");
  }
  if (false){//Nuclear PDF N8
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526N8";
    TString bintree = "Neutron/sidisbin_8.root";
    Lanalysis N8(datadir);
    double lumi = 1.0e+10 * pow(0.197327, 2);
    double days = 21.0;
    double ST = 0.6;
    N8.SetSimInfo(lumi, days, ST, 5);
    N8.ENuclearPDF(bintree, "Neutron/sidisnpdf_8.root");
  }

  /************************************************/
  if (false){//Nuclear effect N11
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526N11";
    TString bintree = "Neutron/sidisbin_11.root";
    Lanalysis N11(datadir);
    double lumi = 1.0e+10 * pow(0.197327, 2);
    double days = 3.5/48.0;
    double ST = 0.6;
    N11.SetSimInfo(lumi, days, ST, 5);
    N11.ENuclearNeutron(bintree, "Neutron/sidisnucl_11.root");
  }
  if (false){//Nuclear effect N8
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526N8";
    TString bintree = "Neutron/sidisbin_8.root";
    Lanalysis N8(datadir);
    double lumi = 1.0e+10 * pow(0.197327, 2);
    double days = 1.5/21.0;
    double ST = 0.6;
    N8.SetSimInfo(lumi, days, ST, 5);
    N8.ENuclearNeutron(bintree, "Neutron/sidisnucl_8.root");
  }
  /************************************************/
  if (false){//Dilution factor P11
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526P11";
    TString bintree = "Proton/sidisbin_11.root";
    Lanalysis P11(datadir);
    double lumi = 1.0e+9 * pow(0.197327, 2);
    double days = 55.0;
    double ST = 0.7;
    P11.SetSimInfo(lumi, days, ST, 5);
    P11.EDilutionProton(bintree, "Proton/sidisdilu_11.root");
  }
  if (false){//Dilution factor P8
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526P8";
    TString bintree = "Proton/sidisbin_8.root";
    Lanalysis P8(datadir);
    double lumi = 1.0e+9 * pow(0.197327, 2);
    double days = 27.5;
    double ST = 0.7;
    P8.SetSimInfo(lumi, days, ST, 5);
    P8.EDilutionProton(bintree, "Proton/sidisdilu_8.root");
  }
  /************************************************/
  if (false){//Total Neutron
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526N11";
    TString bintree = "Neutron/sidisbin_11.root";
    Lanalysis Neu(datadir);
    double lumi = 1.0e+10 * pow(0.197327, 2);
    double days = 3.5/48.0;
    double ST = 0.6;
    Neu.SetSimInfo(lumi, days, ST, 5);
    Neu.ETotalNeutron("Neutron/", "Neutron/solidneutron.root");
  }
  /************************************************/
  if (false){//Total Proton
    TString datadir = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0526P11";
    TString bintree = "Proton/sidisbin_11.root";
    Lanalysis Pro(datadir);
    double lumi = 1.0e+9 * pow(0.197327, 2);
    double days = 55.0;
    double ST = 0.7;
    Pro.SetSimInfo(lumi, days, ST, 5);
    Pro.ETotalProton("Proton/", "Proton/solidproton.root");
  }
  
  return 0;
}
