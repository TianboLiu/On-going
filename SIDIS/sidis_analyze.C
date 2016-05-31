#include "Lanalysis.h"

using namespace std;

int main(int argc, char* argv[]){
  //Configuration: N11, N8, P11, P8
  TString Acontrol = argv[1];

  /*********************************/
  if (true){//make bin info tree}
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
  




  return 0;
}
