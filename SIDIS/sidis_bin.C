#include <iostream>
#include <fstream>

#include "Lanalysis.h"

using namespace std;

int main(){
  if (false){
    TString bininfofile11p = "/var/phy/project/mepg/tl190/SIDIS/11GeV/newBins/bin_info_nh3_11p.dat";
    TString bintreefile11p = "bininfo_NH3_11p.root";
    TString bininfofile11m = "/var/phy/project/mepg/tl190/SIDIS/11GeV/newBins/bin_info_nh3_11m.dat";
    TString bintreefile11m = "bininfo_NH3_11m.root";
    TString datadir11 = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0229E11Proton";
    TString sidis11 = "sidisbin_11.root";
    lan->makebininfotree(bininfofile11p, bintreefile11p, 11.0);
    lan->makebininfotree(bininfofile11m, bintreefile11m, 11.0);
    lan->binanalysis(bintreefile11p, bintreefile11m, sidis11, datadir11);
  }

  if (true){
    TString bininfofile8p = "/var/phy/project/mepg/tl190/SIDIS/11GeV/newBins/bin_info_nh3_8p.dat";
    TString bintreefile8p = "bininfo_NH3_8p.root";
    TString bininfofile8m = "/var/phy/project/mepg/tl190/SIDIS/11GeV/newBins/bin_info_nh3_8m.dat";
    TString bintreefile8m = "bininfo_NH3_8m.root";
    TString datadir8 = "/var/phy/project/mepg/tl190/SoLID-cluster/RUN0229E8Proton";
    TString sidis8 = "sidisbin_8.root";
    lan->makebininfotree(bininfofile8p, bintreefile8p, 8.8);
    lan->makebininfotree(bininfofile8m, bintreefile8m, 8.8);
    lan->binanalysis(bintreefile8p, bintreefile8m, sidis8, datadir8);
  }

  return 0;
}
