#include "Lsample.h"

using namespace std;

int main(){
  if (false){
    Lsample a2(2);
    a2.GetData("Neutron/solidneutron.root", 2);
    a2.GetData("Proton/solidproton.root", 2);
    a2.InitialA2para();
    a2.UniformSampleA2(1, 100000.0, "Paraset/collinssample1.root");
    a2.UniformSampleA2(0, 100000.0, "Paraset/collinssample0.root");
  }
  if (true){
    Lsample re(2);
    cout << re.FindChi2Limit("Paraset/collinssample0.root", 0.955, 0) << endl;
    cout << re.FindChi2Limit("Paraset/collinssample1.root", 0.955, 1) << endl;
  }

  return 0;
}
