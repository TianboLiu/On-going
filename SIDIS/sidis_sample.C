#include "Lsample.h"

using namespace std;

int main(){
  Lsample a2(2);
  a2.GetData("Neutron/solidneutron.root", 2);
  a2.GetData("Proton/solidproton.root", 2);
  a2.PrintData(23);
  a2.InitialA2para();
  cout << a2.Chi2A2(0) << endl;
  
  return 0;
}
