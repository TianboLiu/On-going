#include "Lsample.h"

using namespace std;

int main(){
  Lsample a2(2);
  a2.GetData("Neutron/solidneutron.root", 2);
  a2.PrintData(23);
  
  return 0;
}
