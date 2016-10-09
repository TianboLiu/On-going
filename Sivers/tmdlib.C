#include "sivers.h"

using namespace std;

//test code
int main(){
  double kin[4] = {0.15, 2.4, 0.5, 0.3};
  cout << AUT(kin, "proton", "pip") << endl;
  return 0;
}
