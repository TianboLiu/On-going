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
    cout << re.FindChi2Limit("Paraset/collinssample0.root", 0.955, 1) << endl;
    cout << re.FindChi2Limit("Paraset/collinssample1.root", 0.955, 1) << endl;
  }
  if (true){//get tensor charge
    Lsample re(2);
    double tc[3][6];
    re.GetFith1charge(2.41, tc, "Paraset/collinssample0.root", 0.955, 1);
    printf("%.4f  %.4f  %.4f\n", tc[0][0], tc[1][0] - tc[0][0], tc[2][0] - tc[0][0]);
    printf("%.4f  %.4f  %.4f\n", tc[0][1], tc[1][1] - tc[0][1], tc[2][1] - tc[0][1]);
    re.GetFith1charge(2.41, tc, "Paraset/collinssample1.root", 0.955, 1);
    printf("%.4f  %.4f  %.4f\n", tc[0][0], tc[1][0] - tc[0][0], tc[2][0] - tc[0][0]);
    printf("%.4f  %.4f  %.4f\n", tc[0][1], tc[1][1] - tc[0][1], tc[2][1] - tc[0][1]);
  }

  return 0;
}
