#include <Matrix.h>
#include <iostream>

using namespace std;
using namespace SPTH;

int main() {
    A4DMatrix A = Translation(0.2, 0.2, 1).Transpose();
    A4DMatrix B = Translation(0.5, 0.5, 0.25).Transpose();
    A4DMatrix C = ZRotation4D(-90).Transpose();
    A4DMatrix D = XRotation4D(-180).Transpose();
    cout << A * B * C * D << endl;
    cout << "==========" << endl;
    cout << D * C * B * A << endl;

    A = Translation(-0.3, 1.5, 0.8).Transpose();
    B = XRotation4D(180).Transpose();
    C = ZRotation4D(70).Transpose();
    D = YRotation4D(40).Transpose();
    cout << A * B * C * D << endl;
    cout << "==========" << endl;
    cout << D * C * B * A << endl;
    return 0;
}

