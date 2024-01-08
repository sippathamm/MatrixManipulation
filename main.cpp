#include <Matrix.h>
#include <iostream>

using namespace std;
using namespace SPTH;

int main() {
    A4DMatrix Z = ZRotation4D(60);
    A4DMatrix Y = YRotation4D(45);
    A4DMatrix X = XRotation4D(30);
    cout << Z.Transpose() * Y.Transpose() * X.Transpose() << endl;
    return 0;
}

