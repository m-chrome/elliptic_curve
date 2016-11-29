#include <ec.hpp>
#include <iostream>
#include <gcrypt.h>

using namespace std;

int main()
{
    EllipticCurve curve;
    curve.build_point(0);
    Point DP = curve.doubling_point(curve.P);
    if (curve.check_projective_point_belongs(DP) == 0)
        cout << "Операция удвоения работает.\n";
    else
        cout << "Операция удвоения не работает.\n";
    return 0;
}

