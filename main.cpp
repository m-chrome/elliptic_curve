#include <ec.hpp>
#include <iostream>
#include <gcrypt.h>

#define MODE_GOST 0
#define MODE_RANDOM 1

using namespace std;

int main()
{
    EllipticCurve curve;
    curve.build_point(MODE_GOST);

    Point DP = curve.doubling_point(curve.P);
    if (curve.check_projective_point_belongs(DP) == 0)
        cout << "Удвоенная точка P(x,y) принадлежит кривой.\n\n";
    else
        cout << "Удвоенная точка P(x,y) не принадлежит кривой.\n\n";

    Point DQ = curve.doubling_point(curve.Q);
    if (curve.check_projective_point_belongs(DQ) == 0)
        cout << "Удвоенная точка Q(x,y) принадлежит кривой.\n\n";
    else
        cout << "Удвоенная точка Q(x,y) не принадлежит кривой.\n\n";

    Point SP = curve.add_points(DQ, DP);
    if (curve.check_projective_point_belongs(SP) == 0)
        cout << "Точка SP = P + Q принадлежит кривой.\n\n";
    else
        cout << "Точка SP = P + Q не принадлежит кривой.\n\n";
    return 0;
}

