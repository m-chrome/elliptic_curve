#include <ec.hpp>
#include <iostream>
#include <gcrypt.h>

using namespace std;

int main()
{
    EllipticCurve curve;
    curve.build_point(0);
    Point DupP = curve.doubling_point(curve.P);
    if (curve.check_point_belongs(DupP) == 0)
    {
        cout << "Точка DupP(x,y) принадлежит кривой.\n";
        cout << "x = ";
        show_mpi(DupP.x);
        cout << "y = ";
        show_mpi(DupP.y);
        cout << endl;
        return 0;
    }
    else
    {
        cout << "Точка DupP(x,y) не принадлежит кривой.\n";
        cout << "x = ";
        show_mpi(DupP.x);
        cout << "y = ";
        show_mpi(DupP.y);
        cout << endl;
        return 1;
    }
    return 0;
}

