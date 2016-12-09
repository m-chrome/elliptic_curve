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

    cout << "Тестирование операции удвоения:\n";
    Point DP;
    curve.doubling_point(DP, curve.P);
    DP.print();
    if (curve.check_projective_point_belongs(DP) == 0)
    {
        cout << "Удвоенная точка принадлежит кривой.\n";
    }
    else
    {
        cout << "Удвоенная точка не принадлежит кривой.\n";
        return 1;
    }
    cout << endl;

    cout << "Тестирование операции сложения точек:\n";
    Point SumP;
    curve.add_points(SumP, curve.P, DP);
    SumP.print();
    if (curve.check_projective_point_belongs(SumP) == 0)
    {
        cout << "Полученная точка принадлежит кривой.\n";
    }
    else
    {
        cout << "Полученная точка не принадлежит кривой.\n";
        return 1;
    }
    cout << endl;

    cout << "Тестирование операции нахождения кратной точки 1:\n";
    Point KP1;
    curve.comp_mult_point(KP1, curve.P, curve.k);
    KP1.print();
    if (curve.check_projective_point_belongs(KP1) == 0)
    {
        cout << "Кратная точка принадлежит кривой.\n";
    }
    else
    {
        cout << "Кратная точка не принадлежит кривой.\n";
        return 1;
    }
    cout << endl;

    cout << "Тестирование операции нахождения кратной точки 2:\n";
    Point KP2;
    curve.comp_mult_point(KP2, curve.P, curve.l);
    KP2.print();
    if (curve.check_projective_point_belongs(KP2) == 0)
    {
        cout << "Кратная точка принадлежит кривой.\n";
    }
    else
    {
        cout << "Кратная точка не принадлежит кривой.\n";
        return 1;
    }
    cout << endl;

    cout << "Тестирование операции нахождения кратной точки 3:\n";
    Point KP3;
    curve.comp_mult_point(KP3, curve.P, curve.q);
    KP3.print();
    if (curve.check_projective_point_belongs(KP3) == 0)
    {
        cout << "Кратная точка принадлежит кривой.\n";
    }
    else
    {
        cout << "Кратная точка не принадлежит кривой.\n";
        return 1;
    }
    cout << endl;
    return 0;
}

