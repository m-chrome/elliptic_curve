#include <ec.hpp>
#include <iostream>
#include <gcrypt.h>

#define MODE_GOST 0
#define MODE_RANDOM 1

using namespace std;

struct
{
    const char *x = "34F0AC4BDA94D3767280F1C0613B67AA55DEBA50CF8389D31747497FE3E65EA6";
    const char *y = "6BD97C34CDEA128F611A66D6CB7257A7174624EEF8B541587D17F83B0697299F";
    const char *z = "01";
} test;

int main()
{
    EllipticCurve curve;
    curve.build_point(MODE_GOST);

    Point O;
    gcry_mpi_set_ui(O.x, 0);
    gcry_mpi_set_ui(O.y, 1);
    gcry_mpi_set_ui(O.z, 0);

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
    curve.comp_mult_point(KP3, curve.P, curve.m);
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

    cout << "Extra проверка правильности всего происходящего: ";
    if (curve.extra_check())
        cout << "[OK]\n";
    else
        cout << "[NOT OK]\n";
    return 0;
}

