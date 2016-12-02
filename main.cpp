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

    Point Q1, Q2;

    gcry_mpi_set_ui(Q1.x, 0);
    gcry_mpi_set_ui(Q1.y, 1);
    gcry_mpi_set_ui(Q1.z, 1);

    curve.doubling_point(Q1, Q1);

    Q1.print();

    cout << "Тестирование операции удвоения:\n";
    Point DP;
    curve.doubling_point(DP, curve.P);
    if (curve.check_projective_point_belongs(DP) == 0)
    {
        cout << "Удвоенная точка принадлежит кривой.\n";
    }
    else
    {
        cout << "Удвоенная точка не принадлежит кривой.\n";
        return 1;
    }
    DP.print();
    cout << endl;

    cout << "Тестирование операции сложения точек:\n";
    Point SumP;
    curve.add_points(SumP, curve.P, DP);
    if (curve.check_projective_point_belongs(DP) == 0)
    {
        cout << "Полученная точка принадлежит кривой.\n";
    }
    else
    {
        cout << "Полученная точка не принадлежит кривой.\n";
        return 1;
    }
    SumP.print();
    cout << endl;

    cout << "Тестирование операции нахождения кратной точки 1:\n";
    Point KP1;
    curve.comp_mult_point(KP1, curve.P, curve.k);
    if (curve.check_projective_point_belongs(DP) == 0)
    {
        cout << "Кратная точка принадлежит кривой.\n";
    }
    else
    {
        cout << "Кратная точка не принадлежит кривой.\n";
    }
    KP1.print();
    cout << endl;

    cout << "Тестирование операции нахождения кратной точки 2:\n";
    Point KP2;
    curve.comp_mult_point(KP2, curve.P, curve.l);
    if (curve.check_projective_point_belongs(DP) == 0)
    {
        cout << "Кратная точка принадлежит кривой.\n";
    }
    else
    {
        cout << "Кратная точка не принадлежит кривой.\n";
    }
    KP2.print();
    cout << endl;

    cout << "Extra проверка правильности всего происходящего: ";
    if (curve.extra_check())
        cout << "[OK]\n";
    else
        cout << "[NOT OK]\n";
    return 0;
}

