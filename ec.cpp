#include <ec.hpp>
#include <iostream>

using namespace std;

EC::EC()
{
    // Параметры из ГОСТ-34.10-2012, пример 1
    cout << "Эллиптическая кривая" << endl;
    cout << "Форма Вейерштрасса: y^2 = x^3 + ax + b (mod p)" << endl;
    cout << "Параметры кривой:\n" << endl;
    cout << "\tp = 57896044618658097711785492504343953926634992332820282019728792003956564821041" << endl;
    cout << "\ta = 7" << endl;
    cout << "\tb = 43308876546767276905765904595650931995942111794451039583252968842033849580414" << endl;
    cout << "\tq = 57896044618658097711785492504343953927082934583725450622380973592137631069619" << endl;

    p = gcry_mpi_new(0);
    a = gcry_mpi_new(0);
    b = gcry_mpi_new(0);
    q = gcry_mpi_new(0);
    gcry_mpi_set(p, 57896044618658097711785492504343953926634992332820282019728792003956564821041);
    gcry_mpi_set(a, 7);
    gcry_mpi_set(b, 43308876546767276905765904595650931995942111794451039583252968842033849580414);
    gcry_mpi_set(q, 57896044618658097711785492504343953927082934583725450622380973592137631069619);
    cout << "[System] Присвоение параметров." << endl;
}

EC::~EC()
{
    cout << "[System] Очистка памяти." << endl;
    gcry_mpi_release(p);
    gcry_mpi_release(a);
    gcry_mpi_release(b);
    gcry_mpi_release(q);
    gcry_mpi_point_release(P0);
    gcry_mpi_point_release(Q);
}

void EC::build_point();
void EC::generate_k_number();
void EC::compute_miltiple_point();
void EC::check_Q_belongs();
bool EC::check_correction();
