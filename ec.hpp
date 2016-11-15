#ifndef EC_HPP
#define EC_HPP

#include <gcrypt.h>
#include <iostream>

class EC
{
private:

    // Параметры эллиптической кривой
    gcry_mpi_t p;
    gcry_mpi_t a;
    gcry_mpi_t b;
    gcry_mpi_t q;
    gcry_mpi_t k;     // Случайное число 1<k<q

    // Точка на кривой
    gcry_mpi_point_t P0;    // Точка P(x0,y0)
    gcry_mpi_point_t Q;     // Кратная точка Q(x,y)=k*P=P+..+P

    // S-exp для кривой
    gcry_sexp_t weierstrass;

public:
    EC();
    ~EC();

    // Функции
    void build_point();
    void generate_k_number();
    void compute_miltiple_point();
    bool check_Q_belongs();

    // Проверка корректности программы
    bool check_correction();

};

#endif // EC_HPP
