#ifndef EC_HPP
#define EC_HPP

#include <gcrypt.h>
#include <iostream>

struct
{
    const char *p = "57896044618658097711785492504343953926634992332820282019728792003956564821041";
    const char *a = "7";
    const char *b = "43308876546767276905765904595650931995942111794451039583252968842033849580414";
    const char *q = "57896044618658097711785492504343953927082934583725450622380973592137631069619";
} ec_param;

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
