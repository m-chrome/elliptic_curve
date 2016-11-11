#ifndef EC_HPP
#define EC_HPP

#include <gcrypt.h>

class EC
{
private:

    // Параметры эллиптической кривой
    gcry_mpi_t p;
    gcry_mpi_t a;
    gcry_mpi_t b;

    // Прочие числа
    gcry_mpi_point_t k;     // Случайное число 1<k<q
    gcry_mpi_point_t q;     // Порядок циклической подгруппы точек ЭК

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
