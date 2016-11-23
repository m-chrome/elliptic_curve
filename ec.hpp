#ifndef EC_HPP
#define EC_HPP

#include <gcrypt.h>
#include <iostream>

struct
{
    const char *p = "115792089237316195423570985008687907853269984665640564039457584007913129639319";
    const char *a = "87789765485885808793369751294406841171614589925193456909855962166505018127157";
    const char *b = "18713751737015403763890503457318596560459867796169830279162511461744901002515";
    const char *q = "28948022309329048855892746252171976963338560298092253442512153408785530358887";
} ec_param;

struct
{
    // Координаты в проективной форме
    const char *x = "65987350182584560790308640619586834712105545126269759365406768962453298326056";
    const char *y = "22855189202984962870421402504110399293152235382908105741749987405721320435292";
} point_param;

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
    void check_p_point();
    void generate_k_number();
    void compute_k_point();
    bool check_Q_belongs();

    // Проверка корректности программы
    bool check_correction();

};


#endif // EC_HPP
