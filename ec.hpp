#ifndef EC_HPP
#define EC_HPP

#include <gcrypt.h>
#include <iostream>

// Набор параметров точки и кривой (HEX)
struct
{
    const char *x = "91E38443A5E82C0D880923425712B2BB658B9196932E02C78B2582FE742DAA28";
    const char *y = "32879423AB1A0375895786C4BB46E9565FDE0B5344766740AF268ADB32322E5C";
} point_param;

struct
{
    const char *p = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD97";
    const char *a = "C2173F1513981673AF4892C23035A27CE25E2013BF95AA33B22C656F277E7335";
    const char *b = "295F9BAE7428ED9CCC20E7C359A9D41A22FCCD9108E17BF7BA9337A6F8AE9513";
    const char *q = "400000000000000000000000000000000FD8CDDFC87B6635C115AF556C360C67";
} ec_param;


class Point
{
public:
    gcry_mpi_t x;
    gcry_mpi_t y;
    gcry_mpi_t z;

    Point();
    ~Point();
};

// Вейерштрасс: y^2 = x^3 + ax + b (p)
class EllipticCurve
{
public:
    gcry_mpi_t p;
    gcry_mpi_t a;
    gcry_mpi_t b;
    gcry_mpi_t q;
    gcry_mpi_t k;

    Point P;
    Point Q;

    EllipticCurve();
    ~EllipticCurve();

    // Посчитать y^2 = f(x) (p)
    gcry_mpi_t comp_fx0(gcry_mpi_t x0);

    // Критерий Эйлера на то, f(x) - квадратичный вычет в p или нет?
    int euler_criteria(gcry_mpi_t fx0);

    // Вычисление y = f(x)^((p+1)/4) (p) при условии, что p = 3 (4)
    gcry_mpi_t comp_y0(gcry_mpi_t x0);

    // Проверка принадлежности точки к кривой
    int check_point_belongs(const Point &point);

    // Построение точки на кривой
    int build_point(int mode);

    // Операция удвоения точки в аффинной форме
    Point doubling_point(const Point &point);

    // Операция сложения двух точек в аффинной форме
    Point add_points(const Point &point1, const Point &point2);

    // Вычисление кратной точки Q = k*P
    void comp_mult_point(const Point &point);
};

#endif // EC_HPP
