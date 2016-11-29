#include <ec.hpp>
#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <gcrypt.h>
#include <gmp.h>

#define MAX_MPI_BUF 200

using namespace std;

// -----------------------------------------------------
// Общие методы
void show_mpi(gcry_mpi_t mpi)
{
    unsigned char *check_buf = (unsigned char *)malloc(MAX_MPI_BUF);
    memset(check_buf, '\0', MAX_MPI_BUF);
    gcry_mpi_print(GCRYMPI_FMT_HEX, check_buf, MAX_MPI_BUF, NULL, mpi);
    printf("%s\n", check_buf);
}

// -----------------------------------------------------
// Методы точки

Point::Point()
{
    x = gcry_mpi_new(0);
    y = gcry_mpi_new(0);
    z = gcry_mpi_new(0);
}

Point::~Point()
{
    gcry_mpi_release(x);
    gcry_mpi_release(y);
    gcry_mpi_release(z);
}

// -----------------------------------------------------
// Методы эллиптической кривой

EllipticCurve::EllipticCurve()
{
    cout << "Эллиптическая кривая" << endl;
    cout << "Форма Вейерштрасса: y^2 = x^3 + ax + b (mod p)" << endl;
    cout << "Параметры кривой:\n" << endl;
    p = gcry_mpi_new(0);
    a = gcry_mpi_new(0);
    b = gcry_mpi_new(0);
    q = gcry_mpi_new(0);
    if (gcry_mpi_scan(&p, GCRYMPI_FMT_HEX, ec_param.p, 0, NULL) != 0)
        cout << "Ошибка записи p.\n";
    if (gcry_mpi_scan(&a, GCRYMPI_FMT_HEX, ec_param.a, 0, NULL) != 0)
        cout << "Ошибка записи a.\n";
    if (gcry_mpi_scan(&b, GCRYMPI_FMT_HEX, ec_param.b, 0, NULL) != 0)
        cout << "Ошибка записи b.\n";
    if (gcry_mpi_scan(&q, GCRYMPI_FMT_HEX, ec_param.q, 0, NULL) != 0)
        cout << "Ошибка записи q.\n";
    cout << "p = ";
    show_mpi(p);
    cout << "a = ";
    show_mpi(a);
    cout << "b = ";
    show_mpi(b);
    cout << "q = ";
    show_mpi(q);

    // Генерация случайного числа k
    k = gcry_mpi_new(0);
    /*gcry_mpi_randomize(k, gcry_mpi_get_nbits(q), GCRY_WEAK_RANDOM);
    gcry_mpi_mod(k, k, q);
    cout << "k = ";
    show_mpi(k);*/
    cout << endl;
}

EllipticCurve::~EllipticCurve()
{
    cout << "Очистка памяти.\n";
    gcry_mpi_release(p);
    gcry_mpi_release(a);
    gcry_mpi_release(b);
    gcry_mpi_release(q);
    gcry_mpi_release(k);
}

gcry_mpi_t EllipticCurve::comp_fx0(gcry_mpi_t x0)
{
    gcry_mpi_t fx0 = gcry_mpi_new(0);
    gcry_mpi_t temp = gcry_mpi_new(0);
    gcry_mpi_set_ui(temp, 3);
    gcry_mpi_powm(fx0, x0, temp, p);
    gcry_mpi_mulm(temp, a, x0, p);
    gcry_mpi_addm(temp, temp, b, p);
    gcry_mpi_addm(fx0, fx0, temp, p);
    gcry_mpi_release(temp);
    return fx0;
}

int EllipticCurve::euler_criteria(gcry_mpi_t fx0)
{
    gcry_mpi_t exp = gcry_mpi_new(0);
    gcry_mpi_t div = gcry_mpi_new(0);
    gcry_mpi_t X = gcry_mpi_new(0);
    gcry_mpi_sub_ui(exp, p, 1);
    gcry_mpi_set_ui(div, 2);
    gcry_mpi_div(exp, NULL, exp, div, 0);
    gcry_mpi_powm(X, fx0, exp, p);
    gcry_mpi_release(exp);
    gcry_mpi_release(div);
    int cmp = gcry_mpi_cmp_ui(X,1);
    gcry_mpi_release(X);
    return cmp;
}

gcry_mpi_t EllipticCurve::comp_y0(gcry_mpi_t x0)
{
    gcry_mpi_t fx0 = comp_fx0(x0);
    gcry_mpi_t y0 = gcry_mpi_new(0);
    if (euler_criteria(fx0)==0)
    {
        cout << "f(x) - квадратичный вычет в р.\n";
        gcry_mpi_t exp = gcry_mpi_new(0);
        gcry_mpi_t div = gcry_mpi_new(0);
        gcry_mpi_add_ui(exp, p, 1);
        gcry_mpi_set_ui(div, 4);
        gcry_mpi_div(exp, NULL, exp, div, 0);
        gcry_mpi_powm(y0, fx0, exp, p);
        gcry_mpi_release(exp);
        gcry_mpi_release(div);
    }
    else
    {
        cout << "f(x) - не квадратичный вычет в р.\n";
        gcry_mpi_set_ui(y0, 0);
    }
    return y0;
}

int EllipticCurve::check_affine_point_belongs(const Point &point)
{
    // Подстановка в уравнение кривой координат точки point и проверка на принадлежность
    gcry_mpi_t left = gcry_mpi_new(0);
    gcry_mpi_t exp = gcry_mpi_new(0);
    gcry_mpi_set_ui(exp, 2);
    gcry_mpi_powm(left, point.y, exp, p);
    gcry_mpi_t right = comp_fx0(point.x);
    int cmp = gcry_mpi_cmp(left, right);
    gcry_mpi_release(left);
    gcry_mpi_release(right);
    return cmp;
}

int EllipticCurve::build_point(int mode)
{
    switch(mode)
    {
    case 0:
    {
        cout << "Будут использованы ГОСТ параметры для точки.\n";
        if (gcry_mpi_scan(&P.x, GCRYMPI_FMT_HEX, point_param.x, 0, NULL) != 0)
        {
            cout << "Ошибка записи координаты x.\n";
            return 1;
        }
        cout << "x = ";
        show_mpi(P.x);
        if (gcry_mpi_scan(&P.y, GCRYMPI_FMT_HEX, point_param.y, 0, NULL) != 0)
        {
            cout << "Ошибка записи координаты x.\n";
            return 1;
        }
        cout << "y = ";
        show_mpi(P.y);
        gcry_mpi_set_ui(P.z, 1);
        cout << "z = ";
        show_mpi(P.z);
        cout << endl;
        break;
    }
    case 1:
    {
        cout << "Будет сгенерирована случайная точка.\n";
        gcry_mpi_randomize(P.x, gcry_mpi_get_nbits(p), GCRY_WEAK_RANDOM);
        P.y = comp_y0(P.x);
        cout << "x = ";
        show_mpi(P.x);
        cout << "y = ";
        show_mpi(P.y);
        gcry_mpi_set_ui(P.z, 1);
        cout << "z = ";
        show_mpi(P.z);
        cout << endl;
        break;
    }
    default:
    {
        cout << "Неверный режим.\n";
        return 1;
        break;
    }
    }
    if (check_affine_point_belongs(P) == 0)
    {
        cout << "Точка P(x,y) принадлежит кривой.\n";
        cout << "x = ";
        show_mpi(P.x);
        cout << "y = ";
        show_mpi(P.y);
        cout << "z = ";
        show_mpi(P.z);
        cout << endl;
        return 0;
    }
    else
    {
        cout << "Точка P(x,y) не принадлежит кривой.\n";
        cout << "x = ";
        show_mpi(P.x);
        cout << "y = ";
        show_mpi(P.y);
        cout << "z = ";
        show_mpi(P.z);
        cout << endl;
        return 1;
    }
}

int EllipticCurve::check_projective_point_belongs(const Point &point)
{
    // point в проективных координатах
    // Подстановка в уравнение
    // Y^2*Z = X^3 + a*X*Z^3 + b*Z^3 (p)

    // Левая часть
    gcry_mpi_t left = gcry_mpi_new(0);
    gcry_mpi_t exp = gcry_mpi_new(0);
    gcry_mpi_set_ui(exp, 2);
    gcry_mpi_powm(left, point.y, exp, p);
    gcry_mpi_mulm(left, left, point.z, p);

    // Правая часть
    gcry_mpi_t right = gcry_mpi_new(0);
    gcry_mpi_t first = gcry_mpi_new(0);
    gcry_mpi_t second = gcry_mpi_new(0);
    gcry_mpi_t third = gcry_mpi_new(0);
    gcry_mpi_powm(second, point.z, exp, p);
    gcry_mpi_mulm(second, second, point.x, p);
    gcry_mpi_mulm(second, second, a, p);

    gcry_mpi_set_ui(exp, 3);
    gcry_mpi_powm(first, point.x, exp, p);

    gcry_mpi_powm(third, point.z, exp, p);
    gcry_mpi_mulm(third, third, b, p);
    gcry_mpi_addm(right, first, second, p);
    gcry_mpi_addm(right, right, third, p);

    int cmp = gcry_mpi_cmp(left, right);
    gcry_mpi_release(left);
    gcry_mpi_release(right);
    gcry_mpi_release(exp);
    gcry_mpi_release(first);
    gcry_mpi_release(second);
    gcry_mpi_release(third);
    return cmp;
}


Point EllipticCurve::doubling_point(const Point &point)
{
    Point DoubleP;
    gcry_mpi_t t[11];
    for(int i = 0; i < 11; ++i)
        t[i] = gcry_mpi_new(0);
    gcry_mpi_t xx = gcry_mpi_new(0);
    gcry_mpi_t B = gcry_mpi_new(0);
    gcry_mpi_t zz = gcry_mpi_new(0);
    gcry_mpi_t w = gcry_mpi_new(0);
    gcry_mpi_t r = gcry_mpi_new(0);
    gcry_mpi_t rr = gcry_mpi_new(0);
    gcry_mpi_t h = gcry_mpi_new(0);
    gcry_mpi_t s = gcry_mpi_new(0);
    gcry_mpi_t sss = gcry_mpi_new(0);
    gcry_mpi_t exp = gcry_mpi_new(0);
    gcry_mpi_set_ui(exp, 2);
    gcry_mpi_powm(xx, point.x, exp, p);
    gcry_mpi_powm(zz, point.z, exp, p);
    gcry_mpi_mul_ui(t[0], xx, 3);
    gcry_mpi_mulm(t[1], a, zz, p);
    gcry_mpi_addm(w, t[0], t[1], p);
    gcry_mpi_mulm(t[2], point.z, point.y, p);
    gcry_mpi_mul_ui(s, t[2], 2);
    gcry_mpi_set_ui(exp, 3);
    gcry_mpi_powm(sss, s, exp, p);
    gcry_mpi_mulm(r, point.y, s, p);
    gcry_mpi_mulm(rr, r, r, p);
    gcry_mpi_addm(t[3], r, point.x, p);
    gcry_mpi_mulm(t[4], t[3], t[3], p);
    gcry_mpi_subm(t[5], t[4], xx, p);
    gcry_mpi_subm(B, t[5], rr, p);
    gcry_mpi_mulm(t[6], w, w, p);
    gcry_mpi_mul_ui(t[7], B, 2);
    gcry_mpi_subm(h, t[6], t[7], p);
    gcry_mpi_mulm(DoubleP.x, h, s, p);
    gcry_mpi_subm(t[8], B, h, p);
    gcry_mpi_mul_ui(t[9], rr, 2);
    gcry_mpi_mulm(t[10], t[8], w, p);
    gcry_mpi_subm(DoubleP.y, t[10], t[9], p);
    gcry_mpi_set(DoubleP.z, sss);
    return DoubleP;
}

