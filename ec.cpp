#include <ec.hpp>
#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <gcrypt.h>
#include <gmp.h>

#define MAX_MPI_BUF 200

using namespace std;

struct
{
    const char *x = "38461FA70F752EA90E080EEDDFA7ECF64C97CA72AB56851012BBC3956CFF3C4A";
    const char *y = "B214F2F87D1CEFD8DCF842514D13A499A9D640D4AA5FB150FB3BF4B84010FA92";
} test;

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

void Point::print()
{
    cout << "x = ";
    show_mpi(x);
    cout << "y = ";
    show_mpi(y);
    cout << "z = ";
    show_mpi(z);
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
        cout << "Будут использованы ГОСТ параметры для точки P(x,y,z).\n";
        if (gcry_mpi_scan(&P.x, GCRYMPI_FMT_HEX, point_param.x, 0, NULL) != 0)
        {
            cout << "Ошибка записи координаты x.\n";
            return 1;
        }
        if (gcry_mpi_scan(&P.y, GCRYMPI_FMT_HEX, point_param.y, 0, NULL) != 0)
        {
            cout << "Ошибка записи координаты x.\n";
            return 1;
        }
        gcry_mpi_set_ui(P.z, 1);
        P.print();
        cout << endl;
        //break;
    }
    case 1:
    {
        cout << "Будет сгенерирована случайная точка.\n";
        //gcry_mpi_randomize(P.x, gcry_mpi_get_nbits(p), GCRY_WEAK_RANDOM);
        gcry_mpi_scan(&Q.x, GCRYMPI_FMT_HEX, test.x, 0, NULL);
        Q.y = comp_y0(Q.x);
        gcry_mpi_set_ui(Q.z, 1);
        Q.print();
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
    if (check_affine_point_belongs(P) == 0 && check_affine_point_belongs(Q)==0)
    {
        cout << "Точка P(x,y) принадлежит кривой.\n";
        cout << "Точка Q(x,y) принадлежит кривой.\n";
        cout << endl;
        return 0;
    }
    else
    {
        cout << "Точка P(x,y) не принадлежит кривой.\n";
        cout << "Точка Q(x,y) не принадлежит кривой.\n";
        cout << endl;
        return 1;
    }
}

int EllipticCurve::check_projective_point_belongs(const Point &point)
{
    // point в проективных координатах
    // Подстановка в уравнение
    // Y^2*Z = X^3 + a*X*Z^2 + b*Z^3 (p)

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

    cout << "Double P:\n";
    DoubleP.print();
    return DoubleP;
}

Point EllipticCurve::add_points(const Point &p1, const Point &p2)
{
    Point p3;
    /*gcry_mpi_t t[8];
    for(int i = 0; i < 8; ++i)
        t[i] = gcry_mpi_new(0);
    gcry_mpi_t Y1Z2 = gcry_mpi_new(0);
    gcry_mpi_t X1Z2 = gcry_mpi_new(0);
    gcry_mpi_t Z1Z2 = gcry_mpi_new(0);
    gcry_mpi_mulm(Y1Z2, p1.y, p2.z, p);
    gcry_mpi_mulm(X1Z2, p1.x, p2.z, p);
    gcry_mpi_mulm(Z1Z2, p1.z, p2.z, p);

    gcry_mpi_mulm(t[0], p2.y, p1.z, p);
    gcry_mpi_t u = gcry_mpi_new(0);
    gcry_mpi_t uu = gcry_mpi_new(0);
    gcry_mpi_subm(u, t[0], Y1Z2, p);
    gcry_mpi_mulm(uu, u, u, p);
    gcry_mpi_mulm(t[1], p2.x, p1.z, p);
    gcry_mpi_t v = gcry_mpi_new(0);
    gcry_mpi_t vv = gcry_mpi_new(0);
    gcry_mpi_t vvv = gcry_mpi_new(0);
    gcry_mpi_subm(v, X1Z2, t[1], p);
    gcry_mpi_mulm(vv, v, v, p);
    gcry_mpi_mulm(vvv, v, vv, p);
    gcry_mpi_t R = gcry_mpi_new(0);
    gcry_mpi_mulm(R, vv, X1Z2, p);
    gcry_mpi_mul_ui(t[2], R, 2);
    gcry_mpi_mulm(t[3], uu, Z1Z2, p);
    gcry_mpi_subm(t[4], t[3], vvv, p);
    gcry_mpi_t A = gcry_mpi_new(0);
    gcry_mpi_subm(A, t[4], t[2], p);
    gcry_mpi_mulm(p3.z, v, A, p);
    gcry_mpi_subm(t[5], R, A, p);
    gcry_mpi_mulm(t[6], vvv, Y1Z2, p);
    gcry_mpi_mulm(t[7], u, t[5], p);
    gcry_mpi_subm(p3.y, t[7], t[6], p);
    gcry_mpi_mulm(p3.z, vvv, Z1Z2, p);

    for(int i = 0; i < 8; ++i)
        gcry_mpi_release(t[i]);
    gcry_mpi_release(Y1Z2);
    gcry_mpi_release(X1Z2);
    gcry_mpi_release(Z1Z2);
    gcry_mpi_release(u);
    gcry_mpi_release(uu);
    gcry_mpi_release(v);
    gcry_mpi_release(vv);
    gcry_mpi_release(vvv);
    gcry_mpi_release(R);
    gcry_mpi_release(A);*/

    /*gcry_mpi_t A = gcry_mpi_new(0);
    gcry_mpi_t B = gcry_mpi_new(0);
    gcry_mpi_t C = gcry_mpi_new(0);
    gcry_mpi_t D = gcry_mpi_new(0);
    gcry_mpi_t E = gcry_mpi_new(0);
    gcry_mpi_t temp1 = gcry_mpi_new(0);
    gcry_mpi_t temp2 = gcry_mpi_new(0);

    // Вычисляем A
    gcry_mpi_mulm(temp1, p2.y, p1.z, p);
    gcry_mpi_mulm(temp2, p1.y, p2.z, p);
    gcry_mpi_subm(A, temp1, temp2, p);

    // Вычисляем B
    gcry_mpi_mulm(temp1, p2.x, p1.z, p);
    gcry_mpi_mulm(temp2, p1.x, p2.z, p);
    gcry_mpi_subm(B, temp1, temp2, p);

    // Вычисляем C
    gcry_mpi_mulm(temp1, p2.x, p1.z, p);
    gcry_mpi_mulm(temp2, p1.x, p2.z, p);
    gcry_mpi_addm(A, temp1, temp2, p);

    // Вычисляем D
    gcry_mpi_mulm(temp1, p1.x, p2.z, p);
    gcry_mpi_mul_ui(temp1, temp1, 2);
    gcry_mpi_mulm(temp2, p2.x, p1.z, p);
    gcry_mpi_addm(A, temp1, temp2, p);

    // Вычисляем E
    gcry_mpi_mulm(E, p1.z, p2.z, p);

    // Вычисляем X3 = B(A^2E - B^2C)
    gcry_mpi_t exp = gcry_mpi_new(0);
    gcry_mpi_set_ui(exp, 2);
    gcry_mpi_powm(temp1, A, exp, p);
    gcry_mpi_mulm(temp1, temp1, E, p);
    gcry_mpi_powm(temp2, B, exp, p);
    gcry_mpi_mulm(temp2, temp2, C, p);
    gcry_mpi_subm(p3.x, temp1, temp2, p);
    gcry_mpi_mulm(p3.x, p3.x, B, p);

    // Вычисляем Y3 = A(B^2D - A^2E) - Y1Z2B^3
    gcry_mpi_t temp3 = gcry_mpi_new(0);
    gcry_mpi_powm(temp1, B, exp, p);
    gcry_mpi_mulm(temp1, temp1, D, p);
    gcry_mpi_powm(temp2, A, exp, p);
    gcry_mpi_mulm(temp2, temp2, E, p);
    gcry_mpi_subm(temp1, temp1, temp2, p);
    gcry_mpi_mulm(temp1, temp1, A, p);
    gcry_mpi_set_ui(exp, 3);
    gcry_mpi_powm(temp3, B, exp, p);
    gcry_mpi_mulm(temp3, temp3, p1.y, p);
    gcry_mpi_mulm(temp3, temp3, p2.z, p);
    gcry_mpi_subm(p3.y, temp1, temp3, p);

    // Вычисляем Z3 = B^3E
    gcry_mpi_powm(p3.z, B, exp, p);
    gcry_mpi_mulm(p3.z, p3.z, E, p);*/

    cout << "P+Q:\n";
    p3.print();
    return p3;
}

Point EllipticCurve::comp_mult_point(const Point &point)
{
    gcry_mpi_set_ui(Q.x, 0);
    gcry_mpi_set_ui(Q.y, 1);
    gcry_mpi_set_ui(Q.z, 0);
}

