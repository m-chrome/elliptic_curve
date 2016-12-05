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

    // Генерация k
    k = gcry_mpi_new(0);
    if (gcry_mpi_scan(&k, GCRYMPI_FMT_HEX, "3", 0, NULL) != 0)
        cout << "Ошибка записи k.\n";
    cout << "k = ";
    show_mpi(k);

    // Генерация l
    l = gcry_mpi_new(0);
    if (gcry_mpi_scan(&l, GCRYMPI_FMT_HEX, "4", 0, NULL) != 0)
        cout << "Ошибка записи l.\n";
    cout << "l = ";
    show_mpi(l);

    // Генерация m
    m = gcry_mpi_new(0);
    gcry_mpi_addm(m, k, l, q);
    cout << "m = ";
    show_mpi(m);
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
    gcry_mpi_t exp = gcry_mpi_new(0);
    gcry_mpi_t div = gcry_mpi_new(0);
    gcry_mpi_add_ui(exp, p, 1);
    gcry_mpi_set_ui(div, 4);
    gcry_mpi_div(exp, NULL, exp, div, 0);
    gcry_mpi_powm(y0, fx0, exp, p);
    gcry_mpi_release(exp);
    gcry_mpi_release(div);
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
        break;
    }
    case 1:
    {
        cout << "Будет сгенерирована случайная точка P(x,y,z).\n";
        gcry_mpi_randomize(P.x, 254, GCRY_WEAK_RANDOM);
        while (euler_criteria(comp_fx0(P.x))!=0)
            gcry_mpi_randomize(P.x, 254, GCRY_WEAK_RANDOM);
        P.y = comp_y0(P.x);
        gcry_mpi_set_ui(P.z, 1);
        P.print();
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
    if (check_projective_point_belongs(P) == 0)
    {
        cout << "Точка P(x,y,z) принадлежит кривой.\n";
        cout << endl;
        return 0;
    }
    else
    {
        cout << "Точка P(x,y,z) не принадлежит кривой.\n";
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

void EllipticCurve::doubling_point(Point &dupPoint, const Point &point)
{
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
    gcry_mpi_mulm(dupPoint.x, h, s, p);
    gcry_mpi_subm(t[8], B, h, p);
    gcry_mpi_mul_ui(t[9], rr, 2);
    gcry_mpi_mulm(t[10], t[8], w, p);
    gcry_mpi_subm(dupPoint.y, t[10], t[9], p);
    gcry_mpi_set(dupPoint.z, sss);
}

void EllipticCurve::add_points(Point &p3, const Point &p1, const Point &p2)
{
    // http://hyperelliptic.org/EFD/g1p/auto-code/shortw/projective/addition/add-1998-cmo-2.op3
    gcry_mpi_t t[8];
    for(int i = 0; i < 8; ++i)
        t[i] = gcry_mpi_new(0);
    gcry_mpi_t Y1Z2 = gcry_mpi_new(0);
    gcry_mpi_t X1Z2 = gcry_mpi_new(0);
    gcry_mpi_t Z1Z2 = gcry_mpi_new(0);
    gcry_mpi_mul(Y1Z2, p1.y, p2.z);
    gcry_mpi_mul(X1Z2, p1.x, p2.z);
    gcry_mpi_mul(Z1Z2, p1.z, p2.z);
    gcry_mpi_mul(t[0], p2.y, p1.z);

    gcry_mpi_t u = gcry_mpi_new(0);
    gcry_mpi_sub(u, t[0], Y1Z2);

    gcry_mpi_t uu = gcry_mpi_new(0);
    gcry_mpi_mul(uu, u, u);
    gcry_mpi_mul(t[1], p2.x, p1.z);

    gcry_mpi_t v = gcry_mpi_new(0);
    gcry_mpi_t vv = gcry_mpi_new(0);
    gcry_mpi_t vvv = gcry_mpi_new(0);
    gcry_mpi_sub(v, t[1], X1Z2);
    gcry_mpi_mul(vv, v, v);
    gcry_mpi_mul(vvv, v, vv);

    gcry_mpi_t R = gcry_mpi_new(0);
    gcry_mpi_mul(R, vv, X1Z2);
    gcry_mpi_mul_ui(t[2], R, 2);
    gcry_mpi_mul(t[3], uu, Z1Z2);
    gcry_mpi_sub(t[4], t[3], vvv);

    gcry_mpi_t A = gcry_mpi_new(0);
    gcry_mpi_sub(A, t[4], t[2]);

    gcry_mpi_mul(p3.x, v, A);

    gcry_mpi_sub(t[5], R, A);
    gcry_mpi_mul(t[6], vvv, Y1Z2);

    gcry_mpi_mul(t[7], u, t[5]);

    gcry_mpi_sub(p3.y, t[7], t[6]);

    gcry_mpi_mul(p3.z, vvv, Z1Z2);

    gcry_mpi_mod(p3.x, p3.x, p);
    gcry_mpi_mod(p3.y, p3.y, p);
    gcry_mpi_mod(p3.z, p3.z, p);

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
    gcry_mpi_release(A);
}

void EllipticCurve::comp_mult_point(Point &kp, const Point &p, const gcry_mpi_t m)
{
    gcry_mpi_set_ui(kp.x, 0);
    gcry_mpi_set_ui(kp.y, 0);
    gcry_mpi_set_ui(kp.z, 0);

    Point tempP;
    gcry_mpi_set(tempP.x, p.x);
    gcry_mpi_set(tempP.y, p.y);
    gcry_mpi_set(tempP.z, p.z);

    int last = 0;
    for(int i = 0; i < (int)gcry_mpi_get_nbits(m); i++)
    {
        if (gcry_mpi_test_bit(m, i))
        {
            for(int j = last; j<i; ++j)
            {
                doubling_point(tempP, tempP);
            }
            last = i;
            if (!(gcry_mpi_cmp_ui(kp.x, 0) && gcry_mpi_cmp_ui(kp.y, 0) && gcry_mpi_cmp_ui(kp.z, 0)))
            {
                kp.x = gcry_mpi_copy(tempP.x);
                kp.y = gcry_mpi_copy(tempP.y);
                kp.z = gcry_mpi_copy(tempP.z);
            }
            else
                add_points(kp, kp, tempP);
        }
    }
}

bool EllipticCurve::extra_check()
{
    Point Q, Q1, Q2;
    comp_mult_point(Q1, P, k);
    comp_mult_point(Q2, P, l);

    /*cout << "Q1:\n";
    Q1.print();
    cout << "Q2:\n";
    Q2.print();*/

    add_points(Q, Q1, Q2);
    /*cout << "Q:\n";
    Q.print();*/

    Point R;
    comp_mult_point(R, P, m);
    /*cout << "R:\n";
    R.print();*/
    gcry_mpi_release(m);

    if (!(gcry_mpi_cmp(Q.x, R.x) && gcry_mpi_cmp(Q.y, R.y) && gcry_mpi_cmp(Q.z, R.z)))
        return true;
    else
        return false;
}
