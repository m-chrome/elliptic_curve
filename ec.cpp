#include <ec.hpp>
#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <gcrypt.h>
#include <gmp.h>

#define MAX_MPI_BUF 200

using namespace std;

// Обработка большого числа
unsigned char* process_mpi(const char *raw_number, gcry_mpi_t &mpi)
{
    mpi = gcry_mpi_new(MAX_MPI_BUF);
    char str_number[MAX_MPI_BUF] = {0};
    strcpy(str_number, raw_number);
    mpi = gcry_mpi_new(MAX_MPI_BUF);
    gcry_mpi_scan(&mpi, GCRYMPI_FMT_USG, str_number, sizeof(str_number), NULL);
    unsigned char *check_buf = (unsigned char *)malloc(MAX_MPI_BUF);
    if (gcry_mpi_print(GCRYMPI_FMT_USG, check_buf, MAX_MPI_BUF, NULL, mpi) == 0)
    {
        free(check_buf);
        return check_buf;
    }
    else
    {
        free(check_buf);
        return NULL;
    }
}

void show_mpi(gcry_mpi_t mpi)
{
    unsigned char *check_buf = (unsigned char *)malloc(MAX_MPI_BUF);
    memset(check_buf, '\0', MAX_MPI_BUF);
    gcry_mpi_print(GCRYMPI_FMT_USG, check_buf, MAX_MPI_BUF, NULL, mpi);
    printf("%s\n", check_buf);
}

EC::EC()
{
    cout << "Эллиптическая кривая" << endl;
    cout << "Форма Вейерштрасса: y^2 = x^3 + ax + b (mod p)" << endl;
    cout << "Параметры кривой:\n" << endl;

    // Заполнение параметрами
    printf("p = %s\n", process_mpi(ec_param.p, p));
    printf("a = %s\n", process_mpi(ec_param.a, a));
    printf("b = %s\n", process_mpi(ec_param.b, b));
    printf("q = %s\n", process_mpi(ec_param.q, q));

    // Точка из ГОСТ
    printf("\nТочка из ГОСТ:\n");
    gcry_mpi_t x, y;
    printf("x0 = %s\n", process_mpi(point_param.x, x));
    printf("y0 = %s\n\n", process_mpi(point_param.y, y));
    P0 = gcry_mpi_point_snatch_set(P0, x, y, NULL);
}

EC::~EC()
{
    cout << "Очистка памяти." << endl;
}

void EC::check_p_point()
{
    // Проверка точки на принадлежность к кривой
    gcry_mpi_t x0, y0, A, B, C, exp, X, Y, div;
    A = gcry_mpi_new(0);
    B = gcry_mpi_new(0);
    C = gcry_mpi_new(0);
    X = gcry_mpi_new(0);
    Y = gcry_mpi_new(0);
    x0 = gcry_mpi_new(0);
    y0 = gcry_mpi_new(0);
    exp = gcry_mpi_new(0);
    div = gcry_mpi_new(0);
    gcry_mpi_point_get(x0, y0, NULL, P0);

    // X: x^3 + ax + b (p)
    gcry_mpi_scan(&exp, GCRYMPI_FMT_USG, "3", 1, NULL);
    gcry_mpi_powm(A, x0, exp, p);
    gcry_mpi_mul(B, a, x0);
    gcry_mpi_add(C, B, b);
    gcry_mpi_add(X, A, C);
    gcry_mpi_mod(X, X, p);

    // Критерий Эйлера X^((p-1)/2) = 1 (p)
    gcry_mpi_sub_ui(exp, p, 1);
    gcry_mpi_scan(&div, GCRYMPI_FMT_USG, "2", 1, NULL);
    gcry_mpi_div(exp, NULL, exp, div, 0);
    gcry_mpi_powm(A, X, exp, p);

    if (gcry_mpi_cmp_ui(A,1) == 0)
        printf("Критерий Эйлера выполнен.\n");

    // Находим Y = X^((p+1)/4) (p)
    gcry_mpi_add_ui(exp, p, 1);
    gcry_mpi_scan(&div, GCRYMPI_FMT_USG, "4", 1, NULL);
    gcry_mpi_div(exp, NULL, exp, div, 0);
    gcry_mpi_powm(Y, X, exp, p);

    if (gcry_mpi_cmp(y0,Y) == 0)
        printf("Точка принадлежит кривой.\n");
    else
        printf("Точка не принадлежит кривой.\n");

    gcry_mpi_release(A);
    gcry_mpi_release(B);
    gcry_mpi_release(C);
    gcry_mpi_release(exp);
}

/*void EC::build_point()
{
    srand(time(NULL));
    cout << "Построение случайной точки:\n";
    mpz_t p, a, b, fx0, x0;
    mpz_init_set_str(p, ec_param.p, 10);
    mpz_init_set_str(a, ec_param.a, 10);
    mpz_init_set_str(b, ec_param.b, 10);

    unsigned long int x0_ui = rand() % RAND_MAX;
    while (true)
    {
        cout << "Выбрано случайное число x0: " << x0_ui << endl;
        mpz_init_set_ui(x0, x0_ui);
        if (residue(x0_ui) == 1)
        {
            cout << "f(x0) является квадратичным вычетом y по модулю p." << endl;
            mpz_set(fx0, count_fx0(x0_ui));
            break;
        }
        else
        {
            cout << "f(x0) является невычетом по модулю p. Генерирую x0 заново." << endl;
            x0_ui = rand() % RAND_MAX;
        }
    }


    mpz_t y0, exp;
    mpz_inits(y0, exp, NULL);

    // Рассчёт степени (p+1)/4
    mpz_add_ui(exp, p, 1);
    mpz_div_ui(exp, exp, 4);

    // Рассчитать y0
    char *str1={0}, *str2={0};
    mpz_powm(y0, fx0, exp, p);
    cout << "x0: " << mpz_get_str(str1, 10, x0) << endl;
    cout << "y0: " << mpz_get_str(str2, 10, y0) << endl;
    mpz_clears(p,a,b,exp,x0,fx0);
}

// Посчитать f(x0)
mpz_t& count_fx0(unsigned long x0_uint)
{
    mpz_t a, b, x0, fx0, first, second;
    mpz_inits(fx0, first, second, NULL);
    mpz_init_set_ui(x0, x0_uint);
    mpz_init_set_str(a, ec_param.a, 10);
    mpz_init_set_str(b, ec_param.b, 10);

    // x^3
    mpz_pow_ui(first, x0, 3);

    // a*x + b
    mpz_mul(second, a, x0);
    mpz_add(second, second, b);

    // x^3 + a*x + b
    mpz_add(fx0, first, second);
    return fx0;
}

// Найти выражение y^2 = f(x0) (mod p) и символ Лежандра
int residue(unsigned long x0_uint)
{
    mpz_t fx0, p;
    mpz_init_set(fx0, count_fx0(x0_uint));
    mpz_init_set_str(p, ec_param.p, 10);

    // Квадратичный вычет
    return mpz_legendre(fx0, p);
}*/
