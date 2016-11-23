#include <ec.hpp>
#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <gcrypt.h>
#include <gmp.h>

#define MAX_MPI_BUF 200

using namespace std;

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
}


// Обработка большого числа
unsigned char* process_mpi(const char *raw_number, gcry_mpi_t mpi)
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

EC::EC()
{
    // Параметры из ГОСТ-34.10-2012, пример 1
    cout << "Эллиптическая кривая" << endl;
    cout << "Форма Вейерштрасса: y^2 = x^3 + ax + b (mod p)" << endl;
    cout << "Параметры кривой:\n" << endl;

    // Заполнение параметрами
    printf("p = %s\n", process_mpi(ec_param.p, p));
    printf("a = %s\n", process_mpi(ec_param.a, a));
    printf("b = %s\n", process_mpi(ec_param.b, b));
    printf("q = %s\n", process_mpi(ec_param.q, q));
}

EC::~EC()
{
    cout << "Очистка памяти." << endl;
    gcry_mpi_release(p);
    gcry_mpi_release(a);
    gcry_mpi_release(b);
    gcry_mpi_release(q);
    gcry_mpi_release(k);
    gcry_mpi_point_release(P0);
    gcry_mpi_point_release(Q);
}

void EC::build_point()
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
