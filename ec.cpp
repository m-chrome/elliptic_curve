#include <ec.hpp>
#include <iostream>
#include <string>
#include <cstdlib>

#define MAX_MPI_BUF 200

using namespace std;

// Обработка большого числа
unsigned char* process_mpi(const char *raw_number, gcry_mpi_t mpi)
{
    mpi = gcry_mpi_new(MAX_MPI_BUF);
    char str_number[MAX_MPI_BUF] = {0};
    strcpy(str_number, raw_number);
    mpi = gcry_mpi_new(MAX_MPI_BUF);
    gcry_mpi_scan(&mpi, GCRYMPI_FMT_USG, str_number, sizeof(str_number), NULL);
    //unsigned char check_buf[MAX_MPI_BUF]={0};
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
    printf("p = %s\n", process_mpi("57896044618658097711785492504343953926634992332820282019728792003956564821041", p));
    printf("a = %s\n", process_mpi("7", a));
    printf("b = %s\n", process_mpi("43308876546767276905765904595650931995942111794451039583252968842033849580414", b));
    printf("q = %s\n", process_mpi("57896044618658097711785492504343953927082934583725450622380973592137631069619", q));

    // Создание S-expression
    //gcry_sexp_new(&weierstrass, , 1024, 0);
    //gcry_sexp_build(&weierstrass, NULL,"(ecc (p %m) (a %m) (b %m) (n %m))", p, a, b, q);
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
    gcry_sexp_release(weierstrass);
}
