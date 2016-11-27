#include <ec.hpp>
#include <iostream>
#include <gcrypt.h>

using namespace std;

int main()
{
    EC test;
    if (test.check_point(test.P0))
        printf("Точка Р принадлежит кривой.\n");
    else
        printf("Точка Р не принадлежит кривой.\n");
    return 0;
}

