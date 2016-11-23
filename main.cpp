#include <ec.hpp>
#include <iostream>
#include <gcrypt.h>

using namespace std;

int main()
{
    EC test;
    test.build_point();
    return 0;
}

