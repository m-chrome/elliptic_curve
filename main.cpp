#include <ec.hpp>
#include <iostream>
#include <gcrypt.h>

using namespace std;

int main()
{
    EC test;
    test.check_p_point();
    return 0;
}

