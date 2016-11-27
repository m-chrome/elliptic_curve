#include <ec.hpp>
#include <iostream>
#include <gcrypt.h>

using namespace std;

int main()
{
    EllipticCurve curve;
    curve.build_point(0);
    return 0;
}

