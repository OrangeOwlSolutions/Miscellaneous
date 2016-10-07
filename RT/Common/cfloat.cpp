#include "cfloat.h"

#include <iostream>

std::ostream & operator<<(std::ostream &os, const cfloat &cf){
    return os << REAL(cf) << " " << IMAG(cf) ;
}

