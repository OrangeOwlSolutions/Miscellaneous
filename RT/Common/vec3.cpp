#include "vec3.h"

std::ostream &operator << (std::ostream &os ,const vec3 & v){
    return os << v.x << " " << v.y <<" " << v.z ; 
}

