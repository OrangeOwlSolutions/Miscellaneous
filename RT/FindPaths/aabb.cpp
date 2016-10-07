#include "aabb.h"

/*******************/
/* PRINTS THE AABB */
/*******************/
std::ostream &operator << (std::ostream &os, const aabb & bb) {
    return os << bb.min_vec << "  " << bb.max_vec;
}
