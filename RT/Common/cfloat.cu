#include "cfloat.h"

#ifdef __linux__
HOST_DEVICE_INLINE
#endif
float &REAL(cfloat &a){
    return a.x;
}

