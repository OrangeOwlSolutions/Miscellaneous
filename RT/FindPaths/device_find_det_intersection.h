#ifndef __DEVICE_FIND_DET_INTERSECTION_H__
#define __DEVICE_FIND_DET_INTERSECTION_H__

#include "ray_list.h"
#include "geointer.h"
#include "RXTX.h"


void device_find_det_inters(
    device_ray_list::R0_R1_iterator begin,
    device_ray_list::R0_R1_iterator end,
    thrust::device_vector<float>::iterator dinter_begin ,
    thrust::device_vector<geointer>::iterator geointer_begin,
    const RX &rx
);



#endif
