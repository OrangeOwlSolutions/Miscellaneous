//#ifndef __TRANSPORT_RAY_FIELD__
//#define __TRANSPORT_RAY_FIELD__
//
//#include <vector>
//
//#include "vec3.h"
//#include "cfloat.h"
//#include "RXTX.h"
//
//extern "C" cfloat transport_ray_field(std::vector<vec3> positions, TX &tx, RX &rx, double lambda);
//
//#endif
#ifndef __TRANSPORT_RAY_FIELD__
#define __TRANSPORT_RAY_FIELD__

#include <vector>

#include "vec3.h"
#include "cfloat.h"
#include "RXTX.h"
#include "load_nbin.h"
#include "ray_path.h"
#include "edge.h"

//extern "C" cfloat transport_ray_field(std::vector<vec3> positions, TX &tx, RX &rx, double lambda);
extern "C" cfloat transport_ray_field(const path &, TX &, RX &, const nvmesh &, const std::vector<edge> &, double lambda);

#endif
