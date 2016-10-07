#ifndef __INTERPOLATE_PATTERN_H__
#define __INTERPOLATE_PATTERN_H__

#include "cfloat.h"
#include "pattern.h"
#include "vec3.h"

void interpolate_pattern(const double phi, const double theta, const pattern *patt, cfloat &h_theta, cfloat &h_phi, vec3 &theta_v, vec3 &phi_v);

#endif
