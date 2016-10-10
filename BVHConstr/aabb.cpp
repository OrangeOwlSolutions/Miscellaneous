#include "aabb.h"
#include <cmath>

#define EPS_R 2.22044604925031308e-16;

#define maximum(a, b) (((a) > (b)) ? (a) : (b))
#define minimum(a, b) (((a) < (b)) ? (a) : (b))

/*************************/
/* ALMOST EQUAL FUNCTION */
/*************************/
static bool almostEqual(double x, double y) {
#ifdef __linux__
    double max_x_y_one = fmax(fmax(1.0, fabs(x)), fabs(y)) ;
#elif _WIN32 || _WIN64
    double max_x_y_one = maximum(maximum(1.0, fabs(x)), fabs(y));
#endif
    return fabs(x - y) <= max_x_y_one * EPS_R ;
}

/************************/
/* DISPLAY BOUNDING BOX */
/************************/
std::ostream & operator << (std::ostream &os, const aabb & bb) { 
	return os << "Lower left point = " << bb.lowerLeftPoint  << " Upper right point = " << bb.upperRightPoint;
}

/****************************/
/* AABB DEFAULT CONSTRUCTOR */
/****************************/
aabb::aabb() : empty(true), lowerLeftPoint(make_vec3d(0., 0., 0.)), upperRightPoint(make_vec3d(0., 0., 0.)) {};

/********************/
/* AABB CONSTRUCTOR */
/********************/
aabb::aabb(const vec3d &lowerLeftPoint_, const vec3d &upperLeftPoint_) : empty(false), lowerLeftPoint(lowerLeftPoint_), upperRightPoint(upperLeftPoint_) {}

/************************/
/* AABB GROWTH FUNCTION */
/************************/
void aabb::grow(const vec3d &v) {
	if (!empty) {
		// --- If the AABB is not empty, then enlarges the AABB with v and the lower left and upper right points
        lowerLeftPoint	= cwisemin(lowerLeftPoint, v);
        upperRightPoint = cwisemax(upperRightPoint, v);
    } else {
        // --- If the AABB is empty, then the largest AABB "coincides" with v
		lowerLeftPoint	= v;
        upperRightPoint = v;
        empty			= false;
    }     
}

void aabb::grow(const aabb &bb) { grow(bb.lowerLeftPoint); grow(bb.upperRightPoint); }

/*******************/
/* GET AABB CENTER */
/*******************/
vec3d aabb::get_center() const { return (upperRightPoint + lowerLeftPoint) * 0.5; }

/**********************/
/* GET AABB HALF-SIZE */
/**********************/
vec3d aabb::get_half_size() const { return (upperRightPoint - lowerLeftPoint) * 0.5; }

/*************************/
/* GET AABB SURFACE AREA */
/*************************/
double aabb::get_surface_area() const {
    vec3d side = upperRightPoint - lowerLeftPoint;
    return 2.0 * (side.x * side.y + side.x * side.z + side.y * side.z);
}
