#ifndef __AABB_H__
#define __AABB_H__

#include "vec3.h"
#include "vec4.h"

/***************/
/* AABB STRUCT */
/***************/
struct aabb {
    
	bool		empty;
    vec3d		lowerLeftPoint;
    vec3d		upperRightPoint;
    
    aabb();
    aabb(const vec3d &, const vec3d &);
    void    grow(const vec3d &);
    void    grow(const aabb &);
    vec3d   get_center() const;
    vec3d   get_half_size() const ;
    double  get_surface_area() const;
};

std::ostream & operator << (std::ostream &, const aabb &);

#endif
